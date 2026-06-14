"""
scATrans tl module.

Primary functions:
- `active_score`: composite active transcription scoring from velocity (spliced/unspliced)
  + differential expression (supports multiple backends including Memento).
- `differential_expression`: standalone DE (supports scanpy methods, PyDESeq2 pseudobulk,
  mixed models, and Memento as a first-class Cell 2024 method-of-moments backend).
  Useful when you have no velocity layers and only want DE + downstream enrichment/plotting.

Downstream tools (`filter_active_genes`, `run_enrichment`, `scat.pl.*`) work on the
results DataFrames from either function.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse  # for type hints in signatures (e.g. spmatrix)
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests

# qc is imported lazily inside active_score to keep startup light, but exposed at package level
from . import qc as _qc  # for unspliced_global integration
from ._bias import fit_huber_bias_correction
from ._de import _run_de_wrapper
from ._permutation import _single_permutation_task
from ._utils import (
    _get_exponential_scale_lambda,
    _is_integer_counts_like,
    _normalize_velocity_layers_by_size_factor,
    _pseudobulk_with_layers,
    _soft_scale,
    comb,  # for small-n permutation space calculation
)
from ._velocity import _compute_moments_velocity_delta, _compute_velocity_delta

try:
    from . import _version

    VERSION = _version.version
except (ImportError, AttributeError):
    VERSION = "0.7.0.dev0"

logger = logging.getLogger(__name__)


def _validate_de_common_options(
    *,
    de_preprocess: str,
    pseudobulk_de_backend: str,
    n_jobs: int,
    use_permutation: bool,
    n_perm: int,
    use_mixed_model: bool,
    mixed_model_pval: str,
    use_memento_de: bool,
    memento_capture_rate: float,
    memento_num_boot: int,
    min_cells: int | None = None,
    min_counts: int | None = None,
) -> None:
    """Shared early validation for DE-related options used by both
    active_score() and differential_expression()."""
    if de_preprocess not in {"auto", "normalize_log1p", "none"}:
        raise ValueError("de_preprocess must be one of {'auto', 'normalize_log1p', 'none'}.")

    if pseudobulk_de_backend not in {"pydeseq2", "scanpy"}:
        raise ValueError("pseudobulk_de_backend must be 'pydeseq2' or 'scanpy'.")

    if not isinstance(n_jobs, int):
        raise ValueError("n_jobs must be an integer.")

    if use_permutation and (not isinstance(n_perm, int) or n_perm < 1):
        raise ValueError("n_perm must be a positive integer when use_permutation=True.")

    if not isinstance(use_mixed_model, bool):
        raise ValueError("use_mixed_model must be boolean.")

    if mixed_model_pval not in ("wald", "lrt"):
        raise ValueError("mixed_model_pval must be 'wald' or 'lrt'.")

    if not isinstance(use_memento_de, bool):
        raise ValueError("use_memento_de must be boolean.")

    if not (0 < memento_capture_rate < 1):
        raise ValueError(
            "memento_capture_rate must be in (0, 1). Typical values: ~0.07 for 10x v1, ~0.15 for v2."
        )

    if memento_num_boot < 100:
        raise ValueError(
            "memento_num_boot should be reasonably large (>=100) for stable estimates."
        )

    if min_cells is not None and min_cells < 1:
        raise ValueError("min_cells must be >= 1.")

    if min_counts is not None and min_counts < 0:
        raise ValueError("min_counts must be non-negative.")

    # Extra type guards to prevent bool/int confusion and weird values
    if not isinstance(use_permutation, bool):
        raise ValueError("use_permutation must be boolean.")

    if not isinstance(use_mixed_model, bool):
        raise ValueError("use_mixed_model must be boolean.")

    if not isinstance(use_memento_de, bool):
        raise ValueError("use_memento_de must be boolean.")

    if not isinstance(memento_num_boot, int) or isinstance(memento_num_boot, bool):
        raise ValueError("memento_num_boot must be an integer.")

    if not isinstance(n_perm, int) or isinstance(n_perm, bool):
        raise ValueError("n_perm must be an integer.")

    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool):
        raise ValueError("n_jobs must be an integer.")

    if min_cells is not None and (not isinstance(min_cells, int) or isinstance(min_cells, bool)):
        raise ValueError("min_cells must be an integer >= 1.")

    if min_counts is not None and (not isinstance(min_counts, int) or isinstance(min_counts, bool)):
        raise ValueError("min_counts must be a non-negative integer.")


def _coerce_memento_de_preprocess(use_memento_de: bool, de_preprocess: str) -> str:
    """If use_memento_de, force de_preprocess to 'none' (Memento needs raw counts)."""
    if use_memento_de and de_preprocess != "none":
        logger.info(
            "use_memento_de=True: forcing de_preprocess='none' "
            "(Memento method-of-moments works on raw counts)."
        )
        return "none"
    return de_preprocess


def active_score(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "GA",
    reference_group: str = "Ctrl",
    subset_col: str | None = None,
    subset_values: str | list[str] | tuple[str, ...] | None = None,
    weight_fc: float = 1.0,
    weight_unspliced: float = 1.0,
    weight_pval: float = 1.0,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.5,
    active_fdr_cutoff: float = 0.05,
    de_method: str = "t-test_overestim_var",  # freely switchable basic option, e.g. "wilcoxon"
    pseudobulk_de_backend: str = "pydeseq2",  # "pydeseq2" or "scanpy" when use_pseudobulk=True
    use_permutation: bool = False,
    perm_de_backend: str = "fast",
    n_perm: int = 100,
    n_jobs: int = -1,
    de_preprocess: str = "auto",
    gene_type_filter: str | None = None,
    use_pseudobulk: bool = False,
    sample_col: str | None = None,
    min_cells: int = 10,
    min_counts: int = 1000,
    pb_x_layer: str = "spliced",
    pb_use_total_for_x: bool = True,
    min_total_counts: int = 50,
    random_seed: int = 42,
    show_plot: bool = False,
    auto_adjust_n_perm: bool = True,
    prior_weight: float = 5.0,
    strict_pydeseq2_counts: bool = True,
    # Dual-track
    mode: str = "heuristic",
    advanced_fallback: bool = True,
    advanced_n_neighbors: int = 30,
    advanced_n_pcs: int = 30,
    advanced_use_precomputed: bool = False,
    allow_advanced_pseudobulk: bool = False,
    advanced_recompute_neighbors: bool = True,
    spliced_layer: str = "spliced",
    unspliced_layer: str = "unspliced",
    # Mixed model (dreamlet/variancePartition-style LMM + delta variance via statsmodels)
    use_mixed_model: bool = False,
    use_delta_variance_pval: bool = False,
    delta_var_pval_cutoff: float = 0.05,
    mixed_model_pval: str = "wald",  # "wald" or "lrt" - which p-value to use for the DE part of active_score when use_mixed_model=True
    # Memento (independent cell-level method-of-moments backend, Cell 2024)
    # Third parallel DE path (alongside scanpy-style and pseudobulk). Only used for the main DE statistics.
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    # Advanced control for permutation (option B): default keeps fast/cheap DE in perms for speed.
    # Set to True only if you really want the null to be generated with the exact same (expensive) Memento backend.
    perm_use_memento_de: bool = False,
    # Bias correction control (default on for basic pipeline cleanliness; opt-out for exploration)
    bias_correction: str = "huber_length_intron",
    # Opt-in transparency for per-gene reference gamma (keeps default output clean)
    show_effective_gamma: bool = False,
    # Advanced convenience for users primarily interested in nascent RNA excess
    prioritize_velocity: bool = False,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Identify genes with condition-wise differences in unspliced (nascent) RNA abundance
    relative to a reference group, using a composite score that also incorporates
    differential expression statistics.

    The function computes:
    - logFC and p_adj between target and reference (via scanpy or PyDESeq2).
    - A velocity delta = U_target − (gamma_ref × S_target), where gamma_ref is a
      shrunk U/S ratio estimated in the reference group.
    - (by default) A Huber regression correction of the delta on log(gene length) and
      log(intron number); the residuals become velocity_residual.
    - A soft-scaled, weighted combination of the three signals, scaled to 0–100.

    Several extensions are available as explicit options (see the README section
    "Optional advanced features"):
    - show_effective_gamma
    - bias_correction="none"
    - use_mixed_model
    - use_permutation
    - prioritize_velocity (convenience for analyses focused on the unspliced excess term)

    Diagnostics (including global unspliced fraction and bias fit details) are stored
    under adata.uns["scatrans"]["diagnostics"]. The full ranked table (all_results)
    is the main output; the built-in significant list is produced by a strict
    conjunction of thresholds and is often small or empty.

    A separate function diagnose_design is available to summarize the experimental
    design and surface relevant warnings before analysis.

    Full usage, recommended workflow, and result interpretation are documented in
    the package README.
    """
    # ==================== EARLY VALIDATION (kept identical) ====================
    if mode not in {"heuristic", "advanced"}:
        raise ValueError("mode must be either 'heuristic' or 'advanced'")

    if not isinstance(advanced_fallback, bool):
        raise ValueError("advanced_fallback must be boolean.")
    if not isinstance(advanced_n_neighbors, int) or advanced_n_neighbors < 2:
        raise ValueError("advanced_n_neighbors must be integer >= 2.")
    if not isinstance(advanced_n_pcs, int) or advanced_n_pcs < 2:
        raise ValueError("advanced_n_pcs must be integer >= 2.")
    if not isinstance(advanced_use_precomputed, bool):
        raise ValueError("advanced_use_precomputed must be boolean.")
    if not isinstance(allow_advanced_pseudobulk, bool):
        raise ValueError("allow_advanced_pseudobulk must be boolean.")
    if not isinstance(advanced_recompute_neighbors, bool):
        raise ValueError("advanced_recompute_neighbors must be boolean.")

    if not isinstance(use_delta_variance_pval, bool):
        raise ValueError("use_delta_variance_pval must be boolean.")
    if not (0 < delta_var_pval_cutoff < 1):
        raise ValueError("delta_var_pval_cutoff must be in (0, 1).")
    if not isinstance(perm_use_memento_de, bool):
        raise ValueError("perm_use_memento_de must be boolean.")

    # Memento requires count data; force no log-norm preprocess for the DE leg
    # (do this BEFORE the single validation call)
    de_preprocess = _coerce_memento_de_preprocess(use_memento_de, de_preprocess)

    # Shared validation for DE options (deduplicated with differential_expression)
    # Only one call in the entire function.
    _validate_de_common_options(
        de_preprocess=de_preprocess,
        pseudobulk_de_backend=pseudobulk_de_backend,
        n_jobs=n_jobs,
        use_permutation=use_permutation,
        n_perm=n_perm,
        use_mixed_model=use_mixed_model,
        mixed_model_pval=mixed_model_pval,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        min_cells=min_cells,
        min_counts=min_counts,
    )

    if not isinstance(use_delta_variance_pval, bool):
        raise ValueError("use_delta_variance_pval must be boolean.")
    if not (0 < delta_var_pval_cutoff < 1):
        raise ValueError("delta_var_pval_cutoff must be in (0, 1).")

    if not (0 < active_fdr_cutoff <= 1):
        raise ValueError("active_fdr_cutoff must be in (0, 1].")

    if min_total_counts < 0:
        raise ValueError("min_total_counts must be non-negative.")

    if perm_de_backend not in {"fast", "same"}:
        raise ValueError("perm_de_backend must be 'fast' or 'same'.")

    # Apply prioritize_velocity convenience (only if user left the default equal weights)
    if prioritize_velocity and weight_fc == 1.0 and weight_unspliced == 1.0 and weight_pval == 1.0:
        weight_unspliced = 3.0
        weight_fc = 0.5
        weight_pval = 0.5
        logger.info(
            "prioritize_velocity=True (advanced option): emphasizing the nascent RNA excess / "
            "velocity_residual component in the active_score. "
            "This is intended for users whose primary interest is condition-wise differences "
            "in unspliced abundance after reference correction."
        )

    if mode == "advanced" and use_pseudobulk and not allow_advanced_pseudobulk:
        raise ValueError(
            "mode='advanced' is not supported with use_pseudobulk=True by default. "
            "Set allow_advanced_pseudobulk=True if you really want to try it."
        )

    if mode == "advanced" and use_pseudobulk and allow_advanced_pseudobulk:
        logger.warning(
            "Advanced mode on pseudobulk was explicitly enabled. "
            "This is experimental and may over-smooth sample-level replicates."
        )

    if mode == "advanced":
        logger.info(
            "Advanced mode is experimental and uses scVelo moments for smoothing, "
            "not scVelo's stochastic/dynamical velocity model."
        )

    logger.info("scATrans %s Analysis started. Mode: %s", VERSION, mode)

    # ==================== SUBSET & BASIC VALIDATION ====================
    if subset_col is not None:
        if subset_col not in adata_input.obs.columns:
            raise ValueError(f"subset_col='{subset_col}' not found in adata.obs.columns")
        if subset_values is None:
            raise ValueError("subset_values must be provided when subset_col is specified")
        if isinstance(subset_values, (str, int, float)):
            subset_values_list = [str(subset_values)]
        else:
            subset_values_list = [str(v) for v in subset_values]
        subset_mask = adata_input.obs[subset_col].astype(str).isin(subset_values_list)
        n_before = adata_input.n_obs
        adata_input = adata_input[subset_mask].copy()
        n_after = adata_input.n_obs
        if n_after == 0:
            raise ValueError(f"No cells remain after subsetting {subset_col}")
        logger.info("Subsetted by %s (%d/%d cells remaining)", subset_col, n_after, n_before)

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    target_group = str(target_group)
    reference_group = str(reference_group)

    if target_group == reference_group:
        raise ValueError("target_group and reference_group must be different.")

    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")

    if target_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"target_group '{target_group}' not found.")
    if reference_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"reference_group '{reference_group}' not found.")

    # Automatic design guidance for small-sample or replicate-structured data
    if sample_col or use_pseudobulk:
        from contextlib import suppress

        with suppress(Exception):
            _ = diagnose_design(
                adata_input,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                sample_col=sample_col,
            )  # never let diagnosis break the main analysis

    # ====================== LAYER NAME HANDLING (kb_python support) ======================
    available_layers = list(adata_input.layers.keys())

    if spliced_layer not in available_layers or unspliced_layer not in available_layers:
        if "mature" in available_layers and "nascent" in available_layers:
            logger.warning(
                "Standard 'spliced'/'unspliced' layers not found in adata.layers. "
                "Auto-detected kb_python layers → using 'mature' as spliced (mature mRNA) "
                "and 'nascent' as unspliced (nascent pre-mRNA). "
                "All internal processing will use standard names after remapping. "
                "You can override with spliced_layer= / unspliced_layer= if needed."
            )
            spliced_layer = "mature"
            unspliced_layer = "nascent"
        else:
            raise ValueError(
                f"Required layers not found. "
                f"Expected '{spliced_layer}' + '{unspliced_layer}' (or kb_python 'mature' + 'nascent'). "
                f"Available layers: {available_layers}"
            )

    keep_mask = adata_input.obs[groupby].astype(str).isin([target_group, reference_group])
    adata = adata_input[keep_mask].copy()

    # Perform layer remapping only on the copied adata to avoid mutating the caller's original object.
    if (
        (spliced_layer != "spliced" or unspliced_layer != "unspliced")
        and spliced_layer in adata.layers
        and unspliced_layer in adata.layers
    ):
        adata.layers["spliced"] = adata.layers[spliced_layer].copy()
        adata.layers["unspliced"] = adata.layers[unspliced_layer].copy()
        logger.info(
            "Layer remapping applied: '%s' → 'spliced', '%s' → 'unspliced' (internal use only)",
            spliced_layer,
            unspliced_layer,
        )

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError("Both 'spliced' and 'unspliced' layers are required after layer handling.")

    if gene_type_filter:
        if "gene_type" not in adata.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata = adata[:, adata.var["gene_type"] == gene_type_filter].copy()

    if adata.n_vars == 0:
        raise ValueError("No genes remain after filtering.")

    # Stricter raw counts alignment check (protects against HVG subset / reordering after store_raw_counts)
    raw_gene_list = adata.uns.get("scatrans", {}).get("raw_gene_list")
    if raw_gene_list is not None:
        raw_gene_list = np.asarray(raw_gene_list)
        current_genes = adata.var_names.to_numpy()
        if len(raw_gene_list) == adata.n_vars:
            if not np.array_equal(raw_gene_list, current_genes):
                logger.warning(
                    "Stored raw_gene_list has same length as current genes but different order. "
                    "Raw count column alignment may be invalid for Memento / PyDESeq2 / counts-based methods. "
                    "Consider re-running store_raw_counts() on the current object before analysis."
                )
        else:
            logger.warning(
                "Stored raw_gene_list length (%d) differs from current adata.n_vars (%d). "
                "This likely indicates gene subsetting after store_raw_counts(). "
                "Enrichment universe and raw counts passed to DE backends may be misaligned. "
                "Call store_raw_counts() again on the current (subsetted) object if you want the HVG set as the new universe, "
                "or use the original full-gene adata for DE/enrichment while using a copy for visualization."
            )

    # Mixed model requirements (cell-level RE; pseudobulk + count DE is separate path)
    if use_mixed_model:
        if sample_col is None:
            raise ValueError(
                "sample_col must be provided when use_mixed_model=True (for the random effect grouping)."
            )
        if use_pseudobulk:
            raise ValueError(
                "use_mixed_model=True and use_pseudobulk=True are incompatible. "
                "use_mixed_model (LMM with (1|sample)) is for cell-level data to account for sample correlation. "
                "For pseudobulk, keep use_pseudobulk=True (with pydeseq2 or scanpy) which already aggregates to the sample level. "
                "See README for guidance and dreampy/NEBULA references."
            )

    if "gene_length" not in adata.var.columns:
        adata.var["gene_length"] = np.nan
    if "intron_number" not in adata.var.columns:
        adata.var["intron_number"] = np.nan

    gene_length = pd.to_numeric(adata.var["gene_length"], errors="coerce").to_numpy()
    intron_number = pd.to_numeric(adata.var["intron_number"], errors="coerce").to_numpy()

    adata.var["gene_length"] = gene_length
    adata.var["intron_number"] = intron_number

    valid_feat = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length >= 0)
        & (intron_number >= 0)
    )

    # ==================== PSEUDOBULK (optional) ====================
    is_pseudobulk = False
    if use_pseudobulk:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_pseudobulk=True")
        logger.info("Performing pseudobulk aggregation...")
        adata = _pseudobulk_with_layers(
            adata,
            sample_col,
            groupby,
            x_layer=pb_x_layer,
            use_total_for_x=pb_use_total_for_x,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        is_pseudobulk = True
        adata.obs[groupby] = pd.Categorical(
            adata.obs[groupby].astype(str), categories=[reference_group, target_group]
        )

        n_pb = adata.n_obs
        if n_pb < 5:
            logger.warning(
                "Only %d pseudobulk samples remain after filtering. "
                "Velocity delta estimation and permutation testing will have very low statistical power. "
                "Results (especially active_score and FDR) should be interpreted with extreme caution. "
                "Consider providing more biological replicates or falling back to single-cell mode if appropriate.",
                n_pb,
            )
        elif n_pb < 8:
            logger.info(
                "Only %d pseudobulk samples. Power for detecting differential velocity is limited; "
                "expect many genes to have near-zero velocity_residual and active_score.",
                n_pb,
            )

    # DE preprocess
    if de_preprocess == "normalize_log1p":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif de_preprocess == "auto" and not (is_pseudobulk and pseudobulk_de_backend == "pydeseq2"):
        if "log1p" not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    elif de_preprocess == "none":
        pass

    X_features = (
        np.column_stack([np.log1p(gene_length[valid_feat]), np.log1p(intron_number[valid_feat])])
        if valid_feat.sum() >= 50
        else None
    )

    effective_n_jobs = joblib.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    # ==================== DE ====================
    logger.info("Performing differential expression analysis...")
    if use_memento_de:
        logger.info(
            "Memento (Cell 2024 method-of-moments) selected as main DE backend. "
            "capture_rate=%.4f, num_boot=%d. This replaces scanpy rank_genes_groups for the DE leg of active_score.",
            memento_capture_rate,
            memento_num_boot,
        )

    # Auto-resolve preserved raw counts for count-based DE backends.
    # We warn above if raw_gene_list does not align with current adata.var_names.
    # The wrapper receives the matrix, but remains responsible for safe use/reindexing.
    resolved_counts = None
    needs_raw_counts = use_memento_de or (is_pseudobulk and pseudobulk_de_backend == "pydeseq2")
    if needs_raw_counts:
        if "scatrans" in adata.uns and "raw_gene_list" in adata.uns.get("scatrans", {}):
            if "counts" in adata.layers:
                resolved_counts = adata.layers["counts"]
            elif getattr(adata, "raw", None) is not None:
                resolved_counts = adata.raw.X
        elif "counts" in getattr(adata, "layers", {}):
            resolved_counts = adata.layers["counts"]

    de_df = _run_de_wrapper(
        adata,
        groupby,
        target_group,
        reference_group,
        de_method=de_method,
        is_pseudobulk=is_pseudobulk,
        pb_backend=pseudobulk_de_backend,
        n_jobs=effective_n_jobs,
        strict_pydeseq2_counts=strict_pydeseq2_counts,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        mixed_model_pval=mixed_model_pval,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        memento_n_cpus=memento_n_cpus,
        counts=resolved_counts,
    )

    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]
    if "delta_variance" in de_df.columns:
        adata.var["delta_variance"] = de_df["delta_variance"]
    if "delta_var_pval" in de_df.columns:
        adata.var["delta_var_pval"] = de_df["delta_var_pval"]

    # Surface Memento-specific columns when the memento backend was used (for variability etc.)
    for extra_col in ["memento_de_se", "memento_dv_coef", "memento_dv_se", "memento_dv_pval"]:
        if extra_col in de_df.columns:
            adata.var[extra_col] = de_df[extra_col]

    # ==================== QC: global unspliced fraction (integrated high-value diagnostic) ====================
    unspliced_fraction = np.nan
    try:
        unspliced_fraction = _qc.unspliced_global(
            adata, spliced_key="spliced", unspliced_key="unspliced", warn_threshold=0.5
        )
    except Exception as _e:
        logger.debug("Could not compute global unspliced fraction: %s", _e)

    uns_layer_raw = adata.layers["unspliced"]
    spl_layer_raw = adata.layers["spliced"]

    if is_pseudobulk:
        uns_layer, spl_layer, _, _ = _normalize_velocity_layers_by_size_factor(
            uns_layer_raw, spl_layer_raw
        )
    else:
        uns_layer, spl_layer = uns_layer_raw, spl_layer_raw

    obs_labels = adata.obs[groupby].astype(str).values
    t_mask = obs_labels == target_group
    r_mask = obs_labels == reference_group

    # ==================== VELOCITY DELTA (dual track) ====================
    moments_info: dict[str, Any] = {}
    velocity_layer_for_perm_uns = uns_layer
    velocity_layer_for_perm_spl = spl_layer
    gamma_ref = np.full(adata.n_vars, np.nan)  # will be overwritten in all branches

    if mode == "heuristic":
        delta_velocity, total_us_velocity, gamma_ref = _compute_velocity_delta(
            uns_layer, spl_layer, t_mask, r_mask, prior_weight
        )
        velocity_source = "heuristic_global_ratio"

    elif mode == "advanced":
        adata_comp = adata.copy()
        if is_pseudobulk:
            adata_comp.layers["unspliced"] = uns_layer.copy()
            adata_comp.layers["spliced"] = spl_layer.copy()

        try:
            delta_velocity, total_us_velocity, gamma_ref, moments_info = (
                _compute_moments_velocity_delta(
                    adata_comp,
                    t_mask,
                    r_mask,
                    prior_weight=prior_weight,
                    n_neighbors=advanced_n_neighbors,
                    n_pcs=advanced_n_pcs,
                    use_precomputed=advanced_use_precomputed,
                    recompute_neighbors=advanced_recompute_neighbors,
                    random_state=random_seed,
                )
            )
            velocity_source = "scvelo_moments_groupwise_ratio"
            velocity_layer_for_perm_uns = adata_comp.layers["Mu"].copy()
            velocity_layer_for_perm_spl = adata_comp.layers["Ms"].copy()
            moments_info["advanced_failed"] = False
        except Exception as e:
            if advanced_fallback:
                logger.warning("Advanced mode failed: %s. Falling back to heuristic.", e)
                delta_velocity, total_us_velocity, gamma_ref = _compute_velocity_delta(
                    uns_layer, spl_layer, t_mask, r_mask, prior_weight
                )
                velocity_source = "heuristic_fallback_from_advanced"
                moments_info = {"advanced_failed": True, "failure_reason": str(e)}
            else:
                raise

    # ==================== BIAS CORRECTION (now uses shared implementation) ====================
    total_us_raw = (
        np.asarray(uns_layer_raw.sum(axis=0)).ravel()
        + np.asarray(spl_layer_raw.sum(axis=0)).ravel()
    )
    total_us_raw = np.nan_to_num(total_us_raw)

    valid_expr = total_us_raw >= min_total_counts

    residual, bias_info = fit_huber_bias_correction(
        delta_velocity,
        gene_length,
        intron_number,
        total_us_raw,
        valid_feat,
        valid_expr,
        X_features,
        bias_correction=bias_correction,
    )

    adata.var["velocity_delta_raw"] = delta_velocity
    adata.var["velocity_residual"] = residual
    adata.var["total_us_counts"] = total_us_raw
    adata.var["total_us_counts_raw"] = total_us_raw
    adata.var["total_us_counts_velocity_layer"] = total_us_velocity
    adata.var["valid_expr"] = valid_expr
    adata.var["velocity_source"] = velocity_source

    # effective_gamma is the per-gene reference-group gamma used internally for the delta.
    # It is only exposed to the user when explicitly requested (keeps default output clean
    # and avoids information overload for the basic pipeline).
    if show_effective_gamma:
        adata.var["effective_gamma"] = gamma_ref

    # ==================== DIAGNOSTICS (high priority for usability & paper rigor) ====================
    # (bias fit info is now recorded in diagnostics and passed to the finalizer)
    diagnostics: dict[str, Any] = {
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "n_genes_with_valid_features": int(valid_feat.sum()),
        "unspliced_global_fraction": float(unspliced_fraction)
        if unspliced_fraction is not None
        else np.nan,
        "bias_correction": bias_info,
        "velocity": {
            "source": velocity_source,
            "n_genes_with_finite_delta": int(np.isfinite(delta_velocity).sum()),
            "effective_gamma_exposed": bool(show_effective_gamma),
        },
        "mixed_model": {
            "used": bool(use_mixed_model),
            "sample_col": sample_col if use_mixed_model else None,
            "n_samples": int(adata.obs[sample_col].nunique())
            if (use_mixed_model and sample_col and sample_col in adata.obs.columns)
            else None,
            "delta_variance_available": "delta_variance" in adata.var.columns,
            "median_delta_variance": float(np.nanmedian(adata.var["delta_variance"]))
            if "delta_variance" in adata.var.columns
            else np.nan,
        },
    }
    if mode == "advanced" and moments_info:
        diagnostics["velocity"]["moments"] = {
            k: moments_info.get(k)
            for k in (
                "n_neighbors_effective",
                "n_pcs_effective",
                "used_precomputed_moments",
                "neighbors_source",
            )
            if k in moments_info
        }

    # ==================== SCORING ====================
    lambda_fc = max(_get_exponential_scale_lambda(adata.var["logFC"].values), 0.25)
    lambda_res = max(_get_exponential_scale_lambda(residual), 1e-8)
    lambda_pval = max(
        _get_exponential_scale_lambda(-np.log10(adata.var["p_adj"].values + 1e-300)), 1.0
    )

    s1 = _soft_scale(adata.var["logFC"].values, lambda_fc)
    s2 = _soft_scale(residual, lambda_res)
    s3 = _soft_scale(-np.log10(adata.var["p_adj"].values + 1e-300), lambda_pval)

    total_w = weight_fc + weight_unspliced + weight_pval
    real_score = (weight_fc * s1 + weight_unspliced * s2 + weight_pval * s3) / total_w * 100.0
    adata.var["active_score"] = real_score

    # ==================== PERMUTATION ====================
    current_max_perm = None
    use_fdr_for_significance = True
    # active_fdr_disabled_reason is recorded only for very small permutation spaces
    # (currently informational; not exposed in the public result)

    if use_permutation:
        if is_pseudobulk:
            n_t, n_r = t_mask.sum(), r_mask.sum()
            current_max_perm = float("inf") if n_t + n_r > 30 else max(1, comb(n_t + n_r, n_t) - 1)

        if perm_de_backend == "fast":
            perm_pb_backend, perm_de_method = "scanpy", "t-test_overestim_var"
        elif perm_de_backend == "same":
            perm_pb_backend, perm_de_method = pseudobulk_de_backend, de_method
        else:
            raise ValueError("perm_de_backend must be 'fast' or 'same'")

        # For permutation, default to fast path even when main DE is Memento (performance).
        # Only use Memento in perms if user explicitly requests the advanced (slow) consistent null.
        perm_memento_de = bool(use_memento_de and perm_use_memento_de)
        if perm_memento_de:
            logger.info(
                "perm_use_memento_de=True: permutations will also use Memento (this is slow; "
                "consider leaving it False unless you have strong reasons for a fully consistent null)."
            )

        if (
            is_pseudobulk
            and auto_adjust_n_perm
            and np.isfinite(current_max_perm or 0)
            and (current_max_perm or 0) < n_perm
        ):
            n_perm = int(current_max_perm)

        # Use the extracted single-task (parallel loop stays here for clarity / progress reporting)
        logger.info("Running parallel permutation testing (%d iterations)...", n_perm)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            perm_results = Parallel(n_jobs=effective_n_jobs, backend="loky")(
                delayed(_single_permutation_task)(
                    i + random_seed,
                    obs_labels,
                    target_group,
                    reference_group,
                    adata,
                    X_features,
                    valid_feat,
                    velocity_layer_for_perm_uns,
                    velocity_layer_for_perm_spl,
                    total_us_raw,
                    min_total_counts,
                    weight_fc,
                    weight_unspliced,
                    weight_pval,
                    lambda_fc,
                    lambda_res,
                    lambda_pval,
                    is_pseudobulk,
                    perm_pb_backend,
                    perm_de_method,
                    prior_weight,
                    de_preprocess,
                    strict_pydeseq2_counts,
                    bias_correction=bias_correction,
                    use_memento_de=perm_memento_de,
                    memento_capture_rate=memento_capture_rate,
                    memento_num_boot=memento_num_boot,
                    memento_n_cpus=memento_n_cpus,
                )
                for i in range(n_perm)
            )

        perm_scores_matrix = np.vstack(perm_results)
        exceed_count = np.sum(perm_scores_matrix >= real_score.reshape(1, -1), axis=0)
        pvals = (1.0 + exceed_count) / (n_perm + 1.0)
        adata.var["active_score_pval"] = pvals

        adata.var["active_score_fdr"] = np.ones(adata.n_vars)
        if valid_expr.sum() > 0:
            adata.var.loc[valid_expr, "active_score_fdr"] = multipletests(
                pvals[valid_expr], method="fdr_bh"
            )[1]

        if current_max_perm is not None and current_max_perm < 100:
            use_fdr_for_significance = False
            _ = "small_permutation_space"

    return _finalize_active_score_results(
        adata,
        diagnostics,
        mode=mode,
        velocity_source=velocity_source,
        use_permutation=use_permutation,
        n_perm=n_perm,
        show_plot=show_plot,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        advanced_fallback=advanced_fallback,
        advanced_use_precomputed=advanced_use_precomputed,
        allow_advanced_pseudobulk=allow_advanced_pseudobulk,
        advanced_recompute_neighbors=advanced_recompute_neighbors,
        bias_correction=bias_correction,
        show_effective_gamma=show_effective_gamma,
        is_pseudobulk=is_pseudobulk,
        pval_cutoff=pval_cutoff,
        logfc_cutoff=logfc_cutoff,
        active_fdr_cutoff=active_fdr_cutoff,
        use_fdr_for_significance=use_fdr_for_significance,
        use_delta_variance_pval=use_delta_variance_pval,
        delta_var_pval_cutoff=delta_var_pval_cutoff,
        de_method=de_method,
        pseudobulk_de_backend=pseudobulk_de_backend,
        perm_de_backend=perm_de_backend if use_permutation else None,
        use_memento_de=use_memento_de,
        perm_use_memento_de=perm_use_memento_de if use_permutation else None,
        memento_capture_rate=memento_capture_rate if use_memento_de else None,
        prior_weight=prior_weight,
        min_total_counts=min_total_counts,
        random_seed=random_seed,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        prioritize_velocity=prioritize_velocity,
        weight_fc=weight_fc,
        weight_unspliced=weight_unspliced,
        weight_pval=weight_pval,
    )


# =============================================================================
# Internal helper: output assembly for active_score (keeps the main function cleaner)
# =============================================================================
def _finalize_active_score_results(
    adata: ad.AnnData,
    diagnostics: dict[str, Any],
    *,
    mode: str,
    velocity_source: str,
    use_permutation: bool,
    n_perm: int,
    show_plot: bool,
    **extra_metadata: Any,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Assemble result tables, write rich scatrans metadata (merging safely),
    emit run summary, optionally plot, and return.
    This is extracted to keep active_score() focused on orchestration.
    """
    # Build result columns. By default we keep the output focused on the basic pipeline.
    cols = [
        "active_score",
        "velocity_delta_raw",
        "velocity_residual",
        "logFC",
        "p_val",
        "p_adj",
        "total_us_counts",
        "total_us_counts_raw",
        "total_us_counts_velocity_layer",
        "valid_expr",
        "gene_length",
        "intron_number",
    ]
    if (
        "show_effective_gamma" in extra_metadata
        and extra_metadata.get("show_effective_gamma")
        and "effective_gamma" in adata.var.columns
    ):
        cols.append("effective_gamma")
    if use_permutation:
        cols.extend(["active_score_pval", "active_score_fdr"])
    if "delta_variance" in adata.var.columns:
        cols.append("delta_variance")
    if "delta_var_pval" in adata.var.columns:
        cols.append("delta_var_pval")
    for mc in ["memento_de_se", "memento_dv_coef", "memento_dv_se", "memento_dv_pval"]:
        if mc in adata.var.columns and mc not in cols:
            cols.append(mc)
    cols = [c for c in cols if c in adata.var.columns]

    # (Re-apply significant mask for the final tables — caller may have passed pre-filtered or not)
    # We re-compute here for cleanliness; in practice caller usually passes the already-filtered significant.
    mask = (
        (adata.var["p_adj"] < extra_metadata.get("pval_cutoff", 0.05))
        & (adata.var["logFC"] > extra_metadata.get("logfc_cutoff", 0.5))
        & (adata.var["velocity_residual"] > 0)
        & (adata.var["valid_expr"])
        & (adata.var["active_score"] > 0)
    )
    if (
        use_permutation
        and extra_metadata.get("use_fdr_for_significance", True)
        and "active_score_fdr" in adata.var.columns
    ):
        mask = mask & (
            adata.var["active_score_fdr"] < extra_metadata.get("active_fdr_cutoff", 0.05)
        )

    if (
        "use_delta_variance_pval" in extra_metadata
        and extra_metadata.get("use_delta_variance_pval")
        and "delta_var_pval" in adata.var.columns
    ):
        mask = mask & (
            adata.var["delta_var_pval"] < extra_metadata.get("delta_var_pval_cutoff", 0.05)
        )

    significant = adata.var.loc[mask, cols].copy().sort_values("active_score", ascending=False)
    all_results = adata.var.loc[:, cols].copy().sort_values("active_score", ascending=False)

    logger.info(
        "Analysis completed in %s mode! Significant active genes: %d", mode, len(significant)
    )

    # --- Rich metadata (merge to protect raw_gene_list etc.) ---
    velocity_delta_layer = (
        "scvelo_Mu_Ms_moments"
        if velocity_source.startswith("scvelo_moments")
        else (
            "size_factor_normalized_spliced_unspliced"
            if extra_metadata.get("is_pseudobulk")
            else "raw_spliced_unspliced"
        )
    )

    existing = dict(adata.uns.get("scatrans", {}))
    meta = {
        "version": VERSION,
        "analysis": "active_score",
        "groupby": extra_metadata.get("groupby"),
        "target_group": extra_metadata.get("target_group"),
        "reference_group": extra_metadata.get("reference_group"),
        "mode": mode,
        "velocity_source": velocity_source,
        "velocity_delta_layer": velocity_delta_layer,
        "advanced_fallback": extra_metadata.get("advanced_fallback"),
        "use_permutation": use_permutation,
        "n_perm": int(n_perm) if use_permutation else 0,
        "bias_correction": extra_metadata.get("bias_correction"),
        "show_effective_gamma": extra_metadata.get("show_effective_gamma", False),
        "pval_cutoff": extra_metadata.get("pval_cutoff"),
        "logfc_cutoff": extra_metadata.get("logfc_cutoff"),
        "active_fdr_cutoff": extra_metadata.get("active_fdr_cutoff"),
        "use_fdr_for_significance": extra_metadata.get("use_fdr_for_significance"),
        "use_delta_variance_pval": extra_metadata.get("use_delta_variance_pval"),
        "delta_var_pval_cutoff": extra_metadata.get("delta_var_pval_cutoff"),
        "de_method": extra_metadata.get("de_method"),
        "pseudobulk_de_backend": extra_metadata.get("pseudobulk_de_backend"),
        "perm_de_backend": extra_metadata.get("perm_de_backend"),
        "use_memento_de": extra_metadata.get("use_memento_de"),
        "perm_use_memento_de": extra_metadata.get("perm_use_memento_de"),
        "memento_capture_rate": extra_metadata.get("memento_capture_rate"),
        "prior_weight": extra_metadata.get("prior_weight"),
        "min_total_counts": extra_metadata.get("min_total_counts"),
        "random_seed": extra_metadata.get("random_seed"),
        "use_mixed_model": extra_metadata.get("use_mixed_model"),
        "sample_col": extra_metadata.get("sample_col"),
        "prioritize_velocity": extra_metadata.get("prioritize_velocity"),
        "weight_fc": extra_metadata.get("weight_fc"),
        "weight_unspliced": extra_metadata.get("weight_unspliced"),
        "weight_pval": extra_metadata.get("weight_pval"),
        "diagnostics": diagnostics,
        "unspliced_global_fraction": diagnostics.get("unspliced_global_fraction", np.nan),
    }

    if use_permutation:
        note = (
            "For efficiency, velocity layers and effective_gamma are fixed from the original data. "
            "Group labels are shuffled to recompute the composite active_score."
        )
        if extra_metadata.get("use_memento_de") and not extra_metadata.get("perm_use_memento_de"):
            note += (
                " Memento was used for the observed DE; permutations used the configured fast/same "
                "non-Memento DE backend unless perm_use_memento_de=True."
            )
        elif extra_metadata.get("use_memento_de") and extra_metadata.get("perm_use_memento_de"):
            note += " Memento was used for both observed DE and permutation null."
        meta["permutation_approximation_note"] = note

    existing.update({k: v for k, v in meta.items() if v is not None})
    adata.uns["scatrans"] = existing

    # Display mode defined early (outside try) to avoid UnboundLocalError in plotting
    display_mode = (
        "advanced→heuristic fallback"
        if velocity_source == "heuristic_fallback_from_advanced"
        else mode
    )

    # --- User-facing run summary ---
    try:
        ufrac = diagnostics.get("unspliced_global_fraction", np.nan)
        bias = diagnostics.get("bias_correction", {})
        n_fit = bias.get("n_genes_used_for_fit", 0)
        fb = " (median fallback)" if bias.get("fallback_to_median") else ""
        logger.info(
            "Run summary — cells: %d | unspliced frac: %.1f%% | bias fit genes: %d%s | mode: %s | sig: %d",
            diagnostics.get("n_cells", 0),
            (ufrac * 100.0) if np.isfinite(ufrac) else float("nan"),
            n_fit,
            fb,
            display_mode,
            len(significant),
        )
        if use_permutation:
            logger.info(
                "Permutation used %d iterations (velocity layers fixed from original labeling).",
                n_perm,
            )
    except Exception:
        pass

    # --- Optional plotting (delegated) ---
    if show_plot:
        try:
            from . import pl

            pl.comet_plot(
                all_results,
                top_n=12,
                title=f"scATrans Active Drivers ({display_mode})",
            )
        except Exception:
            logger.debug("show_plot=True but plotting failed (missing optional deps or display).")

    return adata, significant, all_results


def differential_expression(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "GA",
    reference_group: str = "Ctrl",
    subset_col: str | None = None,
    subset_values: str | list[str] | tuple[str, ...] | None = None,
    de_method: str = "t-test_overestim_var",
    pseudobulk_de_backend: str = "pydeseq2",
    use_pseudobulk: bool = False,
    sample_col: str | None = None,
    min_cells: int = 10,
    min_counts: int = 1000,
    pb_x_layer: str = "X",  # for pseudobulk, what to aggregate (usually the count matrix)
    pb_use_total_for_x: bool = True,
    de_preprocess: str = "auto",
    min_total_counts: int = 50,  # reserved for API compatibility / future use; currently not enforced in DE path
    strict_pydeseq2_counts: bool = True,
    use_mixed_model: bool = False,
    use_delta_variance_pval: bool = False,
    delta_var_pval_cutoff: float = 0.05,
    mixed_model_pval: str = "wald",
    # Memento support (first-class, integrated backend)
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    # Advanced permutation support (rarely needed for pure DE)
    use_permutation: bool = False,
    perm_de_backend: str = "fast",
    n_perm: int = 100,
    active_fdr_cutoff: float = 0.05,  # reserved for API compat with active_score; ignored here
    random_seed: int = 42,  # reserved for API compat
    n_jobs: int = -1,
    gene_type_filter: str | None = None,
    # Allow providing raw counts separately when adata.X is already HVG+log (very common)
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData | None = None,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Standalone differential expression (DE) using the same flexible backends
    as scATrans (scanpy methods, PyDESeq2 pseudobulk, mixed linear models,
    and Memento -- the Cell 2024 method-of-moments framework).

    This function does **not** require spliced/unspliced (velocity) layers.
    It is intended for users who want high-quality DE (especially via Memento),
    followed by scATrans' downstream tools:

        candidates = scat.filter_active_genes(de_results, ...)
        # For enrichment, pass adata= (if store_raw_counts was used) so it uses
        # the preserved full measured gene set as universe, not just current HVGs.
        enrich = scat.run_enrichment(candidates.index.tolist(), ..., adata=adata)
        scat.pl.volcano_plot(de_results, ...)
        scat.pl.enrich_dotplot(enrich, ...)

    All DE-related options from `active_score` are supported here
    (pseudobulk, mixed models, Memento, etc.).

    Returns
    -------
    (adata_with_results, results_df)
        - results_df is a ranked DataFrame (by |logFC| or p_adj) containing
          at minimum: logFC, p_val, p_adj, and (when use_memento_de) the
          native memento_de_* / memento_dv_* columns.
        - adata.var is updated with the same columns for convenience.
        - Metadata is stored under adata.uns["scatrans"].
    """
    # --- minimal shared validation (subset + group checks) ---
    if subset_col is not None:
        if subset_col not in adata_input.obs.columns:
            raise ValueError(f"subset_col='{subset_col}' not found in adata.obs.columns")
        if subset_values is None:
            raise ValueError("subset_values must be provided when subset_col is specified")
        if isinstance(subset_values, (str, int, float)):
            subset_values_list = [str(subset_values)]
        else:
            subset_values_list = [str(v) for v in subset_values]
        subset_mask = adata_input.obs[subset_col].astype(str).isin(subset_values_list)
        adata_input = adata_input[subset_mask].copy()
        if adata_input.n_obs == 0:
            raise ValueError("No cells remain after subsetting.")

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    target_group = str(target_group)
    reference_group = str(reference_group)
    if target_group == reference_group:
        raise ValueError("target_group and reference_group must be different.")
    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")
    if target_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"target_group '{target_group}' not found.")
    if reference_group not in adata_input.obs[groupby].astype(str).unique():
        raise ValueError(f"reference_group '{reference_group}' not found.")

    if gene_type_filter:
        if "gene_type" not in adata_input.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata_input = adata_input[:, adata_input.var["gene_type"] == gene_type_filter].copy()

    if adata_input.n_vars == 0:
        raise ValueError("No genes remain after filtering.")

    if use_mixed_model and sample_col is None:
        raise ValueError("sample_col must be provided when use_mixed_model=True")

    if use_pseudobulk and sample_col is None:
        raise ValueError("sample_col must be provided when use_pseudobulk=True")

    # Memento-specific guard (same as in active_score)
    if use_memento_de and use_pseudobulk:
        raise ValueError(
            "use_memento_de=True is not supported with use_pseudobulk=True "
            "(Memento is a cell-level method-of-moments estimator)."
        )

    # Memento requires count data; force no log-norm preprocess for the DE leg
    de_preprocess = _coerce_memento_de_preprocess(use_memento_de, de_preprocess)

    # Shared DE option validation (deduplicated via helper)
    _validate_de_common_options(
        de_preprocess=de_preprocess,
        pseudobulk_de_backend=pseudobulk_de_backend,
        n_jobs=n_jobs,
        use_permutation=use_permutation,
        n_perm=n_perm,
        use_mixed_model=use_mixed_model,
        mixed_model_pval=mixed_model_pval,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        min_cells=min_cells,
        min_counts=min_counts,
    )

    # Auto-resolve preserved raw counts from scatrans metadata if user called store_raw_counts early
    # This makes Memento / PyDESeq2 paths more convenient by auto-supplying preserved
    # raw counts (when store_raw_counts was called) without the user having to pass
    # counts= explicitly. Alignment warnings are emitted earlier if the stored gene
    # list no longer matches the current adata.
    if counts is None and (
        use_memento_de or (use_pseudobulk and pseudobulk_de_backend == "pydeseq2")
    ):
        if "scatrans" in adata_input.uns and "raw_gene_list" in adata_input.uns.get("scatrans", {}):
            # Prefer the layer if present (it has the actual count matrix for those genes)
            if "counts" in adata_input.layers:
                counts = adata_input.layers["counts"]
            elif getattr(adata_input, "raw", None) is not None:
                counts = adata_input.raw.X
            # else the gene list is available but matrix may have been lost; wrapper will warn
        elif "counts" in getattr(adata_input, "layers", {}):
            counts = adata_input.layers["counts"]

    # Note on raw counts: users should call scat.store_raw_counts(adata) early.
    # The DE backends will use layers[layer] or adata.raw when available.

    # --- prepare data (pseudobulk if requested) ---
    adata = adata_input.copy()

    if use_pseudobulk:
        logger.info("Performing pseudobulk aggregation for DE...")
        adata = _pseudobulk_with_layers(
            adata,
            sample_col,
            groupby,
            x_layer=pb_x_layer if pb_x_layer != "X" else None,
            use_total_for_x=pb_use_total_for_x,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        adata.obs[groupby] = pd.Categorical(
            adata.obs[groupby].astype(str), categories=[reference_group, target_group]
        )

    # DE preprocess
    if use_memento_de and de_preprocess != "none":
        logger.info(
            "use_memento_de=True: forcing de_preprocess='none' (Memento works on raw counts)."
        )
        de_preprocess = "none"

    if de_preprocess == "normalize_log1p":
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif de_preprocess == "auto" and not (use_pseudobulk and pseudobulk_de_backend == "pydeseq2"):
        if "log1p" not in adata.uns:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
    elif de_preprocess == "none":
        pass

    effective_n_jobs = joblib.cpu_count() if n_jobs == -1 else max(1, n_jobs)

    # --- run DE via the shared engine (Memento, scanpy, DESeq2, mixedlm all supported) ---
    logger.info("Performing differential expression analysis (differential_expression mode)...")
    de_df = _run_de_wrapper(
        adata,
        groupby,
        target_group,
        reference_group,
        de_method=de_method,
        is_pseudobulk=use_pseudobulk,
        pb_backend=pseudobulk_de_backend,
        n_jobs=effective_n_jobs,
        strict_pydeseq2_counts=strict_pydeseq2_counts,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        mixed_model_pval=mixed_model_pval,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        memento_n_cpus=memento_n_cpus,
        counts=counts,
    )

    # Store results
    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]

    for extra in [
        "delta_variance",
        "delta_var_pval",
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
    ]:
        if extra in de_df.columns:
            adata.var[extra] = de_df[extra]

    # Build clean results table (no velocity columns)
    cols = ["logFC", "p_val", "p_adj"]
    for c in [
        "delta_variance",
        "delta_var_pval",
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
    ]:
        if c in adata.var.columns:
            cols.append(c)

    # Add a simple base expression measure when possible
    if "total_us_counts" in adata.var.columns:
        cols.append("total_us_counts")
    else:
        # fallback: mean of current X (after any preprocess the user chose)
        try:
            means = np.asarray(adata.X.mean(axis=0)).ravel()
            adata.var["baseMean"] = means
            cols.append("baseMean")
        except Exception:
            pass

    results = adata.var[cols].copy()
    results = results.sort_values("p_adj", ascending=True)

    # Optional permutation-based FDR on the DE p-values themselves
    if use_permutation:
        raise NotImplementedError(
            "use_permutation=True is not supported in differential_expression(). "
            "Permutation FDR on pure DE is not implemented. "
            "Use active_score(..., use_permutation=True) for composite permutation-based FDR "
            "when you have velocity layers, or compute your own permutation test on the p-values."
        )

    # Metadata — merge to preserve raw_gene_list etc. from store_raw_counts()
    existing = dict(adata.uns.get("scatrans", {}))
    existing.update(
        {
            "mode": "differential_expression",
            "version": VERSION,
            "groupby": groupby,
            "target_group": target_group,
            "reference_group": reference_group,
            "use_pseudobulk": use_pseudobulk,
            "use_mixed_model": use_mixed_model,
            "use_memento_de": use_memento_de,
            "memento_capture_rate": memento_capture_rate if use_memento_de else None,
            "de_method": de_method,
            "pseudobulk_de_backend": pseudobulk_de_backend,
            "use_permutation": use_permutation,
            "n_perm": n_perm if use_permutation else 0,
            "de_preprocess": de_preprocess,
            "strict_pydeseq2_counts": strict_pydeseq2_counts,
            "min_cells": min_cells,
            "min_counts": min_counts,
            "min_total_counts": min_total_counts,
            "sample_col": sample_col,
            "pb_x_layer": pb_x_layer,
            "pb_use_total_for_x": pb_use_total_for_x,
            "use_delta_variance_pval": use_delta_variance_pval,
            "delta_var_pval_cutoff": delta_var_pval_cutoff,
            "mixed_model_pval": mixed_model_pval if use_mixed_model else None,
            "memento_num_boot": memento_num_boot if use_memento_de else None,
            "memento_n_cpus": memento_n_cpus if use_memento_de else None,
            "perm_de_backend": perm_de_backend if use_permutation else None,
            "active_fdr_cutoff": active_fdr_cutoff if use_permutation else None,
            "random_seed": random_seed,
            "n_jobs": n_jobs,
            "gene_type_filter": gene_type_filter,
        }
    )
    adata.uns["scatrans"] = existing

    logger.info("DE completed. %d genes in results table.", len(results))
    return adata, results


def store_raw_counts(
    adata: Any, layer: str = "counts", save_raw: bool = False, overwrite: bool = False
) -> None:
    """
    Store raw counts and the original spliced/unspliced (or mature/nascent) layers
    early in the analysis, right after loading and basic QC, but BEFORE HVG selection,
    normalization, or log1p.

    This is critical for scATrans because:
    - Memento and PyDESeq2 need raw counts for proper modeling.
    - Velocity / active transcription analysis (active_score) needs the original
      spliced/unspliced matrices on as many genes as possible.

    By default we only save to the given layer (save_raw defaults to False so we do
    not automatically touch adata.raw unless you explicitly ask for it).

    We automatically save any existing velocity layers under "raw_spliced",
    "raw_unspliced" (or "raw_mature", "raw_nascent"). These raw_* layers are
    subject to the normal AnnData behavior: if you later gene-subset the object
    (e.g. to HVGs), the layers are subsetted as well. They do **not** magically
    retain the original full-gene matrices after subsetting.

    If you need the full-gene raw velocity data after HVG-based visualization,
    either:
      - call store_raw_counts() on the full object and keep the full object for
        DE / active_score / enrichment while using a copy for visualization, or
      - use save_raw=True (which sets adata.raw).

    Recommended early call:
        scat.store_raw_counts(adata, layer="counts", save_raw=False)
    """
    if layer in adata.layers and not overwrite:
        if _is_integer_counts_like(adata.layers[layer]):
            logger.debug(f"Layer '{layer}' already exists with integer counts; not overwriting.")
        else:
            logger.warning(
                f"Existing layer '{layer}' does not look like raw counts; "
                "pass overwrite=True to replace it."
            )
    else:
        if not _is_integer_counts_like(adata.X):
            logger.warning(
                "Current adata.X does not look like raw integer counts. "
                "store_raw_counts should be called early (after basic QC, before normalize/log1p/HVG)."
            )
        adata.layers[layer] = adata.X.copy()
        logger.info(f"Saved raw counts to adata.layers['{layer}'].")

    # Save the gene list at this moment as the measured universe for later enrichment
    # and count-based analyses. Unlike layers, this metadata list is kept in .uns and
    # therefore is not automatically subsetted when adata is later subsetted to HVGs.
    # Downstream code must still verify that any stored count matrix is aligned with
    # the current adata.var_names before using it.
    if "scatrans" not in adata.uns:
        adata.uns["scatrans"] = {}
    adata.uns["scatrans"]["raw_gene_list"] = list(adata.var_names)
    logger.info(
        "Saved the current gene list as the measured universe for enrichment (in adata.uns['scatrans']['raw_gene_list'])."
    )

    if save_raw:
        if getattr(adata, "raw", None) is not None and not overwrite:
            logger.debug("adata.raw already exists; skipping (pass overwrite=True to replace).")
        else:
            adata.raw = adata.copy()
            logger.info("Set adata.raw to preserve full data.")

    # Preserve original velocity layers (spliced/unspliced or mature/nascent)
    # under "raw_*" names for the current gene set at the time of the call.
    # Note: these are normal AnnData layers; if the adata is later gene-subsetted
    # (e.g. to HVGs), the layers will be subsetted as well. They do not retain
    # the original full-gene matrices after subsetting. Use save_raw=True or
    # keep a separate full-gene object if full-gene recovery is required.
    for vel_name in ("spliced", "unspliced", "mature", "nascent"):
        if vel_name in adata.layers:
            raw_vel_name = f"raw_{vel_name}"
            if raw_vel_name not in adata.layers or overwrite:
                adata.layers[raw_vel_name] = adata.layers[vel_name].copy()
                logger.info(f"Saved original {vel_name} to adata.layers['{raw_vel_name}'].")


def restore_raw_counts(adata: Any, layer: str = "counts", inplace: bool = False) -> Any | None:
    """
    Restore raw counts from the stored layer (preferred) or adata.raw back into .X.

    This is useful when you have done HVG + log1p on .X for visualization,
    but want to work with (or pass to other tools) the raw counts for the
    genes currently in the adata (or the preserved set).

    It only uses explicitly stored raw data (from store_raw_counts), never
    attempts to recover from log-transformed data.

    Parameters
    ----------
    adata : AnnData
        The AnnData object.
    layer : str
        The layer name where raw counts were stored (default "counts").
    inplace : bool
        If True, modify adata in place and return None.
        If False (default), return a new AnnData with .X set to raw counts.

    Returns
    -------
    AnnData or None
        If not inplace, a copy of adata with raw counts in .X.
    """
    if layer in adata.layers:
        raw = adata.layers[layer].copy()
        source = f"layers['{layer}']"
    elif getattr(adata, "raw", None) is not None:
        raw = adata.raw.X.copy()
        source = "adata.raw"
    else:
        raise ValueError(
            f"No raw counts found in layer '{layer}' or adata.raw. "
            "Call scat.store_raw_counts(adata) early to preserve them."
        )

    if raw.shape[1] != adata.n_vars:
        raise ValueError(
            f"Stored raw counts have {raw.shape[1]} genes, but current adata has {adata.n_vars} genes. "
            "Cannot restore into .X without explicit gene reindexing. "
            "Use the object before gene subsetting, or call store_raw_counts() again on the current object."
        )

    # Additional guard when restoring from adata.raw (same dimension but possibly different order/names)
    if (
        source == "adata.raw"
        and hasattr(adata.raw, "var_names")
        and not np.array_equal(adata.raw.var_names, adata.var_names)
    ):
        raise ValueError(
            "adata.raw has the same number of genes as current adata, but gene names/order differ. "
            "Cannot restore into .X without explicit gene reindexing."
        )

    if inplace:
        adata.X = raw
        logger.info(f"Restored raw counts from {source} into adata.X (inplace).")
        return None
    else:
        adata_restored = adata.copy()
        adata_restored.X = raw
        logger.info(f"Created copy with raw counts from {source} in .X.")
        return adata_restored


_NOT_PROVIDED = object()


def filter_active_genes(
    results: pd.DataFrame,
    *,
    preset: str | None = None,
    active_score_cutoff: Any = _NOT_PROVIDED,
    pval_cutoff: Any = _NOT_PROVIDED,
    velocity_residual_cutoff: Any = _NOT_PROVIDED,
    logfc_cutoff: Any = _NOT_PROVIDED,
    active_score_fdr_cutoff: Any = _NOT_PROVIDED,
    effective_gamma_min: Any = _NOT_PROVIDED,
    effective_gamma_max: Any = _NOT_PROVIDED,
    delta_variance_min: Any = _NOT_PROVIDED,
) -> pd.DataFrame:
    """Apply custom post-filtering to a results DataFrame (from `active_score` or `differential_expression`).

    This helper works for both:
    - Full `active_score` output (has `active_score` + velocity_residual).
    - Pure DE results from `differential_expression` (only logFC / p_adj + optional memento columns).

    It standardizes the common workflow:
    1. Run `active_score(...)` or `differential_expression(...)`.
    2. Use this function on the returned table to derive a final gene list.

    The function supports `preset` to automatically select reasonable default thresholds
    for different analysis modes:

    - preset="heuristic": stricter defaults suitable for typical single-cell data with default weights
      (active_score >= 55, velocity_residual > 1.0, etc.).
    - preset="pseudobulk": more lenient defaults that account for the much smaller
      magnitude of velocity_residual and active_score after sample-level aggregation
      (active_score >= 5, velocity_residual > 0.05, logFC > 0.2, etc.).
    - preset=None (or "permissive"/"none"): apply only explicitly provided cutoffs; this is the most
      permissive / backward-compatible mode and returns nearly the full table (subject to any
      user-supplied thresholds).

    Presets are oriented toward target-group "activated" / upregulated signals (positive
    logFC + positive velocity_residual). For pure two-sided DE gene selection from
    differential_expression() results, use preset=None and supply your own logfc_cutoff
    (e.g. a negative value or use abs() filtering yourself after the call).

    If you explicitly pass any cutoff parameter, it takes precedence over the preset.

    Calling with no arguments (or only the DataFrame) and no preset returns the full
    `all_results` (fully permissive).

    Only filters corresponding to columns present in the DataFrame are applied.
    This is safe whether or not `use_permutation=True` or `use_mixed_model=True` was used.

    Parameters
    ----------
    results : pd.DataFrame
        The `all_results` table returned as the third element of `active_score`.
    preset : str or None
        One of "heuristic", "pseudobulk", "permissive", "none".
        When provided, supplies recommended cutoff values for that analysis style
        for any parameters you did not explicitly pass.
    active_score_cutoff : float
        Minimum composite active transcription score (0-100).
    pval_cutoff : float
        Maximum adjusted p-value if a "p_adj" column is present in the results DataFrame;
        otherwise falls back to the nominal "p_val" column. This matches the behavior
        of the internal significant mask in active_score() and common user expectations
        (FDR control when available).
    velocity_residual_cutoff : float
        Minimum bias-corrected velocity residual.
    logfc_cutoff : float
        Minimum log fold change.
    active_score_fdr_cutoff : float or None
        If the column exists (use_permutation=True), max permutation FDR on active_score.
    effective_gamma_min : float
        Minimum reference-group effective gamma. See README section on effective_gamma.
    effective_gamma_max : float or None
        Optional upper bound on effective_gamma.
    delta_variance_min : float or None
        If the column exists (use_mixed_model=True), minimum variance fraction
        explained by condition.

    Returns
    -------
    pd.DataFrame
        Subset of the input, sorted by active_score descending.
    """
    if not isinstance(results, pd.DataFrame):
        raise ValueError("results must be the all_results DataFrame returned by active_score")

    # Resolve values from preset + explicit overrides
    if preset is not None:
        p = preset.lower()
        if p in ("heuristic", "single_cell", "default"):
            preset_vals = {
                "active_score_cutoff": 55.0,
                "pval_cutoff": 0.05,
                "velocity_residual_cutoff": 1.0,
                "logfc_cutoff": 0.35,
                "active_score_fdr_cutoff": 0.25,
                "effective_gamma_min": 0.05,
                "effective_gamma_max": 1.0,
                "delta_variance_min": None,
            }
        elif p in ("pseudobulk", "bulk"):
            preset_vals = {
                "active_score_cutoff": 5.0,
                "pval_cutoff": 0.05,
                "velocity_residual_cutoff": 0.05,
                "logfc_cutoff": 0.2,
                "active_score_fdr_cutoff": 0.25,
                "effective_gamma_min": 0.05,
                "effective_gamma_max": 1.0,
                "delta_variance_min": None,
            }
        elif p in ("permissive", "none", "all", "no_filter"):
            preset_vals = {
                "active_score_cutoff": 0.0,
                "pval_cutoff": 1.0,
                "velocity_residual_cutoff": float("-inf"),
                "logfc_cutoff": float("-inf"),
                "active_score_fdr_cutoff": 1.0,
                "effective_gamma_min": float("-inf"),
                "effective_gamma_max": None,
                "delta_variance_min": None,
            }
        else:
            raise ValueError(
                f"Unknown preset '{preset}'. "
                "Valid presets: 'heuristic', 'pseudobulk', 'permissive'."
            )
    else:
        preset_vals = {}

    # Apply preset only where user did not explicitly provide a value
    def _resolve(name: str, current: Any, default: Any) -> Any:
        if current is not _NOT_PROVIDED:
            return current
        return preset_vals.get(name, default)

    active_score_cutoff = _resolve("active_score_cutoff", active_score_cutoff, 0.0)
    pval_cutoff = _resolve("pval_cutoff", pval_cutoff, 1.0)
    velocity_residual_cutoff = _resolve(
        "velocity_residual_cutoff", velocity_residual_cutoff, float("-inf")
    )
    logfc_cutoff = _resolve("logfc_cutoff", logfc_cutoff, float("-inf"))
    active_score_fdr_cutoff = _resolve("active_score_fdr_cutoff", active_score_fdr_cutoff, 1.0)
    effective_gamma_min = _resolve("effective_gamma_min", effective_gamma_min, float("-inf"))
    effective_gamma_max = _resolve("effective_gamma_max", effective_gamma_max, None)
    delta_variance_min = _resolve("delta_variance_min", delta_variance_min, None)

    df = results.copy()
    mask = pd.Series(True, index=df.index)

    # Core filters
    # Prefer adjusted p-value when present (consistent with active_score internal significant mask
    # and common user expectation). Fall back to nominal p_val only if p_adj is absent.
    if "active_score" in df.columns:
        mask &= df["active_score"] >= active_score_cutoff
    if "p_adj" in df.columns:
        mask &= df["p_adj"] < pval_cutoff
    elif "p_val" in df.columns:
        mask &= df["p_val"] < pval_cutoff
    if "velocity_residual" in df.columns:
        mask &= df["velocity_residual"] > velocity_residual_cutoff
    if "logFC" in df.columns:
        mask &= df["logFC"] > logfc_cutoff

    # Permutation FDR
    if active_score_fdr_cutoff is not None and "active_score_fdr" in df.columns:
        mask &= df["active_score_fdr"] < active_score_fdr_cutoff

    # effective_gamma
    if "effective_gamma" in df.columns:
        gamma = df["effective_gamma"]
        mask &= gamma.notna() & (gamma > effective_gamma_min)
        if effective_gamma_max is not None:
            mask &= gamma < effective_gamma_max

    # Delta variance
    if delta_variance_min is not None and "delta_variance" in df.columns:
        mask &= df["delta_variance"] >= delta_variance_min

    filtered = df[mask].copy()

    # Sorting: prefer active_score when present (velocity + DE composite),
    # otherwise fall back to p_adj, then positive logFC (for pure DE results,
    # consistent with target-group activation bias of the helper).
    if "active_score" in filtered.columns:
        filtered = filtered.sort_values("active_score", ascending=False)
    elif "p_adj" in filtered.columns:
        sort_cols = ["p_adj"]
        ascending = [True]
        if "logFC" in filtered.columns:
            sort_cols.append("logFC")
            ascending.append(False)
        filtered = filtered.sort_values(sort_cols, ascending=ascending)
    else:
        filtered = filtered.sort_values(filtered.columns[0])
    return filtered


def diagnose_design(
    adata_input: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None = None,
    min_cells_per_sample: int = 10,  # reserved; not yet used to filter samples or adjust power estimates
) -> dict[str, Any]:
    """
    Analyze the experimental design and provide guidance on suitable analysis choices
    and expected power/limitations.

    This is intended as a pre-flight or post-subset diagnostic to help users
    interpret warnings and choose between single-cell, pseudobulk, or mixed-model paths.

    Returns a dictionary with keys:
      - n_cells_target, n_cells_reference
      - n_samples_target, n_samples_reference (if sample_col provided)
      - unspliced_global_fraction
      - recommendations: list of human-readable strings
      - warnings: list of human-readable strings
      - suggested_preset: "heuristic", "pseudobulk", or None
    """
    adata = adata_input.copy()

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs")

    target_mask = adata.obs[groupby].astype(str) == str(target_group)
    ref_mask = adata.obs[groupby].astype(str) == str(reference_group)

    n_t = int(target_mask.sum())
    n_r = int(ref_mask.sum())

    result: dict[str, Any] = {
        "n_cells_target": n_t,
        "n_cells_reference": n_r,
        "n_samples_target": None,
        "n_samples_reference": None,
        "unspliced_global_fraction": None,
        "recommendations": [],
        "warnings": [],
        "suggested_preset": None,
    }

    # Global unspliced fraction (important technical QC)
    try:
        ufrac = _qc.unspliced_global(
            adata, spliced_key="spliced", unspliced_key="unspliced", warn_threshold=0.5
        )
        result["unspliced_global_fraction"] = float(ufrac)
        if ufrac > 0.5:
            result["warnings"].append(
                f"Global unspliced fraction is high ({ufrac:.1%}). "
                "This often indicates nuclear enrichment or gDNA contamination and can "
                "reduce the reliability of velocity-based signals."
            )
    except Exception:
        pass

    # Sample structure
    if sample_col and sample_col in adata.obs.columns:
        n_s_t = adata.obs.loc[target_mask, sample_col].nunique()
        n_s_r = adata.obs.loc[ref_mask, sample_col].nunique()
        result["n_samples_target"] = int(n_s_t)
        result["n_samples_reference"] = int(n_s_r)

        if min(n_s_t, n_s_r) < 3:
            result["warnings"].append(
                f"Very few biological samples per group (target={n_s_t}, reference={n_s_r}). "
                "Pseudobulk aggregation will have extremely low power for velocity delta. "
                "Consider using the cell-level mixed-model path (use_mixed_model=True) "
                "instead of use_pseudobulk=True, or interpret results with extreme caution."
            )
        elif min(n_s_t, n_s_r) < 5:
            result["warnings"].append(
                f"Small number of biological samples per group (target={n_s_t}, reference={n_s_r}). "
                "Power for detecting differential nascent RNA excess will be limited. "
                "Permutation-based FDR (if used) will also have reduced reliability."
            )
            result["suggested_preset"] = "pseudobulk"

        result["recommendations"].append(
            "With multiple samples per group, both pseudobulk (with PyDESeq2 or scanpy) "
            "and cell-level mixed model (use_mixed_model=True) are viable. "
            "See the small-sample guidance in the documentation."
        )
    else:
        result["recommendations"].append(
            "No sample_col provided. The analysis will treat cells as independent. "
            "If cells come from multiple biological replicates, consider providing sample_col "
            "and using either use_pseudobulk=True or use_mixed_model=True to avoid "
            "pseudoreplication."
        )

    # Very small total cell numbers
    if min(n_t, n_r) < 50:
        result["warnings"].append(
            f"Very small number of cells in at least one group (target={n_t}, reference={n_r}). "
            "Velocity delta estimation and any downstream permutation testing will have low power."
        )

    # Suggest preset
    if result["suggested_preset"] is None:
        if sample_col and result.get("n_samples_target", 0) >= 5:
            result["suggested_preset"] = "pseudobulk"
        else:
            result["suggested_preset"] = "heuristic"

    # General advice
    result["recommendations"].append(
        "After running active_score, always inspect adata.uns['scatrans']['diagnostics'] "
        "and the distributions in the returned all_results DataFrame before applying cutoffs."
    )

    # Print a concise user-facing summary
    logger.info("Design diagnosis:")
    logger.info("  Cells — target: %d | reference: %d", n_t, n_r)
    if result["n_samples_target"] is not None:
        logger.info(
            "  Samples — target: %d | reference: %d",
            result["n_samples_target"],
            result["n_samples_reference"],
        )
    if result["unspliced_global_fraction"] is not None:
        logger.info(
            "  Global unspliced fraction: %.1f%%", result["unspliced_global_fraction"] * 100
        )

    for w in result["warnings"]:
        logger.warning("  [WARNING] %s", w)
    for r in result["recommendations"]:
        logger.info("  [RECOMMENDATION] %s", r)

    return result
