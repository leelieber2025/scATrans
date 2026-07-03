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
import math
import warnings
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sparse  # for type hints in signatures (e.g. spmatrix)

# qc is imported lazily inside active_score to keep startup light, but exposed at package level
from . import qc as _qc  # for unspliced_global integration
from ._de import _run_de_wrapper
from ._permutation import run_permutation_test
from ._utils import (
    LEGACY_VELOCITY_RESIDUAL_COL,
    UNSPLICED_EXCESS_DELTA_COL,
    UNSPLICED_EXCESS_FDR_COL,
    UNSPLICED_EXCESS_PVAL_COL,
    UNSPLICED_EXCESS_RESIDUAL_COL,
    _apply_de_preprocess,
    _clear_log_preprocess_metadata,
    _get_exponential_scale_lambda,
    _is_integer_counts_like,
    _normalize_group_label,
    _normalize_velocity_layers_by_size_factor,
    _pseudobulk_with_layers,
    _resolve_aligned_raw_counts,
    _soft_scale,
    _validate_group_contrast,
    _write_unspliced_excess_columns,
    comb,  # for small-n permutation space calculation
)
from ._utils import _fit_huber_bias_correction as fit_huber_bias_correction
from ._velocity import _compute_moments_velocity_delta, _compute_velocity_delta

try:
    from . import _version

    VERSION = _version.version
except (ImportError, AttributeError):
    VERSION = "0.9.2"

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _materialize_if_view(adata: ad.AnnData) -> ad.AnnData:
    """Return a writable AnnData when ``adata`` is an obs/var view (needed before in-place ``.obs`` writes)."""
    if getattr(adata, "is_view", False):
        return adata.copy()
    return adata


def _select_obs(
    adata: ad.AnnData,
    mask: pd.Series,
    *,
    copy_input: bool,
) -> ad.AnnData:
    """Subset to obs ``mask``; call ``AnnData.copy()`` only when ``copy_input=True``."""
    mask = mask.reindex(adata.obs_names, fill_value=False)
    if copy_input:
        return adata[mask].copy()
    if bool(mask.all()):
        return adata
    return adata[mask]


def _select_var(
    adata: ad.AnnData,
    mask: pd.Series,
    *,
    copy_input: bool,
) -> ad.AnnData:
    """Subset to var ``mask``; call ``AnnData.copy()`` only when ``copy_input=True``."""
    if copy_input:
        return adata[:, mask].copy()
    if bool(mask.all()):
        return adata
    return adata[:, mask]


def _resolve_velocity_layer_keys(adata: ad.AnnData) -> tuple[str, str] | None:
    """Return (spliced_key, unspliced_key) for QC, including kb_python mature/nascent."""
    layers = set(adata.layers.keys())
    if "spliced" in layers and "unspliced" in layers:
        return "spliced", "unspliced"
    if "mature" in layers and "nascent" in layers:
        return "mature", "nascent"
    return None


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

    if use_mixed_model and mixed_model_pval not in ("wald", "lrt"):
        raise ValueError("mixed_model_pval must be 'wald' or 'lrt'.")

    if use_memento_de:
        if not (0 < memento_capture_rate < 1):
            raise ValueError(
                "memento_capture_rate must be in (0, 1). "
                "Typical values: ~0.07 for 10x v1, ~0.15 for v2."
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

    if use_memento_de and (
        not isinstance(memento_num_boot, int) or isinstance(memento_num_boot, bool)
    ):
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


def _require_explicit_groups(
    target_group: str | None,
    reference_group: str | None,
    *,
    func_name: str,
) -> None:
    """Require explicit contrast labels; prevents silent GA/Ctrl mismatches."""
    if target_group is None or reference_group is None:
        raise ValueError(
            f"{func_name}() requires explicit target_group and reference_group "
            "matching adata.obs[groupby] values. "
            'Historical defaults "GA"/"Ctrl" were removed to prevent silent mismatches '
            "that yield empty subsets and NaN results. "
            "Use active_score_simple() / differential_expression_simple() for "
            "Disease/Control convenience defaults, or call recommend_workflow() first."
        )


def _resolve_deprecated_active_score_kwargs(kwargs: dict[str, Any]) -> tuple[float, bool]:
    """Pop deprecated active_score kwargs and emit DeprecationWarning when used."""
    active_fdr_cutoff = 0.05
    prioritize_velocity = False
    if "active_fdr_cutoff" in kwargs:
        active_fdr_cutoff = kwargs.pop("active_fdr_cutoff")
        warnings.warn(
            "active_fdr_cutoff is deprecated and no longer used for the built-in "
            "'significant' gene list. Use unspliced_excess_fdr_cutoff instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        if not (0 < active_fdr_cutoff <= 1):
            raise ValueError("active_fdr_cutoff must be in (0, 1].")
    if "prioritize_velocity" in kwargs:
        prioritize_velocity = kwargs.pop("prioritize_velocity")
        warnings.warn(
            "prioritize_velocity is deprecated; use ranking_mode='nascent_excess' instead.",
            DeprecationWarning,
            stacklevel=3,
        )
    if kwargs:
        raise TypeError(
            f"active_score() got unexpected keyword argument(s): {', '.join(sorted(kwargs))}"
        )
    return active_fdr_cutoff, prioritize_velocity


def active_score(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str | None = None,
    reference_group: str | None = None,
    subset_col: str | None = None,
    subset_values: str | list[str] | tuple[str, ...] | None = None,
    weight_fc: float = 1.0,
    weight_unspliced: float = 1.0,
    weight_pval: float = 1.0,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.5,
    unspliced_excess_fdr_cutoff: float = 0.05,
    de_method: str = "t-test_overestim_var",  # freely switchable basic option, e.g. "wilcoxon"
    pseudobulk_de_backend: str = "pydeseq2",  # "pydeseq2" or "scanpy" when use_pseudobulk=True
    # Minimum summed counts (across obs) required to retain a gene when using the pydeseq2 backend.
    # The legacy hard-coded value of 10 can remove too many genes on tiny pseudobulk or low-depth data.
    pydeseq2_min_counts: int = 10,
    use_permutation: bool = False,
    perm_de_backend: str = "same",
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
    # Gamma estimation method for reference group unspliced/spliced ratio.
    # "heuristic_shrink": classic global-ratio + prior_weight shrinkage (default, prior_weight=5.0)
    # "robust_median": use median ratio from reference for better stability with small reference groups
    # "empirical_bayes": robust log-ratio empirical Bayes shrinkage (recommended for small reference)
    # "raw": minimal shrinkage (use observed ratios directly)
    gamma_method: str = "heuristic_shrink",
    ranking_mode: str = "composite",
    copy_input: bool = True,
    **deprecated_kwargs: Any,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Identify genes showing **higher** unspliced (nascent) RNA in the target group
    relative to reference (positive unspliced excess after reference-gamma correction),
    combined with upregulation (positive logFC).

    The returned ``significant`` DataFrame (second return value) is **strictly one-sided**:
    it only contains genes that are upregulated in target (logFC > cutoff) **AND**
    have positive bias-corrected unspliced excess (nascent excess > 0).
    Downregulated genes or genes with negative excess are never included in the
    built-in ``significant`` list even if DE is strong. Use ``all_results`` +
    :func:`filter_active_genes` (with ``logfc_direction="down"`` or ``"both"``) for
    other directions or custom thresholds.

    **Required:** ``target_group`` and ``reference_group`` must match values in
    ``adata.obs[groupby]`` (no implicit defaults). Use :func:`active_score_simple`
    for Disease/Control convenience defaults.

    **Recommended entry points (to avoid the long parameter list):**
    - For new users: :func:`active_score_simple`
    - For guided configuration: :func:`recommend_workflow` then ``active_score(..., **rec["suggested_kwargs"])``
    - Presets are defined in ``WORKFLOW_PRESETS``.

    Deprecated keyword-only arguments (emit ``DeprecationWarning``):
    ``active_fdr_cutoff``, ``prioritize_velocity``.

    The full signature below is for power users and internal composition. Many parameters
    have inter-dependencies that are validated early; hidden interactions exist (e.g.
    ranking_mode affects some weight_* defaults).

    The function computes:
    - logFC and p_adj between target and reference (via scanpy or PyDESeq2).
    - An unspliced (nascent) excess delta = U_target − (gamma_ref × S_target), where
      gamma_ref is a shrunk U/S ratio estimated in the reference group.
    - (by default) A Huber regression correction of the delta on log(gene length) and
      log(intron number); the residuals become ``unspliced_excess_residual``.
    - When ``use_permutation=True``, independent one-sided permutation p-values and
      BH-FDR are computed for the bias-corrected unspliced excess residual
      (``unspliced_excess_pval``, ``unspliced_excess_fdr``).
    - A soft-scaled, weighted combination of the three signals, scaled to 0–100.

    Several extensions are available as explicit options (see the README section
    "Optional advanced features"):
    - show_effective_gamma
    - gamma_method="robust_median" (heuristic variant of additive shrinkage using median per-gene ratio as base; not Bayesian)
    - gamma_method="empirical_bayes" (hierarchical empirical Bayes log-ratio shrinkage; recommended for small reference groups)
    - bias_correction="none"
    - use_mixed_model
    - use_permutation
    - prioritize_velocity (deprecated convenience; prefer ``ranking_mode="nascent_excess"``)
    - ranking_mode: ``"composite"`` (default) or ``"nascent_excess"`` (rank from
      unspliced_excess_residual only)

    Diagnostics (including global unspliced fraction and bias fit details) are stored
    under adata.uns["scatrans"]["diagnostics"]. The full ranked table (all_results)
    is the main output; the built-in significant list is produced by a strict
    conjunction of thresholds and is often small or empty.

    A separate function diagnose_design is available to summarize the experimental
    design and surface relevant warnings before analysis.

    **Important statistical note (reporting boundaries)**:
    - `active_score` is a **heuristic ranking score only**. It is NOT a p-value,
      effect size with calibrated uncertainty, or evidence of causal transcriptional activation.
    - `unspliced_excess_*` columns are **group-contrast proxies** (reference-gamma excess),
      not outputs of a stochastic/dynamical RNA velocity model. Do not treat them as
      literal nascent transcription rates without independent validation.
    - For significance claims use DE ``p_adj`` and/or permutation ``unspliced_excess_fdr``
      (when ``use_permutation=True``). Cross-check with orthogonal methods when possible.
    - The built-in ``significant`` list is a strict conjunction and is frequently empty —
      this is intentional. Use ``all_results`` + ``filter_active_genes`` for exploration.

    Full usage, recommended workflow, and result interpretation are documented in
    the package README ("Statistical interpretation and reporting boundaries").

    copy_input : bool, default True
        If True (default), deep-copy the input once after combining obs filters
        (``subset_col`` + target/reference groups) so the caller's object is not
        mutated. If False, reuse the input in-place when no obs filtering is
        required; otherwise subset without calling ``AnnData.copy()`` (lower memory
        on large objects). The returned AnnData is always the working object and
        may be mutated (new ``.var`` columns, layer remapping, etc.).
    """
    # ==================== EARLY VALIDATION (kept identical) ====================
    if mode not in {"heuristic", "advanced"}:
        raise ValueError("mode must be either 'heuristic' or 'advanced'")

    if gamma_method not in {"heuristic_shrink", "robust_median", "empirical_bayes", "raw"}:
        raise ValueError(
            "gamma_method must be one of "
            "{'heuristic_shrink', 'robust_median', 'empirical_bayes' (hierarchical), 'raw'}."
        )

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
    active_fdr_cutoff, prioritize_velocity = _resolve_deprecated_active_score_kwargs(
        deprecated_kwargs
    )
    _require_explicit_groups(target_group, reference_group, func_name="active_score")

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

    if not (0 < unspliced_excess_fdr_cutoff <= 1):
        raise ValueError("unspliced_excess_fdr_cutoff must be in (0, 1].")

    if min_total_counts < 0:
        raise ValueError("min_total_counts must be non-negative.")

    if perm_de_backend not in {"fast", "same"}:
        raise ValueError("perm_de_backend must be 'fast' or 'same'.")

    if ranking_mode not in {"composite", "nascent_excess"}:
        raise ValueError(
            "ranking_mode must be 'composite' (default) or 'nascent_excess' "
            "(active_score from unspliced_excess_residual only)."
        )

    default_weights = weight_fc == 1.0 and weight_unspliced == 1.0 and weight_pval == 1.0
    if ranking_mode == "nascent_excess":
        if prioritize_velocity:
            logger.warning("prioritize_velocity is ignored when ranking_mode='nascent_excess'.")
        if default_weights:
            weight_fc, weight_pval, weight_unspliced = 0.0, 0.0, 1.0
            logger.info(
                "ranking_mode='nascent_excess': active_score ranks genes by "
                "bias-corrected unspliced_excess_residual only (DE weights set to 0)."
            )
    elif prioritize_velocity and default_weights:
        weight_unspliced = 3.0
        weight_fc = 0.5
        weight_pval = 0.5
        logger.info(
            "prioritize_velocity=True (deprecated): emphasizing nascent excess in composite "
            "active_score. Prefer ranking_mode='nascent_excess' for pure residual ranking."
        )

    if weight_fc + weight_unspliced + weight_pval <= 0:
        raise ValueError(
            "At least one of weight_fc, weight_unspliced, weight_pval must be positive."
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
    obs_filter = pd.Series(True, index=adata_input.obs_names)
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
        n_after = int(subset_mask.sum())
        if n_after == 0:
            raise ValueError(f"No cells remain after subsetting {subset_col}")
        obs_filter &= subset_mask
        logger.info("Subsetted by %s (%d/%d cells remaining)", subset_col, n_after, n_before)

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")

    target_group, reference_group, norm_groups = _validate_group_contrast(
        adata_input.obs[groupby],
        groupby=groupby,
        target_group=str(target_group),
        reference_group=str(reference_group),
    )

    obs_filter &= norm_groups.isin([target_group, reference_group])

    # Automatic design guidance for small-sample or replicate-structured data
    if sample_col or use_pseudobulk:
        try:
            _ = diagnose_design(
                adata_input,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                sample_col=sample_col,
                copy_input=False,  # pure read-only diagnostic; avoid expensive deep copy
            )
        except Exception as e:
            logger.debug("diagnose_design skipped (non-fatal): %s", e)

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

    adata = _select_obs(adata_input, obs_filter, copy_input=copy_input)
    if adata.n_obs == 0:
        raise ValueError(
            "No cells match target/reference groups after filtering. "
            f"Check target_group='{target_group}' and reference_group='{reference_group}' "
            f"against adata.obs['{groupby}'] (missing labels are excluded)."
        )
    adata = _materialize_if_view(adata)
    adata.obs[groupby] = norm_groups.loc[obs_filter].values

    # Perform layer remapping on the working adata (copy_input=True isolates caller's object).
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
        adata = _select_var(
            adata, adata.var["gene_type"] == gene_type_filter, copy_input=copy_input
        )
        adata = _materialize_if_view(adata)

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
            layers=["spliced", "unspliced"],
            x_layer=pb_x_layer,
            use_total_for_x=pb_use_total_for_x,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        is_pseudobulk = True
        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", category=UserWarning
            )  # pandas/ann implicit index str conversion is benign here
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
    _apply_de_preprocess(
        adata,
        de_preprocess,
        skip_auto=is_pseudobulk and pseudobulk_de_backend == "pydeseq2",
    )

    # For permutation tasks we pass *already preprocessed* adata copies (or layers).
    # Re-applying normalize_log1p (or auto) inside _single_permutation_task would double-transform
    # the permuted data while the real data was transformed only once -> biased FDR.
    # "auto" is safe in practice because of the "log1p" in uns check, but we force "none" for perm.
    perm_de_preprocess = "none"

    # Always provide feature matrix for genes with length/intron info when available.
    # The actual decision to use Huber regression (vs median fallback) lives inside
    # _fit_huber_bias_correction and is based on n_fit >= min_fit_obs (30 by default).
    # This removes the prior 50/30 threshold mismatch.
    X_features = (
        np.column_stack([np.log1p(gene_length[valid_feat]), np.log1p(intron_number[valid_feat])])
        if np.any(valid_feat)
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

    # Auto-resolve preserved raw counts for count-based DE backends (alignment-checked).
    resolved_counts = None
    needs_raw_counts = use_memento_de or (is_pseudobulk and pseudobulk_de_backend == "pydeseq2")
    if needs_raw_counts:
        resolved_counts = _resolve_aligned_raw_counts(adata, layer="counts", require_integer=True)
        if resolved_counts is None and ("counts" in adata.layers or getattr(adata, "raw", None)):
            logger.warning(
                "Count-based DE was requested but no safely aligned raw counts were found. "
                "Call scat.store_raw_counts(adata) before HVG/normalize, or pass counts= explicitly."
            )

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
        min_counts_per_gene=pydeseq2_min_counts,
    )

    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]
    if "delta_variance" in de_df.columns:
        adata.var["delta_variance"] = de_df["delta_variance"]
    if "delta_var_pval" in de_df.columns:
        adata.var["delta_var_pval"] = de_df["delta_var_pval"]

    # Surface Memento-specific columns when the memento backend was used (for variability etc.)
    for extra_col in [
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
        "memento_p_adj_native",
    ]:
        if extra_col in de_df.columns:
            adata.var[extra_col] = de_df[extra_col]

    n_mixed_failed = 0
    mixed_failed_rate = 0.0
    if use_mixed_model:
        n_mixed_failed = int(
            de_df.attrs.get("n_genes_failed_fit", 0) if hasattr(de_df, "attrs") else 0
        )
        mixed_failed_rate = float(
            de_df.attrs.get("failed_fit_rate", 0.0) if hasattr(de_df, "attrs") else 0.0
        )

    # ==================== QC: global unspliced fraction (integrated high-value diagnostic) ====================
    unspliced_fraction = np.nan
    if "spliced" in adata.layers and "unspliced" in adata.layers:
        try:
            unspliced_fraction = _qc.unspliced_global(
                adata, spliced_key="spliced", unspliced_key="unspliced", warn_threshold=0.5
            )
        except Exception as _e:
            logger.debug("Could not compute global unspliced fraction: %s", _e)
    else:
        logger.debug(
            "Skipping global unspliced fraction: required layers not present after filtering."
        )

    uns_layer_raw = adata.layers["unspliced"]
    spl_layer_raw = adata.layers["spliced"]

    # Always library-size normalize spliced + unspliced layers (per cell or per pseudobulk sample)
    # before computing group means for velocity delta. This removes per-observation depth
    # confounding from the U/S excess statistic.
    # Previously only done for pseudobulk; cell-level path inherited raw count scale differences
    # between groups (common batch/sequencing-depth effects).
    uns_layer, spl_layer, _row_totals, _factors = _normalize_velocity_layers_by_size_factor(
        uns_layer_raw, spl_layer_raw
    )

    obs_labels = adata.obs[groupby].map(_normalize_group_label).values
    t_mask = obs_labels == target_group
    r_mask = obs_labels == reference_group

    # ==================== VELOCITY DELTA (dual track) ====================
    moments_info: dict[str, Any] = {}
    velocity_layer_for_perm_uns = uns_layer
    velocity_layer_for_perm_spl = spl_layer
    gamma_ref = np.full(adata.n_vars, np.nan)  # will be overwritten in all branches
    gamma_info: dict[str, Any] = {}

    if mode == "heuristic":
        delta_velocity, total_us_velocity, gamma_ref, gamma_info = _compute_velocity_delta(
            uns_layer, spl_layer, t_mask, r_mask, prior_weight, gamma_method=gamma_method
        )
        velocity_source = "heuristic_global_ratio"

    elif mode == "advanced":
        adata_comp = adata.copy()
        # Use the (library-size normalized) velocity layers for moments computation too.
        # This ensures cell-level advanced mode also benefits from depth correction
        # before neighbor graph + scVelo moments.
        adata_comp.layers["unspliced"] = uns_layer.copy()
        adata_comp.layers["spliced"] = spl_layer.copy()

        try:
            delta_velocity, total_us_velocity, gamma_ref, moments_info = (
                _compute_moments_velocity_delta(
                    adata_comp,
                    t_mask,
                    r_mask,
                    prior_weight=prior_weight,
                    gamma_method=gamma_method,
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
            gamma_info = moments_info.get("gamma_info", {})
        except (ValueError, RuntimeError, KeyError, TypeError, AttributeError) as e:
            if advanced_fallback:
                import traceback

                tb = traceback.format_exc()
                logger.warning(
                    "Advanced mode failed (%s). Falling back to heuristic. "
                    "Set advanced_fallback=False (or use mode='heuristic') to see the full traceback.",
                    e,
                )
                logger.debug("Advanced mode full traceback:\n%s", tb)
                delta_velocity, total_us_velocity, gamma_ref, gamma_info = _compute_velocity_delta(
                    uns_layer, spl_layer, t_mask, r_mask, prior_weight, gamma_method=gamma_method
                )
                velocity_source = "heuristic_fallback_from_advanced"
                moments_info = {
                    "advanced_failed": True,
                    "failure_reason": str(e),
                    "traceback": tb,
                }
            else:
                raise
        except Exception as e:
            # For completely unexpected errors, still respect fallback but be more explicit
            if advanced_fallback:
                import traceback

                tb = traceback.format_exc()
                logger.warning(
                    "Advanced mode hit an unexpected error. Falling back. "
                    "Please report this if it persists. Full traceback is in diagnostics.",
                )
                logger.debug("Unexpected advanced error:\n%s", tb)
                delta_velocity, total_us_velocity, gamma_ref, gamma_info = _compute_velocity_delta(
                    uns_layer, spl_layer, t_mask, r_mask, prior_weight, gamma_method=gamma_method
                )
                velocity_source = "heuristic_fallback_from_advanced"
                moments_info = {"advanced_failed": True, "failure_reason": str(e), "traceback": tb}
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

    _write_unspliced_excess_columns(adata.var, delta=delta_velocity, residual=residual)
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

    if gamma_method == "empirical_bayes" and "shrinkage_weights" in gamma_info:
        adata.var["gamma_shrinkage_weight"] = gamma_info["shrinkage_weights"]

    eb_prior_for_perm: dict[str, Any] | None = None
    if gamma_method == "empirical_bayes":
        eb_prior_for_perm = gamma_info.get("eb_prior")

    # ==================== DIAGNOSTICS (high priority for usability & paper rigor) ====================
    gamma_stats = gamma_info.get("effective_gamma_stats")
    if not gamma_stats:
        try:
            gamma_stats = {
                "median": float(np.nanmedian(gamma_ref))
                if np.any(np.isfinite(gamma_ref))
                else np.nan,
                "mean": float(np.nanmean(gamma_ref)) if np.any(np.isfinite(gamma_ref)) else np.nan,
                "min": float(np.nanmin(gamma_ref)) if np.any(np.isfinite(gamma_ref)) else np.nan,
                "max": float(np.nanmax(gamma_ref)) if np.any(np.isfinite(gamma_ref)) else np.nan,
                "n_finite": int(np.isfinite(gamma_ref).sum()),
            }
        except Exception:
            gamma_stats = {"median": np.nan, "mean": np.nan, "min": np.nan, "max": np.nan}

    velocity_diag: dict[str, Any] = {
        "source": velocity_source,
        "n_genes_with_finite_delta": int(np.isfinite(delta_velocity).sum()),
        "effective_gamma_exposed": bool(show_effective_gamma),
        "gamma_method": gamma_method,
        "prior_weight": float(prior_weight),
        "effective_gamma_stats": gamma_stats,
        "gamma_method_detailed": gamma_info.get("gamma_method_detailed"),
        "shrinkage_note": (
            "per-gene shrinkage applied using prior_weight (higher = stronger pull toward reference ratio)"
            if gamma_method in ("heuristic_shrink", "robust_median")
            else (
                "log-ratio empirical Bayes shrinkage with robust trimmed prior"
                if gamma_method == "empirical_bayes"
                else "minimal/no shrinkage"
            )
        ),
    }
    for _gk in (
        "gamma_prior_mean",
        "gamma_prior_scale",
        "gamma_prior_tau_squared",
        "n_genes_used_for_prior",
        "shrinkage_summary",
        "fallback_triggered",
        "count_pseudocount",
        "used_fixed_prior",
    ):
        if _gk in gamma_info:
            velocity_diag[_gk] = gamma_info[_gk]

    diagnostics: dict[str, Any] = {
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "n_genes_with_valid_features": int(valid_feat.sum()),
        "unspliced_global_fraction": float(unspliced_fraction)
        if unspliced_fraction is not None
        else np.nan,
        "bias_correction": bias_info,
        "velocity": velocity_diag,
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
            "n_genes_failed_fit": n_mixed_failed if use_mixed_model else 0,
            "failed_fit_rate": mixed_failed_rate if use_mixed_model else 0.0,
            "note": (
                "Lightweight LMM analogue (log1p + Wald/LRT); not NB-GLMM/voom. "
                "Inspect failed_fit_rate before publication claims."
                if use_mixed_model
                else None
            ),
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
    _perm_disabled_reason = None
    # disabled_reason (e.g. "small_permutation_space") is now returned from run_permutation_test
    # and propagated to metadata/diagnostics for transparency.

    if use_permutation:
        if is_pseudobulk:
            n_t, n_r = t_mask.sum(), r_mask.sum()
            current_max_perm = float("inf") if n_t + n_r > 30 else max(1, comb(n_t + n_r, n_t) - 1)

        if perm_de_backend == "fast":
            perm_pb_backend, perm_de_method = "scanpy", "t-test_overestim_var"
            logger.warning(
                "perm_de_backend='fast': permutations use scanpy t-test regardless of the "
                "main analysis DE backend (%s). Null and observed statistics may be "
                "mismatched — prefer the default perm_de_backend='same' for reporting.",
                de_method,
            )
        elif perm_de_backend == "same":
            perm_pb_backend, perm_de_method = pseudobulk_de_backend, de_method
            logger.info(
                "perm_de_backend='same': permutations use the same DE backend as the "
                "main analysis (%s / pseudobulk=%s).",
                de_method,
                is_pseudobulk,
            )
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

        # Delegate to the canonical implementation in _permutation (eliminates dead-code duplication).
        # Pass explicit valid_expr to avoid .var.get differences and ensure consistent masking.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            (
                active_score_pval_arr,
                active_score_fdr_arr,
                unspliced_excess_pval_arr,
                unspliced_excess_fdr_arr,
                perm_use_fdr,
                _perm_disabled_reason,
            ) = run_permutation_test(
                n_perm=n_perm,
                effective_n_jobs=effective_n_jobs,
                random_seed=random_seed,
                obs_labels=obs_labels,
                target_group=target_group,
                reference_group=reference_group,
                adata=adata,
                X_features=X_features,
                valid_feat=valid_feat,
                velocity_layer_for_perm_uns=velocity_layer_for_perm_uns,
                velocity_layer_for_perm_spl=velocity_layer_for_perm_spl,
                total_us_raw=total_us_raw,
                min_total_counts=min_total_counts,
                weight_fc=weight_fc,
                weight_unspliced=weight_unspliced,
                weight_pval=weight_pval,
                lambda_fc=lambda_fc,
                lambda_res=lambda_res,
                lambda_pval=lambda_pval,
                is_pseudobulk=is_pseudobulk,
                perm_pb_backend=perm_pb_backend,
                perm_de_method=perm_de_method,
                prior_weight=prior_weight,
                gamma_method=gamma_method,
                de_preprocess=perm_de_preprocess,
                strict_pydeseq2_counts=strict_pydeseq2_counts,
                real_score=real_score,
                real_residual=residual,
                eb_prior=eb_prior_for_perm,
                velocity_source=velocity_source,
                bias_correction=bias_correction,
                use_memento_de=perm_memento_de,
                memento_capture_rate=memento_capture_rate,
                memento_num_boot=memento_num_boot,
                memento_n_cpus=memento_n_cpus,
                valid_expr=valid_expr,
            )

        adata.var["active_score_pval"] = active_score_pval_arr
        adata.var[UNSPLICED_EXCESS_PVAL_COL] = unspliced_excess_pval_arr
        adata.var["active_score_fdr"] = active_score_fdr_arr
        adata.var[UNSPLICED_EXCESS_FDR_COL] = unspliced_excess_fdr_arr

        # Use the decision returned by canonical run_permutation_test (now authoritative).
        # Previously always-True + vestigial if-not + separate max_perm check were inconsistent.
        use_fdr_for_significance = bool(perm_use_fdr)

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
        unspliced_excess_fdr_cutoff=unspliced_excess_fdr_cutoff,
        use_fdr_for_significance=use_fdr_for_significance,
        perm_disabled_reason=_perm_disabled_reason,
        use_delta_variance_pval=use_delta_variance_pval,
        delta_var_pval_cutoff=delta_var_pval_cutoff,
        de_method=de_method,
        pseudobulk_de_backend=pseudobulk_de_backend,
        perm_de_backend=perm_de_backend if use_permutation else None,
        use_memento_de=use_memento_de,
        perm_use_memento_de=perm_use_memento_de if use_permutation else None,
        memento_capture_rate=memento_capture_rate if use_memento_de else None,
        prior_weight=prior_weight,
        gamma_method=gamma_method,
        min_total_counts=min_total_counts,
        random_seed=random_seed,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        prioritize_velocity=prioritize_velocity,
        ranking_mode=ranking_mode,
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
    # Build result columns. Primary unspliced-excess names; legacy velocity_* remain in adata.var.
    cols = [
        "active_score",
        UNSPLICED_EXCESS_DELTA_COL,
        UNSPLICED_EXCESS_RESIDUAL_COL,
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
        cols.extend(
            [
                "active_score_pval",
                "active_score_fdr",
                UNSPLICED_EXCESS_PVAL_COL,
                UNSPLICED_EXCESS_FDR_COL,
            ]
        )
    if "delta_variance" in adata.var.columns:
        cols.append("delta_variance")
    if "delta_var_pval" in adata.var.columns:
        cols.append("delta_var_pval")
    for mc in ["memento_de_se", "memento_dv_coef", "memento_dv_se", "memento_dv_pval"]:
        if mc in adata.var.columns and mc not in cols:
            cols.append(mc)
    if "gamma_shrinkage_weight" in adata.var.columns:
        cols.append("gamma_shrinkage_weight")
    cols = [c for c in cols if c in adata.var.columns]

    # Built-in significant list: DE significance + positive unspliced excess + permutation FDR.
    # active_score is for ranking/visualization only (heuristic, not a p-value).
    residual_col = (
        UNSPLICED_EXCESS_RESIDUAL_COL
        if UNSPLICED_EXCESS_RESIDUAL_COL in adata.var.columns
        else LEGACY_VELOCITY_RESIDUAL_COL
    )
    ue_fdr_cutoff = extra_metadata.get("unspliced_excess_fdr_cutoff", 0.05)

    if use_permutation and UNSPLICED_EXCESS_FDR_COL in adata.var.columns:
        mask = (
            (adata.var["p_adj"] < extra_metadata.get("pval_cutoff", 0.05))
            & (adata.var["logFC"] > extra_metadata.get("logfc_cutoff", 0.5))
            & (adata.var[residual_col] > 0)
            & (adata.var["valid_expr"])
        )
        if extra_metadata.get("use_fdr_for_significance", True):
            mask = mask & (adata.var[UNSPLICED_EXCESS_FDR_COL] < ue_fdr_cutoff)
        else:
            logger.warning(
                "Permutation space is very small; unspliced_excess_fdr was not applied "
                "to the built-in significant list."
            )
        # Apply delta_variance pval filter only in the permutation path where a meaningful
        # significant mask can be produced (avoids dead-code application to the all-False mask).
        if extra_metadata.get("use_delta_variance_pval") and "delta_var_pval" in adata.var.columns:
            mask = mask & (
                adata.var["delta_var_pval"] < extra_metadata.get("delta_var_pval_cutoff", 0.05)
            )
    else:
        mask = pd.Series(False, index=adata.var.index)
        if not use_permutation:
            logger.warning(
                "Built-in 'significant' list requires use_permutation=True "
                "(unspliced_excess_fdr_cutoff=%.3f). Returning empty significant; "
                "inspect all_results and use filter_active_genes for custom thresholds.",
                ue_fdr_cutoff,
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
    # Keep a lightweight history so multiple calls don't completely overwrite previous runs
    history = existing.get("history", [])
    if "analysis" in existing:
        # save previous run summary (lightweight)
        prev = {
            k: existing.get(k)
            for k in ("analysis", "mode", "target_group", "reference_group", "timestamp")
            if k in existing
        }
        if prev:
            history.append(prev)
            if len(history) > 5:
                history = history[-5:]
    existing["history"] = history

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
        "unspliced_excess_fdr_cutoff": extra_metadata.get("unspliced_excess_fdr_cutoff"),
        "use_fdr_for_significance": extra_metadata.get("use_fdr_for_significance"),
        "perm_disabled_reason": extra_metadata.get("perm_disabled_reason"),
        "significant_criteria": {
            "logFC": f"> {extra_metadata.get('logfc_cutoff')}",
            "p_adj": f"< {extra_metadata.get('pval_cutoff')}",
            "unspliced_excess_residual": "> 0",
            "unspliced_excess_fdr": (
                f"< {ue_fdr_cutoff} (requires use_permutation=True)"
                if use_permutation
                else "not evaluated (use_permutation=False)"
            ),
            "active_score": "ranking only (not used for significance)",
        },
        "use_delta_variance_pval": extra_metadata.get("use_delta_variance_pval"),
        "delta_var_pval_cutoff": extra_metadata.get("delta_var_pval_cutoff"),
        "de_method": extra_metadata.get("de_method"),
        "use_pseudobulk": extra_metadata.get("is_pseudobulk", False),
        "pseudobulk_de_backend": extra_metadata.get("pseudobulk_de_backend"),
        "perm_de_backend": extra_metadata.get("perm_de_backend"),
        "use_memento_de": extra_metadata.get("use_memento_de"),
        "perm_use_memento_de": extra_metadata.get("perm_use_memento_de"),
        "memento_capture_rate": extra_metadata.get("memento_capture_rate"),
        "prior_weight": extra_metadata.get("prior_weight"),
        "gamma_method": extra_metadata.get("gamma_method", "heuristic_shrink"),
        "min_total_counts": extra_metadata.get("min_total_counts"),
        "random_seed": extra_metadata.get("random_seed"),
        "use_mixed_model": extra_metadata.get("use_mixed_model"),
        "sample_col": extra_metadata.get("sample_col"),
        "prioritize_velocity": extra_metadata.get("prioritize_velocity"),
        "ranking_mode": extra_metadata.get("ranking_mode", "composite"),
        "weight_fc": extra_metadata.get("weight_fc"),
        "weight_unspliced": extra_metadata.get("weight_unspliced"),
        "weight_pval": extra_metadata.get("weight_pval"),
        "diagnostics": diagnostics,
        "unspliced_global_fraction": diagnostics.get("unspliced_global_fraction", np.nan),
    }

    if use_permutation:
        note = (
            "For efficiency, unspliced/spliced layers and reference gamma are fixed from the "
            "original data. Group labels are shuffled to recompute DE, unspliced excess residual, "
            "composite active_score, and unspliced_excess permutation p-values."
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
    target_group: str | None = None,
    reference_group: str | None = None,
    subset_col: str | None = None,
    subset_values: str | list[str] | tuple[str, ...] | None = None,
    de_method: str = "t-test_overestim_var",
    pseudobulk_de_backend: str = "pydeseq2",
    pydeseq2_min_counts: int = 10,
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
    n_jobs: int = -1,
    gene_type_filter: str | None = None,
    # Allow providing raw counts separately when adata.X is already HVG+log (very common)
    counts: str | np.ndarray | sparse.spmatrix | pd.DataFrame | ad.AnnData | None = None,
    copy_input: bool = True,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Standalone differential expression (DE) using the same flexible backends
    as scATrans (scanpy methods, PyDESeq2 pseudobulk, mixed linear models,
    and Memento -- the Cell 2024 method-of-moments framework).

    This function does **not** require spliced/unspliced (velocity) layers.
    It is intended for users who want high-quality DE (especially via Memento),
    followed by scATrans' downstream tools:

        candidates = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3)  # upregulated
        # down or both directions:
        # down_cands = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="down")
        # For enrichment, pass adata= (if store_raw_counts was used) so it uses
        # the preserved full measured gene set as universe, not just current HVGs.
        enrich = scat.run_enrichment(candidates.index.tolist(), ..., adata=adata)
        scat.pl.volcano_plot(de_results, ...)
        scat.pl.enrich_dotplot(enrich, ...)

    All DE-related options from `active_score` are supported here
    (pseudobulk, mixed models, Memento, etc.), except permutation-based FDR
    (use ``active_score(..., use_permutation=True)`` when velocity layers are available).
    For a minimal-parameter entry point see ``active_score_simple`` or ``run_default_pipeline``.

    copy_input : bool, default True
        Same semantics as :func:`active_score`: one combined obs-filter copy when
        True; zero ``AnnData.copy()`` calls when False and no obs filtering is needed.

    Returns
    -------
    (adata_with_results, results_df)
        - results_df is a ranked DataFrame (by |logFC| or p_adj) containing
          at minimum: logFC, p_val, p_adj, and (when use_memento_de) the
          native memento_de_* / memento_dv_* columns.
        - adata.var is updated with the same columns for convenience.
        - Metadata is stored under adata.uns["scatrans"].
    """
    _require_explicit_groups(target_group, reference_group, func_name="differential_expression")

    # --- minimal shared validation (subset + group checks) ---
    obs_filter = pd.Series(True, index=adata_input.obs_names)
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
        if int(subset_mask.sum()) == 0:
            raise ValueError("No cells remain after subsetting.")
        obs_filter &= subset_mask

    if not adata_input.var_names.is_unique:
        raise ValueError("adata.var_names must be unique.")

    if groupby not in adata_input.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found.")

    target_group, reference_group, norm_groups = _validate_group_contrast(
        adata_input.obs[groupby],
        groupby=groupby,
        target_group=str(target_group),
        reference_group=str(reference_group),
    )

    obs_filter &= norm_groups.isin([target_group, reference_group])
    adata_input = _select_obs(adata_input, obs_filter, copy_input=copy_input)
    if adata_input.n_obs == 0:
        raise ValueError(
            "No cells match target/reference groups after filtering. "
            f"Check target_group='{target_group}' and reference_group='{reference_group}' "
            f"against adata.obs['{groupby}'] (missing labels are excluded)."
        )
    adata_input = _materialize_if_view(adata_input)
    adata_input.obs[groupby] = norm_groups.loc[obs_filter].values

    if gene_type_filter:
        if "gene_type" not in adata_input.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata_input = _select_var(
            adata_input, adata_input.var["gene_type"] == gene_type_filter, copy_input=copy_input
        )
        adata_input = _materialize_if_view(adata_input)

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
        use_permutation=False,
        n_perm=0,
        use_mixed_model=use_mixed_model,
        mixed_model_pval=mixed_model_pval,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        min_cells=min_cells,
        min_counts=min_counts,
    )

    if min_total_counts != 50:
        logger.warning(
            "differential_expression: min_total_counts=%s is not enforced in the DE-only path "
            "(reserved for API compatibility). It affects gene filtering in active_score() only. "
            "Use min_cells / min_counts for pseudobulk filtering instead.",
            min_total_counts,
        )

    if use_delta_variance_pval:
        logger.warning(
            "differential_expression: use_delta_variance_pval=True is not enforced in the DE-only "
            "path (this function returns the full ranked results table, not a significant-gene "
            "subset). Use active_score() for delta-variance filtering, or filter manually via "
            "results['delta_var_pval'] < delta_var_pval_cutoff (currently %.4g).",
            delta_var_pval_cutoff,
        )

    # Auto-resolve aligned raw counts for count-based backends (Memento / PyDESeq2).
    if counts is None and (
        use_memento_de or (use_pseudobulk and pseudobulk_de_backend == "pydeseq2")
    ):
        counts = _resolve_aligned_raw_counts(adata_input, layer="counts", require_integer=True)
        if counts is None and ("counts" in adata_input.layers or getattr(adata_input, "raw", None)):
            logger.warning(
                "Count-based DE was requested but no safely aligned raw counts were found. "
                "Call scat.store_raw_counts(adata) before HVG/normalize, or pass counts= explicitly."
            )

    # --- prepare data (pseudobulk if requested) ---
    adata = adata_input

    if use_pseudobulk:
        logger.info("Performing pseudobulk aggregation for DE...")
        available_layers = [layer for layer in ("spliced", "unspliced") if layer in adata.layers]
        adata = _pseudobulk_with_layers(
            adata,
            sample_col,
            groupby,
            layers=available_layers,
            x_layer=pb_x_layer if pb_x_layer != "X" else None,
            use_total_for_x=pb_use_total_for_x and bool(available_layers),
            min_cells=min_cells,
            min_counts=min_counts,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            adata.obs[groupby] = pd.Categorical(
                adata.obs[groupby].astype(str), categories=[reference_group, target_group]
            )

    # DE preprocess
    # (Memento coercion to 'none' already performed early via _coerce_memento_de_preprocess)

    _apply_de_preprocess(
        adata,
        de_preprocess,
        skip_auto=use_pseudobulk and pseudobulk_de_backend == "pydeseq2",
    )

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
        min_counts_per_gene=pydeseq2_min_counts,
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
        "memento_p_adj_native",
    ]:
        if extra in de_df.columns:
            adata.var[extra] = de_df[extra]

    n_mixed_failed_de = (
        int(de_df.attrs.get("n_genes_failed_fit", 0))
        if (use_mixed_model and hasattr(de_df, "attrs"))
        else 0
    )
    mixed_failed_rate_de = (
        float(de_df.attrs.get("failed_fit_rate", 0.0))
        if (use_mixed_model and hasattr(de_df, "attrs"))
        else 0.0
    )

    # Build clean results table (no velocity columns)
    cols = ["logFC", "p_val", "p_adj"]
    for c in [
        "delta_variance",
        "delta_var_pval",
        "memento_de_se",
        "memento_dv_coef",
        "memento_dv_se",
        "memento_dv_pval",
        "memento_p_adj_native",
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

    # Metadata — merge to preserve raw_gene_list etc. from store_raw_counts()
    existing = dict(adata.uns.get("scatrans", {}))
    history = existing.get("history", [])
    if "analysis" in existing:
        prev = {
            k: existing.get(k)
            for k in ("analysis", "mode", "target_group", "reference_group")
            if k in existing
        }
        if prev:
            history.append(prev)
            if len(history) > 5:
                history = history[-5:]
    existing["history"] = history

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
            "n_genes_failed_mixed_fit": n_mixed_failed_de,
            "failed_fit_rate_mixed": mixed_failed_rate_de,
            "memento_has_native_padj": bool(
                use_memento_de and "memento_p_adj_native" in de_df.columns
            ),
            "memento_num_boot": memento_num_boot if use_memento_de else None,
            "memento_n_cpus": memento_n_cpus if use_memento_de else None,
            "n_jobs": n_jobs,
            "gene_type_filter": gene_type_filter,
        }
    )
    adata.uns["scatrans"] = existing

    logger.info("DE completed. %d genes in results table.", len(results))
    return adata, results


def _record_raw_counts_metadata(
    adata: Any, *, save_raw: bool = False, overwrite: bool = False
) -> None:
    """Write scatrans metadata, optional adata.raw, and raw_* velocity layers."""
    if "scatrans" not in adata.uns:
        adata.uns["scatrans"] = {}
    prev_list = adata.uns["scatrans"].get("raw_gene_list")
    n_genes = int(adata.n_vars)
    if prev_list is not None and len(prev_list) != n_genes:
        logger.warning(
            "Updating raw_gene_list (%d → %d genes). If you subsetted after an earlier "
            "store_raw_counts(), enrichment universe will reflect the current gene set only. "
            "For the full pre-HVG universe, keep a separate full-gene AnnData or use save_raw=True "
            "on the original object before subsetting.",
            len(prev_list),
            n_genes,
        )
    adata.uns["scatrans"]["raw_gene_list"] = list(adata.var_names)
    adata.uns["scatrans"]["store_raw_counts_n_genes"] = n_genes
    logger.info(
        "Saved the current gene list as the measured universe for enrichment "
        "(in adata.uns['scatrans']['raw_gene_list'])."
    )

    if save_raw:
        if getattr(adata, "raw", None) is not None and not overwrite:
            logger.debug("adata.raw already exists; skipping (pass overwrite=True to replace).")
        else:
            adata.raw = adata.copy()
            logger.info("Set adata.raw to preserve full data.")

    for vel_name in ("spliced", "unspliced", "mature", "nascent"):
        if vel_name in adata.layers:
            raw_vel_name = f"raw_{vel_name}"
            if raw_vel_name not in adata.layers or overwrite:
                adata.layers[raw_vel_name] = adata.layers[vel_name].copy()
                logger.info(f"Saved original {vel_name} to adata.layers['{raw_vel_name}'].")


def ensure_raw_counts(
    adata: Any, layer: str = "counts", save_raw: bool = False, overwrite: bool = False
) -> None:
    """
    Ensure raw integer counts are available in ``adata.layers[layer]``.

    Convenience wrapper around :func:`store_raw_counts` that also tries to recover
    counts from ``adata.raw`` when ``adata.X`` is already normalized or log-transformed
    (common after HVG + ``sc.pp.log1p``).

    Resolution order:
    1. Existing ``layers[layer]`` if it already looks like integer counts
    2. Current ``adata.X`` if integer counts-like
    3. ``adata.raw.X`` when gene names/order match ``adata.var_names``

    Always updates ``adata.uns['scatrans']['raw_gene_list']`` and velocity ``raw_*`` layers
    via the same metadata path as :func:`store_raw_counts`.
    """
    if (
        layer in adata.layers
        and not overwrite
        and _is_integer_counts_like(adata.layers[layer])
        and adata.layers[layer].shape[1] == adata.n_vars
    ):
        _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)
        logger.debug("Layer '%s' already holds aligned integer counts.", layer)
        return

    if _is_integer_counts_like(adata.X):
        store_raw_counts(adata, layer=layer, save_raw=save_raw, overwrite=overwrite)
        return

    raw = getattr(adata, "raw", None)
    if raw is not None and _is_integer_counts_like(raw.X) and raw.shape[1] == adata.n_vars:
        if hasattr(raw, "var_names") and np.array_equal(raw.var_names, adata.var_names):
            adata.layers[layer] = raw.X.copy()
            logger.info(
                "ensure_raw_counts: recovered raw counts from adata.raw into layers['%s'].",
                layer,
            )
            _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)
            return
        logger.warning(
            "adata.raw exists but gene names/order do not match current adata.var_names. "
            "Cannot recover counts automatically."
        )

    logger.warning(
        "ensure_raw_counts: adata.X does not look like raw counts and adata.raw could not be used. "
        "Falling back to store_raw_counts (may warn again)."
    )
    store_raw_counts(adata, layer=layer, save_raw=save_raw, overwrite=overwrite)


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
        mat = adata.X.copy()
        if mat.shape[1] != adata.n_vars:
            raise ValueError(
                f"Cannot store raw counts: matrix has {mat.shape[1]} columns "
                f"but adata has {adata.n_vars} genes."
            )
        adata.layers[layer] = mat
        logger.info(f"Saved raw counts to adata.layers['{layer}'].")

    if layer in adata.layers and adata.layers[layer].shape[1] != adata.n_vars:
        raise ValueError(
            f"Layer '{layer}' has {adata.layers[layer].shape[1]} columns but adata has "
            f"{adata.n_vars} genes. Pass overwrite=True after fixing alignment."
        )

    _record_raw_counts_metadata(adata, save_raw=save_raw, overwrite=overwrite)


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

    target = adata if inplace else adata.copy()
    target.X = raw
    _clear_log_preprocess_metadata(target)
    if inplace:
        logger.info(f"Restored raw counts from {source} into adata.X (inplace).")
        return None
    logger.info(f"Created copy with raw counts from {source} in .X.")
    return target


_NOT_PROVIDED = object()


def filter_active_genes(
    results: pd.DataFrame,
    *,
    preset: str | None = None,
    active_score_cutoff: Any = _NOT_PROVIDED,
    pval_cutoff: Any = _NOT_PROVIDED,
    unspliced_excess_residual_cutoff: Any = _NOT_PROVIDED,
    logfc_cutoff: Any = _NOT_PROVIDED,
    logfc_direction: str = "up",
    active_score_fdr_cutoff: Any = _NOT_PROVIDED,
    unspliced_excess_fdr_cutoff: Any = _NOT_PROVIDED,
    effective_gamma_min: Any = _NOT_PROVIDED,
    effective_gamma_max: Any = _NOT_PROVIDED,
    delta_variance_min: Any = _NOT_PROVIDED,
    return_mask: bool = False,
    inplace: bool = False,
    **deprecated_kwargs: Any,
) -> pd.DataFrame | pd.Series:
    """Apply custom post-filtering to a results DataFrame (from `active_score` or `differential_expression`).

    This helper works for both:
    - Full `active_score` output (has `active_score` + unspliced_excess_residual).
    - Pure DE results from `differential_expression` (only logFC / p_adj + optional memento columns).

    It standardizes the common workflow:
    1. Run `active_score(...)` or `differential_expression(...)`.
    2. Use this function on the returned table to derive a final gene list.

    The function supports `preset` to automatically select reasonable default thresholds
    for different analysis modes:

    - preset="heuristic": stricter defaults suitable for typical single-cell data with default weights
      (active_score >= 55, unspliced_excess_residual > 1.0, etc.).
    - preset="pseudobulk": more lenient defaults that account for the much smaller
      magnitude of unspliced_excess_residual and active_score after sample-level aggregation
      (active_score >= 5, unspliced_excess_residual > 0.05, logFC > 0.2, etc.).
    - preset=None (or "permissive"/"none"): apply only explicitly provided cutoffs; this is the most
      permissive / backward-compatible mode and returns nearly the full table (subject to any
      user-supplied thresholds).

    Presets are oriented toward target-group "activated" / upregulated signals (positive
    logFC + positive unspliced_excess_residual; direction defaults to "up").
    For downregulated or two-sided selection from differential_expression() results,
    pass preset=None + logfc_direction="down" or "both" (with your desired logfc_cutoff).

    If you explicitly pass any cutoff parameter, it takes precedence over the preset.

    Calling with no arguments (or only the DataFrame) and no preset returns the full
    `all_results` (fully permissive).

    Only filters corresponding to columns present in the DataFrame are applied.
    This is safe whether or not `use_permutation=True` or `use_mixed_model=True` was used.

    New options for power users:
    - return_mask=True: return the boolean mask (pd.Series) instead of the filtered DataFrame.
      Useful to combine with other logic or apply yourself.
    - inplace=True: mutate the input `results` DataFrame in-place (keeps only passing rows
      and re-sorts). Returns the (mutated) DataFrame. Use with caution; the input reference
      will be modified. Ignored when return_mask=True.

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
    unspliced_excess_residual_cutoff : float
        Minimum bias-corrected unspliced (nascent) excess residual.
    logfc_cutoff : float
        Magnitude threshold for logFC (treated as non-negative). See logfc_direction.
    logfc_direction : {"up", "down", "both"}
        Direction filter applied when a "logFC" column is present:
        - "up" (default): keep if logFC > logfc_cutoff (upregulated in target)
        - "down": keep if logFC < -logfc_cutoff (downregulated in target)
        - "both": keep if |logFC| > logfc_cutoff (differentially expressed either way)
        Presets and this helper's design are "active"/up-biased by default.
        Use direction="down" or "both" with preset=None for standalone DE results.
    active_score_fdr_cutoff : float or None
        If the column exists, max permutation FDR on the composite active_score (ranking aid).
    unspliced_excess_fdr_cutoff : float or None
        If the column exists (use_permutation=True), max permutation FDR on
        ``unspliced_excess_residual`` (recommended for final gene lists).
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
        Subset of the input. Sorted by active_score (desc) when present;
        otherwise by p_adj then logFC (direction-aware for pure DE tables).
    """
    if not isinstance(results, pd.DataFrame):
        raise ValueError("results must be the all_results DataFrame returned by active_score")

    velocity_residual_cutoff = _NOT_PROVIDED
    if "velocity_residual_cutoff" in deprecated_kwargs:
        velocity_residual_cutoff = deprecated_kwargs.pop("velocity_residual_cutoff")
        warnings.warn(
            "velocity_residual_cutoff is deprecated; use unspliced_excess_residual_cutoff.",
            DeprecationWarning,
            stacklevel=2,
        )
    if deprecated_kwargs:
        raise TypeError(
            "filter_active_genes() got unexpected keyword argument(s): "
            f"{', '.join(sorted(deprecated_kwargs))}"
        )

    # Resolve values from preset + explicit overrides
    if preset is not None:
        p = preset.lower()
        if p in ("heuristic", "single_cell", "default"):
            preset_vals = {
                "active_score_cutoff": 55.0,
                "pval_cutoff": 0.05,
                "velocity_residual_cutoff": 1.0,
                "unspliced_excess_residual_cutoff": 1.0,
                "logfc_cutoff": 0.35,
                "active_score_fdr_cutoff": 0.25,
                "unspliced_excess_fdr_cutoff": 0.05,
                "effective_gamma_min": 0.05,
                "effective_gamma_max": 1.0,
                "delta_variance_min": None,
            }
        elif p in ("pseudobulk", "bulk"):
            preset_vals = {
                "active_score_cutoff": 5.0,
                "pval_cutoff": 0.05,
                "velocity_residual_cutoff": 0.05,
                "unspliced_excess_residual_cutoff": 0.05,
                "logfc_cutoff": 0.2,
                "active_score_fdr_cutoff": 0.25,
                "unspliced_excess_fdr_cutoff": 0.05,
                "effective_gamma_min": 0.05,
                "effective_gamma_max": 1.0,
                "delta_variance_min": None,
            }
        elif p in ("permissive", "none", "all", "no_filter"):
            preset_vals = {
                "active_score_cutoff": 0.0,
                "pval_cutoff": float("inf"),
                "velocity_residual_cutoff": float("-inf"),
                "unspliced_excess_residual_cutoff": float("-inf"),
                "logfc_cutoff": float("inf"),
                "active_score_fdr_cutoff": float("inf"),
                "unspliced_excess_fdr_cutoff": float("inf"),
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

    def _coerce_numeric_cutoff(val: Any, default: float, name: str) -> float:
        if val is _NOT_PROVIDED:
            return float(default)
        try:
            out = float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be numeric, got {val!r}") from exc
        return out

    active_score_cutoff = _coerce_numeric_cutoff(
        _resolve("active_score_cutoff", active_score_cutoff, 0.0),
        0.0,
        "active_score_cutoff",
    )
    pval_cutoff = _coerce_numeric_cutoff(
        _resolve("pval_cutoff", pval_cutoff, float("inf")), float("inf"), "pval_cutoff"
    )
    if pval_cutoff < 0 or (not math.isfinite(pval_cutoff) and not math.isinf(pval_cutoff)):
        raise ValueError("pval_cutoff must be non-negative, finite, or +inf (permissive).")
    if (
        velocity_residual_cutoff is not _NOT_PROVIDED
        and unspliced_excess_residual_cutoff is _NOT_PROVIDED
    ):
        unspliced_excess_residual_cutoff = velocity_residual_cutoff
    velocity_residual_cutoff = _resolve(
        "velocity_residual_cutoff", velocity_residual_cutoff, float("-inf")
    )
    unspliced_excess_residual_cutoff = _resolve(
        "unspliced_excess_residual_cutoff",
        unspliced_excess_residual_cutoff,
        velocity_residual_cutoff,
    )
    logfc_cutoff = _coerce_numeric_cutoff(
        _resolve("logfc_cutoff", logfc_cutoff, float("inf")), float("inf"), "logfc_cutoff"
    )
    # logfc_direction is not preset-driven (presets remain up-biased); normalize here
    dir_raw = str(logfc_direction).lower() if logfc_direction is not None else "up"
    if dir_raw in {"up", "positive", "pos", "u"}:
        direction = "up"
    elif dir_raw in {"down", "negative", "neg", "d"}:
        direction = "down"
    elif dir_raw in {"both", "two_sided", "twosided", "abs", "absolute", "either", "any", "b"}:
        direction = "both"
    else:
        raise ValueError(
            f'logfc_direction={logfc_direction!r} not recognized. Use one of: "up", "down", "both".'
        )
    active_score_fdr_cutoff = _resolve(
        "active_score_fdr_cutoff", active_score_fdr_cutoff, float("inf")
    )
    unspliced_excess_fdr_cutoff = _resolve(
        "unspliced_excess_fdr_cutoff", unspliced_excess_fdr_cutoff, float("inf")
    )
    effective_gamma_min = _resolve("effective_gamma_min", effective_gamma_min, float("-inf"))
    effective_gamma_max = _resolve("effective_gamma_max", effective_gamma_max, None)
    delta_variance_min = _resolve("delta_variance_min", delta_variance_min, None)

    df = results.copy()
    mask = pd.Series(True, index=df.index)

    # Core filters
    # Prefer adjusted p-value when present (consistent with active_score internal significant mask
    # and common user expectation). Fall back to nominal p_val only if p_adj is absent.
    if "active_score" in df.columns:
        active_vals = pd.to_numeric(df["active_score"], errors="coerce")
        mask &= active_vals.notna() & (active_vals >= active_score_cutoff)
    if "p_adj" in df.columns:
        padj_vals = pd.to_numeric(df["p_adj"], errors="coerce")
        mask &= padj_vals.notna() & (padj_vals < pval_cutoff)
    elif "p_val" in df.columns:
        pval_vals = pd.to_numeric(df["p_val"], errors="coerce")
        mask &= pval_vals.notna() & (pval_vals < pval_cutoff)
    residual_col = (
        UNSPLICED_EXCESS_RESIDUAL_COL
        if UNSPLICED_EXCESS_RESIDUAL_COL in df.columns
        else LEGACY_VELOCITY_RESIDUAL_COL
    )
    if residual_col in df.columns:
        resid_vals = pd.to_numeric(df[residual_col], errors="coerce")
        mask &= resid_vals.notna() & (resid_vals > unspliced_excess_residual_cutoff)
    if "logFC" in df.columns:
        lc = logfc_cutoff
        if isinstance(lc, (int, float)) and math.isinf(lc):
            # permissive preset: no logFC threshold
            pass
        else:
            if not math.isfinite(lc):
                raise ValueError("logfc_cutoff must be finite or +inf (permissive).")
            if lc < 0:
                lc = -lc  # always use positive magnitude
            logfc_vals = pd.to_numeric(df["logFC"], errors="coerce")
            if direction == "up":
                mask &= logfc_vals.notna() & (logfc_vals > lc)
            elif direction == "down":
                mask &= logfc_vals.notna() & (logfc_vals < -lc)
            else:  # both
                mask &= logfc_vals.notna() & (logfc_vals.abs() > lc)

    # Permutation FDR on composite score (optional ranking filter)
    if active_score_fdr_cutoff is not None and "active_score_fdr" in df.columns:
        fdr = pd.to_numeric(df["active_score_fdr"], errors="coerce")
        if math.isinf(active_score_fdr_cutoff):
            # Permissive: NaN FDR = not computed / no permutation — do not drop.
            mask &= fdr.isna() | (fdr < active_score_fdr_cutoff)
        else:
            mask &= fdr.notna() & (fdr < active_score_fdr_cutoff)

    # Permutation FDR on unspliced excess residual (recommended significance filter)
    if unspliced_excess_fdr_cutoff is not None and UNSPLICED_EXCESS_FDR_COL in df.columns:
        fdr = pd.to_numeric(df[UNSPLICED_EXCESS_FDR_COL], errors="coerce")
        if math.isinf(unspliced_excess_fdr_cutoff):
            mask &= fdr.isna() | (fdr < unspliced_excess_fdr_cutoff)
        else:
            mask &= fdr.notna() & (fdr < unspliced_excess_fdr_cutoff)

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
    # otherwise fall back to p_adj, then logFC according to direction
    # (strongest up first, strongest down first, or largest |logFC| for both).
    if "active_score" in filtered.columns:
        filtered = filtered.sort_values("active_score", ascending=False)
    elif "p_adj" in filtered.columns:
        sort_cols = ["p_adj"]
        ascending = [True]
        if "logFC" in filtered.columns:
            if direction == "down":
                sort_cols.append("logFC")
                ascending.append(True)  # most negative first
            elif direction == "both":
                filtered = filtered.assign(_logfc_abs=filtered["logFC"].abs())
                sort_cols.append("_logfc_abs")
                ascending.append(False)
            else:
                sort_cols.append("logFC")
                ascending.append(False)
        filtered = filtered.sort_values(sort_cols, ascending=ascending)
        if "_logfc_abs" in filtered.columns:
            filtered = filtered.drop(columns=["_logfc_abs"])
    else:
        filtered = filtered.sort_values(filtered.columns[0])

    if return_mask:
        return mask

    if inplace:
        # Mutate caller's DataFrame: drop non-matching rows, keep sorted order
        to_drop = results.index.difference(filtered.index)
        if len(to_drop) > 0:
            results.drop(index=to_drop, inplace=True)
        # Re-apply the sort in place on the surviving rows
        if "active_score" in results.columns:
            results.sort_values("active_score", ascending=False, inplace=True)
        elif "p_adj" in results.columns:
            sort_cols = ["p_adj"]
            ascending = [True]
            if "logFC" in results.columns:
                if direction == "down":
                    sort_cols.append("logFC")
                    ascending.append(True)
                elif direction == "both":
                    # Use direct assignment so inplace=True actually mutates the caller's
                    # DataFrame (assign() would rebind the local variable and lose the effect).
                    results["_logfc_abs"] = results["logFC"].abs()
                    sort_cols.append("_logfc_abs")
                    ascending.append(False)
                else:
                    sort_cols.append("logFC")
                    ascending.append(False)
            results.sort_values(sort_cols, ascending=ascending, inplace=True)
            if "_logfc_abs" in results.columns:
                results.drop(columns=["_logfc_abs"], inplace=True)
        return results

    return filtered


WORKFLOW_PRESETS: dict[str, dict[str, Any]] = {
    "explore": {
        "label": "Quick exploration (ranking only, no permutation)",
        "active_score_kwargs": {
            "use_permutation": False,
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "heuristic",
    },
    "report": {
        "label": "Manuscript reporting with permutation FDR on unspliced excess",
        "active_score_kwargs": {
            "use_permutation": True,
            "n_perm": 500,
            "perm_de_backend": "same",
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "heuristic",
    },
    "pseudobulk_report": {
        "label": "Multi-replicate pseudobulk DE + permutation FDR",
        "active_score_kwargs": {
            "use_pseudobulk": True,
            "pseudobulk_de_backend": "pydeseq2",
            "use_permutation": True,
            "n_perm": 200,
            "perm_de_backend": "same",
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "pseudobulk",
    },
    "nascent_focus": {
        "label": "Rank by bias-corrected nascent excess residual only",
        "active_score_kwargs": {
            "ranking_mode": "nascent_excess",
            "use_permutation": False,
            "mode": "heuristic",
        },
        "filter_preset": "heuristic",
    },
}


def _permutation_power_guidance(
    n_cells_target: int,
    n_cells_reference: int,
    n_genes: int,
    *,
    n_samples_target: int | None = None,
    n_samples_reference: int | None = None,
    n_perm: int = 100,
    n_jobs: int = -1,
) -> dict[str, Any]:
    """Rough runtime / permutation-space guidance for diagnose_design."""
    import os

    is_pseudobulk = (
        n_samples_target is not None
        and n_samples_reference is not None
        and min(n_samples_target, n_samples_reference) >= 2
    )
    max_exact: int | None = None
    if is_pseudobulk:
        n_t = max(1, int(n_samples_target or 1))
        n_r = max(1, int(n_samples_reference or 1))
        if n_t + n_r <= 30:
            max_exact = max(1, comb(n_t + n_r, n_t) - 1)

    cpu = os.cpu_count() or 4
    effective_cores = cpu if n_jobs == -1 else max(1, min(n_jobs, cpu))
    sec_per_perm = max(0.03, (n_genes / 5000.0) * 0.10)
    est_minutes = n_perm * sec_per_perm / effective_cores / 60.0

    notes: list[str] = [
        "Estimates assume heuristic mode on a typical laptop/workstation; "
        "advanced mode, Memento, or perm_de_backend='same' with PyDESeq2 can be slower.",
    ]
    if is_pseudobulk and max_exact is not None and max_exact < n_perm:
        notes.append(
            f"Pseudobulk design allows at most {max_exact} exact label permutations; "
            f"auto_adjust_n_perm=True will cap n_perm accordingly."
        )
    if min(n_cells_target, n_cells_reference) < 50:
        notes.append(
            "Very small cell counts: permutation FDR on unspliced excess will have low power "
            "even if runtime is acceptable."
        )

    return {
        "n_genes": int(n_genes),
        "default_n_perm": int(n_perm),
        "is_pseudobulk_context": is_pseudobulk,
        "max_exact_permutations_pseudobulk": max_exact,
        "estimated_runtime_minutes_heuristic": round(est_minutes, 1),
        "effective_cores_assumed": effective_cores,
        "notes": notes,
    }


def diagnose_design(
    adata_input: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None = None,
    _min_cells_per_sample: int = 10,  # reserved/internal; not yet used (will affect sample filtering/power in future)
    *,
    copy_input: bool = True,
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
      - suggested_preset: filter_active_genes preset ("heuristic" or "pseudobulk")
      - power_summary: permutation runtime / power guidance dict
      - workflow_preset: recommended entry from WORKFLOW_PRESETS

    Note: ``_min_cells_per_sample`` is reserved for future use and currently ignored.

    copy_input : bool, default True
        If True (default), a full deep copy of the input AnnData is made before
        reading. Set False for a zero-copy read-only diagnostic when the caller
        guarantees the input will not be mutated (saves large amounts of memory
        and time on big datasets with many layers).
    """
    adata = adata_input.copy() if copy_input else adata_input

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs")

    norm_groups = adata.obs[groupby].map(_normalize_group_label)
    target_mask = norm_groups == _normalize_group_label(target_group)
    ref_mask = norm_groups == _normalize_group_label(reference_group)

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
        "workflow_preset": "explore",
        "power_summary": None,
    }

    # Global unspliced fraction (important technical QC)
    layer_keys = _resolve_velocity_layer_keys(adata)
    if layer_keys is not None:
        spliced_key, unspliced_key = layer_keys
        try:
            ufrac = _qc.unspliced_global(
                adata,
                spliced_key=spliced_key,
                unspliced_key=unspliced_key,
                warn_threshold=0.5,
            )
            result["unspliced_global_fraction"] = float(ufrac)
            if ufrac > 0.5:
                result["warnings"].append(
                    f"Global unspliced fraction is high ({ufrac:.1%}). "
                    "This often indicates nuclear enrichment or gDNA contamination and can "
                    "reduce the reliability of velocity-based signals."
                )
        except Exception as e:
            logger.debug("Could not compute global unspliced fraction: %s", e)

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

    if n_r < 80:
        result["recommendations"].append(
            f"Reference group has {n_r} cells — consider gamma_method='empirical_bayes' "
            "for more stable per-gene reference gamma (robust log-ratio shrinkage)."
        )

    # Suggest filter_active_genes preset
    if result["suggested_preset"] is None:
        if sample_col and result.get("n_samples_target", 0) >= 5:
            result["suggested_preset"] = "pseudobulk"
        else:
            result["suggested_preset"] = "heuristic"

    # Workflow preset (ties to WORKFLOW_PRESETS for active_score kwargs)
    if (
        sample_col
        and min(result.get("n_samples_target") or 0, result.get("n_samples_reference") or 0) >= 3
    ):
        result["workflow_preset"] = "pseudobulk_report"
    elif min(n_t, n_r) >= 100 and not result["warnings"]:
        result["workflow_preset"] = "report"
    else:
        result["workflow_preset"] = "explore"

    result["power_summary"] = _permutation_power_guidance(
        n_t,
        n_r,
        adata.n_vars,
        n_samples_target=result.get("n_samples_target"),
        n_samples_reference=result.get("n_samples_reference"),
    )

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

    ps = result.get("power_summary") or {}
    if ps.get("estimated_runtime_minutes_heuristic") is not None:
        logger.info(
            "  Permutation runtime (heuristic, n_perm=%d): ~%.1f min on %d cores",
            ps.get("default_n_perm", 100),
            ps["estimated_runtime_minutes_heuristic"],
            ps.get("effective_cores_assumed", 1),
        )

    return result


def recommend_workflow(
    adata_input: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None = None,
) -> dict[str, Any]:
    """
    High-level recommendation for analysis path based on experimental design.

    This is a thin, user-friendly wrapper around diagnose_design that returns
    actionable preset + backend suggestions.

    Returns keys:
      - workflow_preset: key into ``WORKFLOW_PRESETS`` (e.g. ``"pseudobulk_report"``)
      - preset_config: full preset dict (label, active_score_kwargs, filter_preset)
      - recommended_preset: ``filter_active_genes`` preset name
      - de_backend: ``"scanpy"`` | ``"pydeseq2"`` | ...
      - suggested_kwargs: merged kwargs for :func:`active_score`
      - warnings, recommendations, power_summary
      - full_diagnosis: raw dict from :func:`diagnose_design`

    Example:
        rec = scat.recommend_workflow(adata, "condition", "GA", "Ctrl", sample_col="sample")
        adata, sig, res = scat.active_score(
            adata, groupby="condition", target_group="GA", reference_group="Ctrl",
            **rec["suggested_kwargs"],
        )
        candidates = scat.filter_active_genes(res, preset=rec["filter_preset"])
    """
    diag = diagnose_design(
        adata_input,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        copy_input=True,  # public entry point: be safe by default
    )

    workflow_key = diag.get("workflow_preset") or "explore"
    preset_config = WORKFLOW_PRESETS.get(workflow_key, WORKFLOW_PRESETS["explore"])
    filter_preset = preset_config.get("filter_preset", diag.get("suggested_preset") or "heuristic")

    suggested_kwargs: dict[str, Any] = dict(preset_config.get("active_score_kwargs", {}))
    if sample_col and sample_col in getattr(adata_input, "obs", pd.DataFrame()).columns:
        suggested_kwargs.setdefault("sample_col", sample_col)
    if suggested_kwargs.get("use_pseudobulk"):
        de_backend = "pydeseq2"
        suggested_kwargs.setdefault("pseudobulk_de_backend", "pydeseq2")
    else:
        de_backend = "scanpy"
        suggested_kwargs.setdefault("de_method", "wilcoxon")

    rec = {
        "workflow_preset": workflow_key,
        "preset_config": preset_config,
        "recommended_preset": filter_preset,
        "filter_preset": filter_preset,
        "de_backend": de_backend,
        "use_permutation": bool(suggested_kwargs.get("use_permutation", False)),
        "suggested_kwargs": suggested_kwargs,
        "warnings": diag.get("warnings", []),
        "recommendations": diag.get("recommendations", []),
        "power_summary": diag.get("power_summary"),
        "full_diagnosis": diag,
    }

    logger.info(
        "Workflow recommendation: workflow_preset=%s, filter_preset=%s, de_backend=%s",
        rec["workflow_preset"],
        rec["filter_preset"],
        rec["de_backend"],
    )
    return rec


def _maybe_add_gene_features(adata: Any, organism: str) -> Any:
    """Attach bundled gene features when length/intron columns are missing or completely empty.

    Semantics:
    - If either column is absent → attach (fill NaNs for missing genes).
    - If both columns present but *all* values are NaN (gl.notna().any() is False
      for both) → attach (user provided empty placeholders).
    - If at least one gene has a real (non-NaN) value in *both* columns → do nothing.
      Partial user data is respected; we never overwrite or "complete" the table.
    This boundary is deliberate but worth knowing: a table with 99% NaN but 1%
    real value will *not* trigger auto-attachment.
    """
    has_length = "gene_length" in adata.var.columns
    has_intron = "intron_number" in adata.var.columns
    needs = not has_length or not has_intron
    if has_length and has_intron:
        gl = pd.to_numeric(adata.var["gene_length"], errors="coerce")
        intr = pd.to_numeric(adata.var["intron_number"], errors="coerce")
        needs = not (gl.notna().any() and intr.notna().any())
    if needs:
        from .pp_bias import add_gene_features

        add_gene_features(adata, organism=organism)
        logger.info("Attached bundled gene features (organism=%s) for bias correction.", organism)
    return adata


def _resolve_simple_backend_kwargs(
    adata: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None,
) -> dict[str, Any]:
    """Pick pseudobulk vs single-cell defaults from replicate structure."""
    kwargs: dict[str, Any] = {
        "de_method": "wilcoxon",
        "use_pseudobulk": False,
        "sample_col": None,
        "pseudobulk_de_backend": "pydeseq2",
    }
    if sample_col and sample_col in adata.obs.columns:
        norm_groups = adata.obs[groupby].map(_normalize_group_label)
        t_mask = norm_groups == _normalize_group_label(target_group)
        r_mask = norm_groups == _normalize_group_label(reference_group)
        n_s_t = int(adata.obs.loc[t_mask, sample_col].nunique())
        n_s_r = int(adata.obs.loc[r_mask, sample_col].nunique())
        if min(n_s_t, n_s_r) >= 3:
            kwargs.update(
                {
                    "use_pseudobulk": True,
                    "sample_col": sample_col,
                    "de_method": "t-test_overestim_var",
                }
            )
            logger.info(
                "Simple path: detected >=3 samples per group — using pseudobulk + PyDESeq2."
            )
        else:
            logger.info(
                "Simple path: few samples per group (target=%d, reference=%d) — "
                "using single-cell Wilcoxon DE.",
                n_s_t,
                n_s_r,
            )
    return kwargs


def active_score_simple(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",  # convenience defaults (distinct from core active_score)
    sample_col: str | None = None,
    organism: str = "mouse",
    *,
    show_plot: bool = False,
    copy_input: bool = True,
    pydeseq2_min_counts: int = 10,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Recommended entry point for new users (minimal parameters).

    Wraps :func:`active_score` with sensible defaults:
    - Uses "Disease"/"Control" as group defaults (unlike core active_score which defaults
      to the historical "GA"/"Ctrl").
    - heuristic mode, no permutation (inspect ``all_results`` + ``filter_active_genes``)
    - auto-attaches bundled gene features when missing
    - pseudobulk + PyDESeq2 when ``sample_col`` has >=3 replicates per group;
      otherwise single-cell Wilcoxon DE

    For full control (permutation, advanced mode, mixed models, etc.) use
    :func:`active_score` directly.
    """
    adata = adata_input.copy() if copy_input else adata_input
    _maybe_add_gene_features(adata, organism)
    backend = _resolve_simple_backend_kwargs(
        adata, groupby, target_group, reference_group, sample_col
    )
    return active_score(
        adata,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        mode="heuristic",
        use_permutation=False,
        show_plot=show_plot,
        copy_input=False,
        pydeseq2_min_counts=pydeseq2_min_counts,
        **backend,
    )


def differential_expression_simple(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",
    sample_col: str | None = None,
    *,
    copy_input: bool = True,
    pydeseq2_min_counts: int = 10,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Minimal-parameter differential expression (no velocity layers required).

    Same backend auto-selection as :func:`active_score_simple`.
    For Memento, mixed models, or custom preprocess use :func:`differential_expression`.
    """
    adata = adata_input.copy() if copy_input else adata_input
    backend = _resolve_simple_backend_kwargs(
        adata, groupby, target_group, reference_group, sample_col
    )
    return differential_expression(
        adata,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        copy_input=False,
        pydeseq2_min_counts=pydeseq2_min_counts,
        **backend,
    )


def run_default_pipeline(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",  # convenience defaults; always prefer explicit values
    sample_col: str | None = None,
    organism: str = "mouse",
    *,
    run_go_enrichment: bool = True,
    gene_sets: str = "GO_Biological_Process",
    filter_preset: str | None = None,
    show_plot: bool = False,
) -> dict[str, Any]:
    """
    End-to-end recommended workflow for first-time users.

    Steps: active scoring → ``filter_active_genes`` → optional GO enrichment.
    Uses "Disease"/"Control" convenience defaults for target/reference.

    The default for ``filter_preset`` is now **auto-detected** from the experimental
    design (via ``_resolve_simple_backend_kwargs``): "pseudobulk" when sample_col
    is provided with >=3 samples per group (which triggers pseudobulk inside
    active_score_simple), otherwise "heuristic". This keeps the thresholds
    consistent with the actual scale of active_score / unspliced_excess_residual
    (see WORKFLOW_PRESETS["pseudobulk_report"]).

    Returns a dict with keys:
      - ``adata``, ``significant``, ``all_results``, ``candidates``
      - ``enrichment`` (DataFrame or None)
      - ``filter_preset``, ``backend`` (kwargs used for DE)
    """
    # Resolve once so we can pick a matching filter_preset (addresses mismatch
    # between auto-pseudobulk in active_score_simple and hardcoded "heuristic").
    backend = _resolve_simple_backend_kwargs(
        adata_input, groupby, target_group, reference_group, sample_col
    )
    if filter_preset is None:
        filter_preset = "pseudobulk" if backend.get("use_pseudobulk") else "heuristic"

    adata_res, significant, all_results = active_score_simple(
        adata_input,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        organism=organism,
        show_plot=show_plot,
    )
    candidates = filter_active_genes(all_results, preset=filter_preset)

    enrichment = None
    if run_go_enrichment and len(candidates) > 0:
        from .enrich import run_enrichment

        enrichment = run_enrichment(
            candidates.index.tolist(),
            gene_sets=gene_sets,
            organism=organism,
            adata=adata_res,
            pval_cutoff=0.05,
        )

    # We already resolved backend above (avoids calling the resolver a second time
    # just to "guess" what active_score_simple decided internally).
    return {
        "adata": adata_res,
        "significant": significant,
        "all_results": all_results,
        "candidates": candidates,
        "enrichment": enrichment,
        "filter_preset": filter_preset,
        "backend": backend,
    }
