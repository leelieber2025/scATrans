"""
scATrans tl module.

Primary analysis function `active_score` for computing active transcription scores
from velocity and differential expression data. Implementation details are
distributed across private supporting modules.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

import anndata as ad
import joblib
import numpy as np
import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests

from ._bias import fit_huber_bias_correction
from ._de import _run_de_wrapper
from ._permutation import _single_permutation_task
from ._utils import (
    _get_exponential_scale_lambda,
    _normalize_velocity_layers_by_size_factor,
    _pseudobulk_with_layers,
    _soft_scale,
    comb,  # for small-n permutation space calculation
)
from ._velocity import _compute_moments_velocity_delta, _compute_velocity_delta

# qc is imported lazily inside active_score to keep startup light, but exposed at package level
from . import qc as _qc  # for unspliced_global integration

try:
    from . import _version

    VERSION = _version.version
except (ImportError, AttributeError):
    VERSION = "0.7.0.dev0"

logger = logging.getLogger(__name__)


def active_score(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "GA",
    reference_group: str = "Ctrl",
    subset_col: Optional[str] = None,
    subset_values: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
    weight_fc: float = 1.0,
    weight_unspliced: float = 1.0,
    weight_pval: float = 1.0,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.5,
    active_fdr_cutoff: float = 0.05,
    de_method: str = "t-test_overestim_var",
    pseudobulk_de_backend: str = "pydeseq2",
    use_permutation: bool = False,
    perm_de_backend: str = "fast",
    n_perm: int = 100,
    n_jobs: int = -1,
    de_preprocess: str = "auto",
    gene_type_filter: Optional[str] = None,
    use_pseudobulk: bool = False,
    sample_col: Optional[str] = None,
    min_cells: int = 10,
    min_counts: int = 1000,
    pb_x_layer: str = "spliced",
    pb_use_total_for_x: bool = True,
    min_total_counts: int = 50,
    random_seed: int = 42,
    show_plot: bool = True,
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
) -> Tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Identify genes that are actively transcribed using a composite "active score".

    The method integrates:
    - Differential expression (logFC and p-value) between target and reference groups.
    - Velocity signal approximated by the difference in unspliced abundance after
      a simple reference-based gamma correction (with optional scVelo moments smoothing
      in "advanced" mode). The per-gene effective gamma is stored in .var["effective_gamma"]
      for transparency.
    - Bias correction for gene length and intron number via Huber regression on the
      velocity delta (or median fallback). Detailed fit diagnostics (coefficients, n_genes_used,
      fallback status) are recorded in adata.uns["scatrans"]["diagnostics"]["bias_correction"].
    - Optional permutation testing for gene-level significance and FDR (note: for speed,
      velocity layers/gammas are computed once on the original labeling and only group
      assignments are permuted).
    - Optional mixed linear model (LMM) DE via statsmodels.mixedlm when use_mixed_model=True
      (requires sample_col). This models cell-level data with sample as random intercept
      (~ condition + (1 | sample)), providing replicate-aware p-values (avoids pseudoreplication).
      Also computes delta_variance (fraction of variance explained by condition, variancePartition-style)
      and a LRT delta_var_pval. These are stored in .var / all_results and can be used as
      supplementary filter. See README for relation to Libra, dreamlet/dreampy, NEBULA.
    - Rich run-time diagnostics (global unspliced fraction, bias fit quality, etc.) are
      always stored under adata.uns["scatrans"]["diagnostics"] and a concise summary is
      logged at completion.

    A soft-scaled weighted combination of the three signals produces the final
    active score (0–100). Results (including the processed AnnData) are returned
    together with tables of significant and all genes.

    Full parameter reference, plotting functions, enrichment, and usage examples
    are documented in the package README.
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

    if not isinstance(use_mixed_model, bool):
        raise ValueError("use_mixed_model must be boolean.")
    if not isinstance(use_delta_variance_pval, bool):
        raise ValueError("use_delta_variance_pval must be boolean.")
    if not (0 < delta_var_pval_cutoff < 1):
        raise ValueError("delta_var_pval_cutoff must be in (0, 1).")

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

    if spliced_layer != "spliced" or unspliced_layer != "unspliced":
        if spliced_layer in adata_input.layers and unspliced_layer in adata_input.layers:
            adata_input.layers["spliced"] = adata_input.layers[spliced_layer].copy()
            adata_input.layers["unspliced"] = adata_input.layers[unspliced_layer].copy()
            logger.info(
                "Layer remapping applied: '%s' → 'spliced', '%s' → 'unspliced' (internal use only)",
                spliced_layer,
                unspliced_layer,
            )

    if "spliced" not in adata_input.layers or "unspliced" not in adata_input.layers:
        raise ValueError("Both 'spliced' and 'unspliced' layers are required after layer handling.")

    keep_mask = adata_input.obs[groupby].astype(str).isin([target_group, reference_group])
    adata = adata_input[keep_mask].copy()

    if gene_type_filter:
        if "gene_type" not in adata.var.columns:
            raise ValueError("'gene_type_filter' provided but 'gene_type' column is missing.")
        adata = adata[:, adata.var["gene_type"] == gene_type_filter].copy()

    if adata.n_vars == 0:
        raise ValueError("No genes remain after filtering.")

    # Mixed model requirements (cell-level RE; pseudobulk + count DE is separate path)
    if use_mixed_model:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_mixed_model=True (for the random effect grouping).")
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
    )

    adata.var["logFC"] = de_df["logFC"]
    adata.var["p_val"] = de_df["p_val"]
    adata.var["p_adj"] = de_df["p_adj"]
    if "delta_variance" in de_df.columns:
        adata.var["delta_variance"] = de_df["delta_variance"]
    if "delta_var_pval" in de_df.columns:
        adata.var["delta_var_pval"] = de_df["delta_var_pval"]

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
    moments_info: Dict[str, Any] = {}
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
            delta_velocity, total_us_velocity, gamma_ref, moments_info = _compute_moments_velocity_delta(
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
    )

    adata.var["velocity_delta_raw"] = delta_velocity
    adata.var["velocity_residual"] = residual
    adata.var["total_us_counts"] = total_us_raw
    adata.var["total_us_counts_raw"] = total_us_raw
    adata.var["total_us_counts_velocity_layer"] = total_us_velocity
    adata.var["valid_expr"] = valid_expr
    adata.var["velocity_source"] = velocity_source
    adata.var["effective_gamma"] = gamma_ref  # per-gene reference gamma used for delta (transparency)

    # ==================== DIAGNOSTICS (high priority for usability & paper rigor) ====================
    n_valid_bias = int(bias_info.get("n_genes_used_for_fit", 0))
    diagnostics: Dict[str, Any] = {
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "n_genes_with_valid_features": int(valid_feat.sum()),
        "unspliced_global_fraction": float(unspliced_fraction) if unspliced_fraction is not None else np.nan,
        "bias_correction": bias_info,
        "velocity": {
            "source": velocity_source,
            "n_genes_with_finite_delta": int(np.isfinite(delta_velocity).sum()),
        },
        "mixed_model": {
            "used": bool(use_mixed_model),
            "sample_col": sample_col if use_mixed_model else None,
            "n_samples": int(adata.obs[sample_col].nunique()) if (use_mixed_model and sample_col and sample_col in adata.obs.columns) else None,
            "delta_variance_available": "delta_variance" in adata.var.columns,
            "median_delta_variance": float(np.nanmedian(adata.var["delta_variance"])) if "delta_variance" in adata.var.columns else np.nan,
        },
    }
    if mode == "advanced" and moments_info:
        diagnostics["velocity"]["moments"] = {
            k: moments_info.get(k)
            for k in ("n_neighbors_effective", "n_pcs_effective", "used_precomputed_moments", "neighbors_source")
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
            # active_fdr_disabled_reason kept for metadata / future use if needed
            _ = "small_permutation_space"

    # ==================== METADATA ====================
    velocity_delta_layer = (
        "scvelo_Mu_Ms_moments"
        if velocity_source.startswith("scvelo_moments")
        else (
            "size_factor_normalized_spliced_unspliced" if is_pseudobulk else "raw_spliced_unspliced"
        )
    )

    adata.uns["scatrans"] = {
        "version": VERSION,
        "groupby": groupby,
        "target_group": target_group,
        "reference_group": reference_group,
        "mode": mode,
        "velocity_source": velocity_source,
        "velocity_delta_layer": velocity_delta_layer,
        "advanced_fallback": advanced_fallback,
        "advanced_use_precomputed": advanced_use_precomputed,
        "allow_advanced_pseudobulk": allow_advanced_pseudobulk,
        "advanced_recompute_neighbors": advanced_recompute_neighbors,
        "advanced_experimental": mode == "advanced",
        "moments_info": moments_info if mode == "advanced" else None,
        "advanced_neighbor_graph_basis": "adata.X_after_de_preprocess"
        if mode == "advanced"
        else None,
        "advanced_layer_preprocessing": (
            "existing_spliced_unspliced_layers_no_scv_filter_and_normalize"
            if mode == "advanced"
            else None
        ),
        "use_permutation": use_permutation,
        "n_perm": int(n_perm) if use_permutation else 0,
        "weight_fc": weight_fc,
        "weight_unspliced": weight_unspliced,
        "weight_pval": weight_pval,
        "pval_cutoff": pval_cutoff,
        "logfc_cutoff": logfc_cutoff,
        "active_fdr_cutoff": active_fdr_cutoff,
        "de_method": de_method,
        "pseudobulk_de_backend": pseudobulk_de_backend,
        "perm_de_backend": perm_de_backend if use_permutation else None,
        "prior_weight": prior_weight,
        "min_total_counts": min_total_counts,
        "random_seed": random_seed,
        "use_mixed_model": bool(use_mixed_model),
        "sample_col": sample_col if use_mixed_model else None,
        "use_delta_variance_pval": bool(use_delta_variance_pval),
        "delta_var_pval_cutoff": float(delta_var_pval_cutoff),
        # New rich diagnostics for usability and reproducibility (high priority)
        "diagnostics": diagnostics,
        # Explicit note on permutation approximation (important for paper & trust)
        "permutation_approximation_note": (
            "For efficiency, velocity layers (raw or Mu/Ms) and effective_gamma are fixed from the original (non-permuted) data. "
            "Only group labels are shuffled when recomputing delta and the composite active_score inside permutations."
        ) if use_permutation else None,
        "unspliced_global_fraction": float(unspliced_fraction) if unspliced_fraction is not None else np.nan,
    }

    # ==================== SIGNIFICANT GENES + RETURN ====================
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
        "effective_gamma",  # new transparency column (may be all-same in heuristic)
    ]
    if use_permutation:
        cols.extend(["active_score_pval", "active_score_fdr"])
    if "delta_variance" in adata.var.columns:
        cols.append("delta_variance")
    if "delta_var_pval" in adata.var.columns:
        cols.append("delta_var_pval")
    cols = [c for c in cols if c in adata.var.columns]

    mask = (
        (adata.var["p_adj"] < pval_cutoff)
        & (adata.var["logFC"] > logfc_cutoff)
        & (adata.var["velocity_residual"] > 0)
        & (adata.var["valid_expr"])
        & (adata.var["active_score"] > 0)
    )

    if use_permutation and use_fdr_for_significance:
        mask = mask & (adata.var["active_score_fdr"] < active_fdr_cutoff)

    # Delta Variance pval as optional supplementary filter (user-controlled)
    if use_delta_variance_pval and "delta_var_pval" in adata.var.columns:
        mask = mask & (adata.var["delta_var_pval"] < delta_var_pval_cutoff)

    significant = adata.var[mask][cols].copy().sort_values("active_score", ascending=False)
    all_results = adata.var[cols].copy().sort_values("active_score", ascending=False)

    logger.info(
        "Analysis completed in %s mode! Significant active genes: %d", mode, len(significant)
    )

    # ==================== USER-FACING RUN SUMMARY (diagnostics for convenience) ====================
    try:
        ufrac = diagnostics.get("unspliced_global_fraction", np.nan)
        bias = diagnostics.get("bias_correction", {})
        n_fit = bias.get("n_genes_used_for_fit", 0)
        fb = " (median fallback)" if bias.get("fallback_to_median") else ""
        disp = locals().get("display_mode", mode)
        logger.info(
            "Run summary — cells: %d | unspliced frac: %.1f%% | bias fit genes: %d%s | mode: %s | sig: %d",
            diagnostics.get("n_cells", 0),
            (ufrac * 100.0) if np.isfinite(ufrac) else float("nan"),
            n_fit,
            fb,
            disp,
            len(significant),
        )
        if use_permutation:
            logger.info("Permutation used %d iterations (velocity layers fixed from original labeling).", n_perm)
    except Exception:
        pass  # never break user run on summary logging

    # ==================== PLOTTING (now delegates to pl – no more inline duplication) ====================
    if show_plot:
        display_mode = mode
        if velocity_source == "heuristic_fallback_from_advanced":
            display_mode = "advanced→heuristic fallback"

        try:
            from . import pl

            # Prefer the rich publication-ready comet plot.
            # This replaces the previous ~35-line ad-hoc scatter that lived inside active_score.
            pl.comet_plot(
                all_results,
                top_n=12,
                title=f"scATrans Active Drivers ({display_mode})",
            )
        except Exception:
            # Never let plotting break the analysis
            logger.debug("show_plot=True but plotting failed (missing optional deps or display).")

    return adata, significant, all_results
