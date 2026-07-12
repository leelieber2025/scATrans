"""scatrans.tl.active — internal package module."""

from __future__ import annotations

import logging
import warnings
from typing import Any

import anndata as ad
import joblib
import numpy as np
import pandas as pd

from .. import qc as _qc
from .._de import _run_de_wrapper
from .._permutation import run_permutation_test
from .._utils import (
    UNSPLICED_EXCESS_DELTA_COL,
    UNSPLICED_EXCESS_FDR_COL,
    UNSPLICED_EXCESS_PVAL_COL,
    UNSPLICED_EXCESS_RESIDUAL_COL,
    _apply_de_preprocess,
    _as_contrast_categorical,
    _as_var_dataframe,
    _composite_active_score_terms,
    _get_exponential_scale_lambda,
    _lambda_pval_for_active_score,
    _matrix_copy,
    _matrix_sum_axis0,
    _normalize_group_label,
    _normalize_velocity_layers_by_size_factor,
    _pseudobulk_with_layers,
    _resolve_aligned_raw_counts,
    _score_direction_effect,
    _subset_obs_mask,
    _validate_group_contrast,
    _warn_if_negative_layer_values,
    _write_unspliced_excess_columns,
    comb,
)
from .._utils import _fit_huber_bias_correction as fit_huber_bias_correction
from .._velocity import _compute_moments_velocity_delta, _compute_velocity_delta
from ._common import (
    HEURISTIC_FILTER_DEFAULTS,
    PSEUDOBULK_FILTER_DEFAULTS,
    VERSION,
    _coerce_memento_de_preprocess,
    _materialize_if_view,
    _require_explicit_groups,
    _resolve_deprecated_active_score_kwargs,
    _select_obs,
    _select_var,
    _validate_de_common_options,
)
from .design import _emit_design_diagnosis_logs, diagnose_design
from .filter import _builtin_significant_mask

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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
    logfc_cutoff: float = float(HEURISTIC_FILTER_DEFAULTS["logfc_cutoff"] or 0.35),
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
    paired_replicates: bool = False,  # mixed model: same sample_col ID = same individual across conditions
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
    # "empirical_bayes": robust log-ratio empirical Bayes shrinkage (recommended for small reference).
    # Note: prior_weight affects empirical_bayes via count_pseudocount (observation precision), so its
    # impact is gentler than on heuristic_shrink/robust_median (direct additive shrinkage).
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
    is the main output; the built-in significant list uses the same default thresholds
    as :func:`filter_active_genes` with ``preset="heuristic"`` (logFC > 0.35,
    unspliced_excess_residual > 1.0, active_score >= 55, etc.) and may still be
    small on low-signal data.

    A separate function diagnose_design is available to summarize the experimental
    design and surface relevant warnings before analysis.

    **Important statistical note (reporting boundaries)**:
    - `active_score` is a **heuristic ranking score only**. It is NOT a p-value,
      effect size with calibrated uncertainty, or evidence of causal transcriptional activation.
      Composite legs for logFC and -log10(p_adj) are upregulation-gated
      (``logFC > 0``, or ``mixedlm_coef > 0`` when MixedLM is used so the gate matches
      what ``p_adj`` tests); downregulated genes do not receive score from the DE
      significance term alone. The p-value soft-scale λ is estimated on
      direction-positive genes only.
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
        If True (default), deep-copy after combining obs filters so the caller's
        object is not mutated. If False, avoid an *extra* intermediate full copy
        during subsetting when possible, but **always isolate** before any write
        (groupby labels, ``de_preprocess`` on ``.X``, layer remap, ``.var`` columns)
        so the caller's AnnData is never modified in place.
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
        paired_replicates=paired_replicates,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        min_cells=min_cells,
        min_counts=min_counts,
    )
    if use_memento_de and use_mixed_model:
        raise ValueError(
            "use_mixed_model=True and use_memento_de=True are incompatible. "
            "Choose one cell-level DE backend."
        )
    if use_memento_de and use_pseudobulk:
        raise ValueError(
            "use_memento_de=True is not supported with use_pseudobulk=True "
            "(Memento is a cell-level method-of-moments estimator)."
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
        # Always residual-only: custom weight_* would otherwise silently keep a composite
        # score while the mode name promises pure unspliced-excess ranking.
        if not default_weights and (
            weight_fc != 0.0 or weight_pval != 0.0 or weight_unspliced != 1.0
        ):
            logger.warning(
                "ranking_mode='nascent_excess' forces weight_fc=0, weight_pval=0, "
                "weight_unspliced=1 (overriding weight_fc=%s, weight_unspliced=%s, weight_pval=%s).",
                weight_fc,
                weight_unspliced,
                weight_pval,
            )
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
        subset_mask = _subset_obs_mask(adata_input.obs[subset_col], subset_values)
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

    # Automatic design guidance for small-sample or replicate-structured data.
    # Capture + re-emit warnings (do not discard the return value) and store under
    # diagnostics so low-replicate / low-power designs surface to the user.
    design_diag: dict[str, Any] | None = None
    if sample_col or use_pseudobulk:
        try:
            design_diag = diagnose_design(
                adata_input,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                sample_col=sample_col,
                copy_input=False,  # pure read-only diagnostic; avoid expensive deep copy
                emit_logs=False,  # active_score owns the log lines below
            )
            _emit_design_diagnosis_logs(design_diag, prefix="active_score")
        except Exception as e:
            logger.warning(
                "diagnose_design failed (non-fatal); continuing without design guidance: %s",
                e,
            )

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
    # Hard isolation: never mutate the caller's object (labels / .X / layers / .var).
    if adata is adata_input:
        adata = adata.copy()
        if not copy_input:
            logger.info(
                "copy_input=False: isolated a working copy before mutation so the "
                "caller's AnnData is left unchanged."
            )
    adata.obs[groupby] = norm_groups.loc[obs_filter].values

    # Perform layer remapping on the working adata (copy_input=True isolates caller's object).
    if (
        (spliced_layer != "spliced" or unspliced_layer != "unspliced")
        and spliced_layer in adata.layers
        and unspliced_layer in adata.layers
    ):
        adata.layers["spliced"] = _matrix_copy(adata.layers[spliced_layer])
        adata.layers["unspliced"] = _matrix_copy(adata.layers[unspliced_layer])
        logger.info(
            "Layer remapping applied: '%s' → 'spliced', '%s' → 'unspliced' (internal use only)",
            spliced_layer,
            unspliced_layer,
        )

    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError("Both 'spliced' and 'unspliced' layers are required after layer handling.")

    _warn_if_negative_layer_values(adata.layers["unspliced"], "unspliced")
    _warn_if_negative_layer_values(adata.layers["spliced"], "spliced")

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

    # gene_length > 0: length 0 is a sentinel for missing annotation (pp_bias fillna/empty
    # exons). log1p(0)=0 is an extreme x-leverage point for Huber (robust to y, not x)
    # and would bias the slope / all residuals. intron_number == 0 is biologically valid.
    valid_feat = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length > 0)
        & (intron_number >= 0)
    )
    n_zero_len = int(np.sum(np.isfinite(gene_length) & (gene_length <= 0)))
    n_nan_len = int(np.sum(~np.isfinite(gene_length)))
    if n_zero_len > 0 or n_nan_len > 0:
        logger.warning(
            "gene_length: %d gene(s) with length <= 0 and %d with missing/non-finite length "
            "are excluded from Huber bias-correction fit (valid_feat requires gene_length > 0). "
            "Their residuals use median centering when expressed. Check gene-feature mapping "
            "if this fraction is large (n_valid_feat=%d / %d).",
            n_zero_len,
            n_nan_len,
            int(valid_feat.sum()),
            int(adata.n_vars),
        )

    # ==================== PSEUDOBULK (optional) ====================
    is_pseudobulk = False
    if use_pseudobulk:
        if sample_col is None:
            raise ValueError("sample_col must be provided when use_pseudobulk=True")
        logger.info("Performing pseudobulk aggregation...")
        # Aggregate velocity layers; also carry counts when present so PyDESeq2 can
        # use raw counts after aggregation (not log1p .X).
        pb_layers = [ly for ly in ("spliced", "unspliced", "counts") if ly in adata.layers]
        x_layer_eff = pb_x_layer if pb_x_layer != "X" else None
        if (
            x_layer_eff is None
            and not pb_use_total_for_x
            and pseudobulk_de_backend == "pydeseq2"
            and "counts" in adata.layers
        ):
            x_layer_eff = "counts"
            logger.info(
                "Pseudobulk: aggregating layers['counts'] into .X for PyDESeq2 "
                "(adata.X may be log-normalized)."
            )
        adata = _pseudobulk_with_layers(
            adata,
            sample_col,
            groupby,
            layers=pb_layers,
            x_layer=x_layer_eff,
            use_total_for_x=pb_use_total_for_x,
            min_cells=min_cells,
            min_counts=min_counts,
        )
        is_pseudobulk = True
        with warnings.catch_warnings():
            warnings.simplefilter(
                "ignore", category=UserWarning
            )  # pandas/ann implicit index str conversion is benign here
            adata.obs[groupby] = _as_contrast_categorical(
                adata.obs[groupby], reference_group, target_group
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
        paired_replicates=paired_replicates,
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

    # Surface MixedLM / Memento extras when present
    for extra_col in [
        "mixedlm_coef",
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
    mixedlm_logfc_method = None
    n_mixed_sign_discordant = 0
    if use_mixed_model:
        n_mixed_failed = int(
            de_df.attrs.get("n_genes_failed_fit", 0) if hasattr(de_df, "attrs") else 0
        )
        mixed_failed_rate = float(
            de_df.attrs.get("failed_fit_rate", 0.0) if hasattr(de_df, "attrs") else 0.0
        )
        if hasattr(de_df, "attrs"):
            mixedlm_logfc_method = de_df.attrs.get("logFC_method")
            n_mixed_sign_discordant = int(
                de_df.attrs.get("n_genes_logFC_mixedlm_sign_discordant", 0) or 0
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
    # Placeholders; every mode branch must overwrite gamma_ref / gamma_info.
    # Unknown modes are rejected early (mode not in {heuristic, advanced}); the
    # final else below is a defensive guard if a new mode is added without wiring.
    gamma_ref = np.full(adata.n_vars, np.nan)
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
        adata_comp.layers["unspliced"] = _matrix_copy(uns_layer)
        adata_comp.layers["spliced"] = _matrix_copy(spl_layer)
        # Precomputed Mu/Ms were almost always built on raw-depth U/S; after rewriting
        # layers to size-factor-normalized values they are inconsistent. Recompute.
        use_precomputed_eff = advanced_use_precomputed
        if advanced_use_precomputed and ("Mu" in adata_comp.layers and "Ms" in adata_comp.layers):
            logger.warning(
                "mode='advanced' with advanced_use_precomputed=True: existing Mu/Ms may "
                "have been computed on raw-depth layers, but scATrans rewrote unspliced/"
                "spliced with size-factor normalization. Recomputing moments from the "
                "normalized layers for consistency. Pass advanced_use_precomputed=False "
                "to silence this message."
            )
            use_precomputed_eff = False
            # Drop stale moments so scVelo recomputes cleanly
            for _ly in ("Mu", "Ms"):
                if _ly in adata_comp.layers:
                    del adata_comp.layers[_ly]

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
                    use_precomputed=use_precomputed_eff,
                    recompute_neighbors=advanced_recompute_neighbors,
                    random_state=random_seed,
                )
            )
            velocity_source = "scvelo_moments_groupwise_ratio"
            velocity_layer_for_perm_uns = _matrix_copy(adata_comp.layers["Mu"])
            velocity_layer_for_perm_spl = _matrix_copy(adata_comp.layers["Ms"])
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
    else:
        # Early validation should already have rejected unknown modes; fail loudly
        # rather than silently keeping all-NaN gamma_ref if a new mode is added.
        raise AssertionError(
            f"Unhandled mode={mode!r}; expected 'heuristic' or 'advanced'. "
            "Add an explicit branch (and overwrite gamma_ref/gamma_info) for new modes."
        )

    # ==================== BIAS CORRECTION (now uses shared implementation) ====================
    total_us_raw = np.nan_to_num(
        _matrix_sum_axis0(uns_layer_raw) + _matrix_sum_axis0(spl_layer_raw)
    )

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

    pydeseq2_diag = {}
    if is_pseudobulk and pseudobulk_de_backend == "pydeseq2" and hasattr(de_df, "attrs"):
        pydeseq2_diag = {
            "used": True,
            "n_genes_filtered_low_count": int(de_df.attrs.get("n_genes_filtered_low_count", 0)),
            "n_genes_nan_from_deseq2": int(de_df.attrs.get("n_genes_nan_from_deseq2", 0)),
            "neutral_fill": bool(de_df.attrs.get("pydeseq2_neutral_fill", True)),
            "note": (
                "Genes filtered by min_counts or marked NaN by DESeq2 independent filtering "
                "appear as logFC=0, p_adj=1 and are not 'tested and non-significant'."
            ),
        }
    n_memento_not_returned = (
        int(de_df.attrs.get("n_genes_not_returned_by_memento", 0))
        if (use_memento_de and hasattr(de_df, "attrs"))
        else 0
    )

    diagnostics: dict[str, Any] = {
        "n_cells": int(adata.n_obs),
        "n_genes_input": int(adata.n_vars),
        "n_genes_with_valid_features": int(valid_feat.sum()),
        "unspliced_global_fraction": float(unspliced_fraction)
        if unspliced_fraction is not None
        else np.nan,
        "bias_correction": bias_info,
        "velocity": velocity_diag,
        # Full diagnose_design payload (warnings/recommendations) when auto-run
        "design": design_diag,
        "mixed_model": {
            "used": bool(use_mixed_model),
            "sample_col": sample_col if use_mixed_model else None,
            "paired_replicates": paired_replicates if use_mixed_model else None,
            "n_samples": int(adata.obs[sample_col].nunique())
            if (use_mixed_model and sample_col and sample_col in adata.obs.columns)
            else None,
            "mixedlm_grouping": (
                de_df.attrs.get("mixedlm_grouping") if hasattr(de_df, "attrs") else None
            )
            if use_mixed_model
            else None,
            "delta_variance_available": "delta_variance" in adata.var.columns,
            "median_delta_variance": float(np.nanmedian(adata.var["delta_variance"]))
            if "delta_variance" in adata.var.columns
            else np.nan,
            "n_genes_failed_fit": n_mixed_failed if use_mixed_model else 0,
            "failed_fit_rate": mixed_failed_rate if use_mixed_model else 0.0,
            "logFC_method": mixedlm_logfc_method if use_mixed_model else None,
            "n_genes_logFC_mixedlm_sign_discordant": (
                n_mixed_sign_discordant if use_mixed_model else 0
            ),
            "note": (
                "Lightweight LMM analogue (log1p + Wald/LRT); not NB-GLMM/voom. "
                "logFC is sample-aware mean-of-means log2FC; p_val tests mixedlm_coef. "
                "Inspect failed_fit_rate and n_genes_logFC_mixedlm_sign_discordant "
                "before publication claims."
                if use_mixed_model
                else None
            ),
        },
        "pydeseq2": pydeseq2_diag or {"used": False},
        "memento": {
            "used": bool(use_memento_de),
            "n_genes_not_returned": n_memento_not_returned,
            "note": (
                "Genes dropped by memento internal filters appear as logFC=0, p_adj=1 "
                "after reindexing and were not tested."
                if use_memento_de
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
    logfc_vals = adata.var["logFC"].values
    mixedlm_coef_vals = (
        adata.var["mixedlm_coef"].values if "mixedlm_coef" in adata.var.columns else None
    )
    direction_effect = _score_direction_effect(logfc_vals, mixedlm_coef=mixedlm_coef_vals)
    lambda_fc = max(_get_exponential_scale_lambda(logfc_vals), 0.25)
    lambda_res = max(_get_exponential_scale_lambda(residual), 1e-8)
    # Scale p-value leg on direction-positive genes only (avoids down-regulated
    # extreme p-values inflating λ and shrinking s3 for true up genes).
    lambda_pval = _lambda_pval_for_active_score(
        adata.var["p_adj"].values, direction_effect, floor=1.0
    )

    s1, s2, s3 = _composite_active_score_terms(
        logfc_vals,
        residual,
        adata.var["p_adj"].values,
        lambda_fc,
        lambda_res,
        lambda_pval,
        direction_effect=direction_effect,
    )

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
            if use_mixed_model:
                logger.warning(
                    "use_mixed_model=True with perm_de_backend='fast': observed active_score "
                    "uses MixedLM logFC/p_adj, but the permutation null uses a naive "
                    "t-test (no (1|sample) random effect). active_score_pval / "
                    "active_score_fdr are therefore **not valid** for the MixedLM-based "
                    "composite score. Prefer perm_de_backend='same' (null refits MixedLM) "
                    "for reporting, or treat only unspliced_excess_fdr as the permutation "
                    "significance (velocity residual null does not use the DE backend)."
                )
        elif perm_de_backend == "same":
            perm_pb_backend, perm_de_method = pseudobulk_de_backend, de_method
            if use_mixed_model:
                logger.info(
                    "perm_de_backend='same': permutations use MixedLM "
                    "(sample_col=%r, mixed_model_pval=%r, paired_replicates=%s) "
                    "with **sample/RE-cluster-level** label permutation "
                    "(not cell-level) so (1|sample) exchangeability holds — "
                    "matching observed DE for a valid active_score null.",
                    sample_col,
                    mixed_model_pval,
                    paired_replicates,
                )
            else:
                logger.info(
                    "perm_de_backend='same': permutations use the same DE backend as the "
                    "main analysis (%s / pseudobulk=%s).",
                    de_method,
                    is_pseudobulk,
                )
            n_genes_for_perm = int(np.sum(valid_feat)) if valid_feat is not None else 0
            if is_pseudobulk and pseudobulk_de_backend == "pydeseq2" and n_genes_for_perm > 5000:
                logger.warning(
                    "use_pseudobulk=True + perm_de_backend='same' (default) refits PyDESeq2 "
                    "from scratch on every permutation. With %d genes and n_perm=%d this can "
                    "take many minutes with no progress output. Consider perm_de_backend='fast' "
                    "(permutations use a scanpy t-test instead; much faster) if you don't need "
                    "the null distribution to use PyDESeq2 exactly.",
                    n_genes_for_perm,
                    n_perm,
                )
            if use_mixed_model and n_genes_for_perm > 2000 and n_perm > 50:
                logger.warning(
                    "use_mixed_model=True + perm_de_backend='same' refits a MixedLM per gene "
                    "on every permutation (%d genes × n_perm=%d). This can be very slow; "
                    "use n_jobs for parallelism, or perm_de_backend='fast' only if you accept "
                    "that active_score_fdr will not match the MixedLM estimator.",
                    n_genes_for_perm,
                    n_perm,
                )
        else:
            raise ValueError("perm_de_backend must be 'fast' or 'same'")

        # MixedLM null: only when perm_de_backend='same' (mirrors "same estimator" contract).
        # 'fast' deliberately uses scanpy t-test (warned above).
        perm_mixed_model = bool(use_mixed_model and perm_de_backend == "same")

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
            and current_max_perm is not None
            and np.isfinite(current_max_perm)
            and current_max_perm < n_perm
        ):
            n_perm = int(current_max_perm)

        # MixedLM sample/cluster-level null: exact permutation space is over
        # RE clusters (≈ biological samples), not cells. Cap n_perm for unpaired.
        if (
            use_mixed_model
            and perm_mixed_model
            and auto_adjust_n_perm
            and sample_col is not None
            and sample_col in adata.obs.columns
            and not paired_replicates
        ):
            try:
                from .._de import _resolve_mixedlm_random_groups

                _obs_tmp = pd.DataFrame(
                    {
                        "condition": np.asarray(obs_labels).astype(str),
                        "sample": adata.obs[sample_col].astype(str).to_numpy(),
                    }
                )
                _g_ids, _ = _resolve_mixedlm_random_groups(
                    _obs_tmp,
                    "condition",
                    "sample",
                    paired_replicates=False,
                    quiet=True,
                )
                _g_ids = np.asarray(_g_ids).astype(str)
                _labs = np.asarray(obs_labels).astype(str)
                # One condition per cluster → count clusters per arm
                _cluster_cond: dict[str, str] = {}
                for _gid in pd.unique(_g_ids):
                    _m = _g_ids == _gid
                    _u, _c = np.unique(_labs[_m], return_counts=True)
                    _cluster_cond[str(_gid)] = str(_u[int(np.argmax(_c))])
                _n_t_cl = sum(1 for _v in _cluster_cond.values() if _v == str(target_group))
                _n_r_cl = sum(1 for _v in _cluster_cond.values() if _v == str(reference_group))
                _n_cl = _n_t_cl + _n_r_cl
                if _n_t_cl > 0 and _n_r_cl > 0 and _n_cl <= 30:
                    _max_sample_perm = max(1, int(comb(_n_cl, _n_t_cl) - 1))
                    if _max_sample_perm < n_perm:
                        logger.info(
                            "MixedLM sample-level permutation space is only %d distinct "
                            "assignments (%d clusters: %d target / %d ref); reducing "
                            "n_perm from %d to %d (auto_adjust_n_perm=True).",
                            _max_sample_perm,
                            _n_cl,
                            _n_t_cl,
                            _n_r_cl,
                            n_perm,
                            _max_sample_perm,
                        )
                        n_perm = _max_sample_perm
                        current_max_perm = float(_max_sample_perm)
            except Exception as e:
                logger.debug("MixedLM n_perm sample-space cap skipped: %s", e)

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
                use_mixed_model=perm_mixed_model,
                sample_col=sample_col if perm_mixed_model else None,
                mixed_model_pval=mixed_model_pval,
                paired_replicates=paired_replicates if perm_mixed_model else False,
                valid_expr=valid_expr,
                min_counts_per_gene=pydeseq2_min_counts,
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
        sample_col=sample_col if (use_mixed_model or is_pseudobulk) else None,
        mixed_model_pval=mixed_model_pval if use_mixed_model else None,
        paired_replicates=paired_replicates if use_mixed_model else None,
        perm_use_mixed_model=(
            bool(use_mixed_model and perm_de_backend == "same") if use_permutation else None
        ),
        prioritize_velocity=prioritize_velocity,
        ranking_mode=ranking_mode,
        weight_fc=weight_fc,
        weight_unspliced=weight_unspliced,
        weight_pval=weight_pval,
    )


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
    if "mixedlm_coef" in adata.var.columns:
        cols.append("mixedlm_coef")
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
    # Residual / score cutoffs scale with single-cell vs pseudobulk (see PSEUDOBULK_FILTER_DEFAULTS).
    scale_defaults = (
        PSEUDOBULK_FILTER_DEFAULTS
        if extra_metadata.get("is_pseudobulk")
        else HEURISTIC_FILTER_DEFAULTS
    )
    ue_fdr_cutoff = extra_metadata.get(
        "unspliced_excess_fdr_cutoff",
        scale_defaults["unspliced_excess_fdr_cutoff"],
    )
    filter_context: dict[str, Any] = {
        "use_permutation": bool(use_permutation),
        "use_fdr_for_significance": bool(extra_metadata.get("use_fdr_for_significance", True))
        if use_permutation
        else False,
        "perm_disabled_reason": extra_metadata.get("perm_disabled_reason"),
        "is_pseudobulk": bool(extra_metadata.get("is_pseudobulk")),
        "pval_cutoff": extra_metadata.get("pval_cutoff", scale_defaults["pval_cutoff"]),
        "logfc_cutoff": extra_metadata.get("logfc_cutoff", scale_defaults["logfc_cutoff"]),
        "unspliced_excess_fdr_cutoff": ue_fdr_cutoff,
        "use_delta_variance_pval": bool(extra_metadata.get("use_delta_variance_pval")),
        "delta_var_pval_cutoff": extra_metadata.get("delta_var_pval_cutoff", 0.05),
    }

    if use_permutation and UNSPLICED_EXCESS_FDR_COL in adata.var.columns:
        if not extra_metadata.get("use_fdr_for_significance", True):
            reason = extra_metadata.get("perm_disabled_reason", "small_permutation_space")
            logger.warning(
                "Permutation space is very small (%s); unspliced_excess_fdr and active_score_fdr "
                "were not applied to the built-in significant list. Use filter_active_genes("
                "preset='significant') on all_results to reproduce this list, or preset="
                "'pseudobulk'/'heuristic' for exploratory cutoffs.",
                reason,
            )
        mask = _builtin_significant_mask(
            adata.var,
            use_permutation=True,
            extra_metadata=filter_context,
        )
    else:
        mask = pd.Series(False, index=adata.var.index)
        if not use_permutation:
            logger.warning(
                "Built-in 'significant' list is empty because use_permutation=False "
                "(permutation FDR is required for that strict list). "
                "This is expected — significant is intentionally strict and often empty. "
                "For exploratory gene lists use: "
                "filter_active_genes(all_results, preset='heuristic')  # or 'pseudobulk' / "
                "'permissive', or pass custom residual/score/FDR cutoffs. "
                "To populate significant, re-run with use_permutation=True (and enough n_perm)."
            )

    var_df = _as_var_dataframe(adata)
    significant = var_df.loc[mask, cols].copy().sort_values("active_score", ascending=False)
    all_results = var_df.loc[:, cols].copy().sort_values("active_score", ascending=False)
    all_results.attrs["scatrans_filter_context"] = filter_context

    logger.info(
        "Analysis completed in %s mode! Significant active genes: %d", mode, len(significant)
    )
    if len(significant) == 0 and use_permutation:
        logger.warning(
            "No genes passed the built-in significant thresholds (strict AND of DE + residual + "
            "score + permutation FDR). This is common on real data. "
            "Use filter_active_genes(all_results, preset=%r) or custom cutoffs for exploration; "
            "inspect all_results score/residual/FDR distributions before claiming no signal.",
            "pseudobulk" if extra_metadata.get("is_pseudobulk") else "heuristic",
        )

    # --- Rich metadata (merge to protect raw_gene_list etc.) ---
    # Both cell-level and pseudobulk paths size-factor-normalize U/S before delta.
    velocity_delta_layer = (
        "scvelo_Mu_Ms_moments"
        if velocity_source.startswith("scvelo_moments")
        else "size_factor_normalized_spliced_unspliced"
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
            "logFC": (f"> {extra_metadata.get('logfc_cutoff', scale_defaults['logfc_cutoff'])}"),
            "p_adj": (f"< {extra_metadata.get('pval_cutoff', scale_defaults['pval_cutoff'])}"),
            "unspliced_excess_residual": (
                f"> {scale_defaults['unspliced_excess_residual_cutoff']}"
            ),
            "active_score": f">= {scale_defaults['active_score_cutoff']}",
            "active_score_fdr": f"< {scale_defaults['active_score_fdr_cutoff']}",
            "unspliced_excess_fdr": (
                f"< {ue_fdr_cutoff} (requires use_permutation=True)"
                if use_permutation
                else "not evaluated (use_permutation=False)"
            ),
            "scale": "pseudobulk" if extra_metadata.get("is_pseudobulk") else "heuristic",
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
        "mixed_model_pval": extra_metadata.get("mixed_model_pval"),
        "paired_replicates": extra_metadata.get("paired_replicates"),
        "perm_use_mixed_model": extra_metadata.get("perm_use_mixed_model"),
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
            "Permutation null: unspliced/spliced (or Mu/Ms) layers are held fixed from the "
            "observed run; group labels are shuffled. For each shuffle we recompute DE, "
            "reference gamma (or EB posterior given a fixed prior), unspliced excess residual, "
            "composite active_score, and unspliced_excess permutation p-values. "
            "Layers are not re-smoothed; gamma is not frozen at the observed per-gene values."
        )
        if extra_metadata.get("use_memento_de") and not extra_metadata.get("perm_use_memento_de"):
            note += (
                " Memento was used for the observed DE; permutations used the configured fast/same "
                "non-Memento DE backend unless perm_use_memento_de=True."
            )
        elif extra_metadata.get("use_memento_de") and extra_metadata.get("perm_use_memento_de"):
            note += " Memento was used for both observed DE and permutation null."
        if extra_metadata.get("use_mixed_model") and extra_metadata.get("perm_use_mixed_model"):
            note += (
                " MixedLM (condition + (1|sample)) was used for both observed DE and the "
                "permutation null (perm_de_backend='same'). Labels are shuffled at the "
                "sample / random-effect cluster level (not independently per cell) so the "
                "null preserves hierarchical exchangeability; active_score_pval/fdr match "
                "the MixedLM-based composite score."
            )
        elif extra_metadata.get("use_mixed_model") and not extra_metadata.get(
            "perm_use_mixed_model"
        ):
            note += (
                " MixedLM was used for observed DE, but permutations used a non-MixedLM DE "
                "backend (typically perm_de_backend='fast' → scanpy t-test). "
                "active_score_pval/fdr are not a valid null for the MixedLM composite; "
                "unspliced_excess_fdr still uses a matching residual estimator."
            )
        meta["permutation_approximation_note"] = note

    from .._utils import _merge_scatrans_uns

    adata.uns["scatrans"] = _merge_scatrans_uns(existing, meta)

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
                "Permutation used %d iterations (layers fixed; gamma/DE/residual recomputed "
                "under shuffled labels).",
                n_perm,
            )
    except Exception:
        pass

    # --- Optional plotting (delegated) ---
    if show_plot:
        try:
            from .. import pl

            pl.comet_plot(
                all_results,
                top_n=12,
                title=f"scATrans Active Drivers ({display_mode})",
            )
        except Exception:
            logger.debug("show_plot=True but plotting failed (missing optional deps or display).")

    return adata, significant, all_results
