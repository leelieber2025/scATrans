"""
Permutation testing support for active_score significance.

The heavy _single_permutation_task (and the orchestration) was one of the
biggest contributors to tl.py line count. Extracted here.
"""

from __future__ import annotations

import logging
import warnings
from typing import Any

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests

from ._de import _run_de_wrapper

# local import to avoid circulars at module load
from ._utils import (
    _apply_de_preprocess,
    _fit_huber_bias_correction,
    _is_integer_counts_like,
    _resolve_aligned_raw_counts,
    _soft_scale,
)
from ._velocity import _compute_velocity_delta

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _shuffle_condition_labels(
    original_labels: np.ndarray,
    *,
    rng: np.random.Generator,
    use_mixed_model: bool = False,
    sample_ids: np.ndarray | None = None,
    paired_replicates: bool = False,
    max_attempts: int = 50,
) -> tuple[np.ndarray, np.ndarray | None] | None:
    """Draw one label permutation under the exchangeability of the DE backend.

    **Non-MixedLM** (scanpy / Memento / pseudobulk re-aggregation): permute
    condition labels at the **cell** level. Pseudobulk nulls re-aggregate under
    the new labels; cell-level i.i.d. exchangeability is the usual null.

    **MixedLM + unpaired** (``paired_replicates=False``): permute conditions at
    the **random-effect cluster** level — the same units
    :func:`scatrans._de._resolve_mixedlm_random_groups` would use under the
    *observed* labels (``sample_col`` when IDs are condition-pure, or
    ``condition::sample`` composite when sample strings are reused across
    conditions). All cells in a cluster share one condition after the shuffle.
    The second return value is the observed RE cluster ID per cell; the caller
    should use it as ``sample_col`` for the null MixedLM fit so cluster
    membership (and ``n_groups``) cannot collapse when recycled sample strings
    would otherwise merge under the new labels.

    **MixedLM + paired** (``paired_replicates=True``): within each ``sample_id``,
    permute cell labels (preserves within-subject condition counts and
    ``(1|sample)`` structure). Second return value is the original sample IDs.

    Returns
    -------
    (shuffled_labels, sample_ids_for_mixedlm) or None
        ``sample_ids_for_mixedlm`` is ``None`` for non-MixedLM shuffles.
    """
    labels = np.asarray(original_labels)
    if labels.size == 0:
        return None

    if not use_mixed_model or sample_ids is None:
        for _ in range(max_attempts):
            shuffled = rng.permutation(labels)
            if not np.array_equal(shuffled, labels):
                return shuffled, None
        return None

    samples = np.asarray(sample_ids)
    if samples.shape[0] != labels.shape[0]:
        raise ValueError(
            f"sample_ids length ({samples.shape[0]}) must match labels ({labels.shape[0]})"
        )

    if paired_replicates:
        return _shuffle_labels_within_samples(labels, samples, rng, max_attempts=max_attempts)

    return _shuffle_labels_by_mixedlm_clusters(labels, samples, rng, max_attempts=max_attempts)


def _shuffle_labels_within_samples(
    labels: np.ndarray,
    samples: np.ndarray,
    rng: np.random.Generator,
    *,
    max_attempts: int = 50,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Paired MixedLM null: permute labels within each sample_id."""
    labels = np.asarray(labels)
    samples = np.asarray(samples)
    for _ in range(max_attempts):
        out = labels.copy()
        changed = False
        for s in pd.unique(samples.astype(str)):
            mask = samples.astype(str) == str(s)
            block = labels[mask]
            if block.size <= 1 or len(np.unique(block)) <= 1:
                continue
            perm_block = rng.permutation(block)
            out[mask] = perm_block
            if not np.array_equal(perm_block, block):
                changed = True
        if changed and not np.array_equal(out, labels):
            return out, samples
    return None


def _shuffle_labels_by_mixedlm_clusters(
    labels: np.ndarray,
    samples: np.ndarray,
    rng: np.random.Generator,
    *,
    max_attempts: int = 50,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Unpaired MixedLM null: reassign conditions to whole RE clusters.

    Clusters match the observed MixedLM random-effect units. Returns the
    observed cluster IDs so the null MixedLM can use them as ``sample_col``
    (stable ``n_groups``, pure within-cluster conditions).
    """
    # Local import avoids a hard cycle at module import time.
    from ._de import _resolve_mixedlm_random_groups

    labels = np.asarray(labels)
    samples = np.asarray(samples)
    obs = pd.DataFrame(
        {
            "condition": pd.Series(labels).astype(str).to_numpy(),
            "sample": pd.Series(samples).astype(str).to_numpy(),
        }
    )
    groups, _meta = _resolve_mixedlm_random_groups(
        obs,
        "condition",
        "sample",
        paired_replicates=False,
        quiet=True,
    )
    groups = np.asarray(groups).astype(str)

    # One condition per cluster under the observed design (composite or pure sample).
    # If a cluster is impure (pathological), split into (group, label) sub-blocks.
    cluster_keys = groups.copy()
    for g in pd.unique(groups):
        m = groups == g
        labs = np.unique(labels[m].astype(str))
        if labs.size > 1:
            cluster_keys[m] = np.array(
                [f"{g}||{lab}" for lab in labels[m].astype(str)],
                dtype=object,
            ).astype(str)

    uniq_clusters = pd.unique(cluster_keys)
    cluster_conds: list[str] = []
    for c in uniq_clusters:
        m = cluster_keys == c
        labs, cnts = np.unique(labels[m].astype(str), return_counts=True)
        cluster_conds.append(str(labs[int(np.argmax(cnts))]))
    cluster_conds_arr = np.asarray(cluster_conds, dtype=object)

    if len(uniq_clusters) < 2 or len(np.unique(cluster_conds_arr)) < 2:
        # Degenerate: fall back to cell-level rather than returning the identity forever.
        for _ in range(max_attempts):
            shuffled = rng.permutation(labels)
            if not np.array_equal(shuffled, labels):
                return shuffled, cluster_keys
        return None

    for _ in range(max_attempts):
        perm_conds = rng.permutation(cluster_conds_arr)
        if np.array_equal(perm_conds, cluster_conds_arr):
            continue
        mapping = dict(zip(uniq_clusters, perm_conds))
        shuffled = np.array([mapping[c] for c in cluster_keys], dtype=object)
        if labels.dtype.kind in {"U", "S", "O"}:
            try:
                shuffled = shuffled.astype(labels.dtype, copy=False)
            except (TypeError, ValueError):
                shuffled = shuffled.astype(str)
        if not np.array_equal(shuffled, labels):
            return shuffled, cluster_keys
    return None


def _compute_perm_velocity_delta(
    *,
    velocity_source: str,
    uns_layer: Any,
    spl_layer: Any,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float,
    gamma_method: str,
    eb_prior: dict[str, Any] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Velocity delta for one permutation replicate.

    Both tracks deliberately share the same groupwise ratio estimator
    (:func:`_compute_velocity_delta`). The difference is only in the **layers**
    the caller passes:

    - ``scvelo_moments*``: pre-smoothed Mu/Ms from the main analysis (label
      permutation keeps smoothed values fixed).
    - heuristic / fallback: library-size-normalized unspliced/spliced layers.

    ``velocity_source`` is accepted for API/diagnostics symmetry with the main
    analysis; it does **not** change the estimator. If a future track needs a
    different smoother or estimator, branch on ``velocity_source`` here (or
    pass distinct layers upstream) — do not reintroduce a no-op if/else.
    """
    _ = velocity_source  # reserved for future track-specific estimators
    return _compute_velocity_delta(
        uns_layer,
        spl_layer,
        t_mask,
        r_mask,
        prior_weight,
        gamma_method=gamma_method,
        eb_prior=eb_prior,
    )


def _single_permutation_task(
    seed: int,
    original_labels: np.ndarray,
    target_group: str,
    reference_group: str,
    adata_subset: Any,
    X_features: np.ndarray | None,
    valid_feat: np.ndarray,
    uns_layer: Any,
    spl_layer: Any,
    total_us_for_filter: np.ndarray,
    min_total_counts: int,
    weight_fc: float,
    weight_unspliced: float,
    weight_pval: float,
    lambda_fc: float,
    lambda_res: float,
    lambda_pval: float,
    is_pseudobulk: bool,
    pb_backend: str,
    de_method: str,
    prior_weight: float,
    gamma_method: str,
    de_preprocess: str,
    strict_pydeseq2_counts: bool,
    eb_prior: dict[str, Any] | None = None,
    bias_correction: str = "huber_length_intron",
    # Memento support for permutation (advanced, usually False for speed)
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    # MixedLM support for permutation (required for valid active_score null when
    # observed DE used use_mixed_model=True + perm_de_backend='same')
    use_mixed_model: bool = False,
    sample_col: str | None = None,
    mixed_model_pval: str = "wald",
    paired_replicates: bool = False,
    velocity_source: str = "heuristic_global_ratio",
    min_counts_per_gene: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """One permutation replicate.

    Returns (active_score_vector, unspliced_excess_residual_vector) for that shuffle.

    bias_correction is forwarded to the shared bias correction routine so that
    permutation scores are computed under the same correction setting the user chose
    for the real data (default = on).

    Count-based DE (PyDESeq2 / Memento) re-resolves ``layers['counts']`` on the
    per-task copy so the null uses the **same count matrix** as the observed DE
    (not post-aggregation ``.X`` which may be spliced+unspliced).

    When ``use_mixed_model=True``, each null DE uses the same LMM
    (``y ~ C(condition) + (1|sample)``) as the observed run, and labels are
    shuffled at the **sample / RE-cluster** level (not independently per cell)
    so the null preserves exchangeability under the hierarchical design.
    """
    rng = np.random.default_rng(seed)
    sample_ids = None
    if use_mixed_model and sample_col is not None:
        if sample_col not in adata_subset.obs.columns:
            logger.warning(
                "use_mixed_model permutation: sample_col=%r missing from obs; "
                "falling back to cell-level label shuffle (invalid MixedLM null).",
                sample_col,
            )
        else:
            sample_ids = adata_subset.obs[sample_col].to_numpy()

    shuffle_result = _shuffle_condition_labels(
        original_labels,
        rng=rng,
        use_mixed_model=use_mixed_model,
        sample_ids=sample_ids,
        paired_replicates=paired_replicates if use_mixed_model else False,
        max_attempts=50,
    )
    if shuffle_result is None:
        logger.warning(
            "Failed to generate a different permutation after 50 attempts. Skipping this replicate."
        )
        # Length must match n_genes (success path), not n_obs — failed reps are
        # discarded via NaN checks, but wrong length would break future partial-keep logic.
        n_genes = int(adata_subset.n_vars)
        return np.full(n_genes, np.nan), np.full(n_genes, np.nan)
    shuffled_labels, sample_ids_for_mixedlm = shuffle_result

    ad_temp = adata_subset.copy()
    # Pin MixedLM RE membership to the *observed* clusters so null n_groups /
    # purity match the observed design (critical when sample strings recycle).
    if (
        use_mixed_model
        and sample_col is not None
        and sample_ids_for_mixedlm is not None
        and sample_col in ad_temp.obs.columns
    ):
        ad_temp.obs[sample_col] = np.asarray(sample_ids_for_mixedlm).astype(str)

    # Preprocess for the *perm* DE backend scale:
    # - PyDESeq2: leave counts untransformed (skip_auto).
    # - Scanpy (incl. perm_de_backend='fast' while main used PyDESeq2): .X is often
    #   still integer pseudobulk sums; rank_genes_groups needs log1p or logFC
    #   overflows (expm1 on large counts). Apply normalize_log1p only when .X
    #   looks like raw counts so we never double-log already-logged matrices.
    if is_pseudobulk and pb_backend == "pydeseq2":
        _apply_de_preprocess(ad_temp, de_preprocess, skip_auto=True)
    elif is_pseudobulk and pb_backend == "scanpy" and _is_integer_counts_like(ad_temp.X):
        _apply_de_preprocess(ad_temp, "normalize_log1p", skip_auto=False)
    else:
        _apply_de_preprocess(
            ad_temp,
            de_preprocess,
            skip_auto=is_pseudobulk and pb_backend == "pydeseq2",
        )

    # Match active_score observed DE: prefer aligned layers['counts'] for count backends.
    perm_counts = None
    if use_memento_de or (is_pseudobulk and pb_backend == "pydeseq2"):
        perm_counts = _resolve_aligned_raw_counts(ad_temp, layer="counts", require_integer=True)

    perm_de_df = _run_de_wrapper(
        ad_temp,
        groupby="_unused_when_labels_provided",
        target_group=target_group,
        reference_group=reference_group,
        de_method=de_method,
        is_pseudobulk=is_pseudobulk,
        pb_backend=pb_backend,
        n_jobs=1,
        labels=shuffled_labels,
        strict_pydeseq2_counts=strict_pydeseq2_counts,
        use_mixed_model=use_mixed_model,
        sample_col=sample_col if use_mixed_model else None,
        mixed_model_pval=mixed_model_pval,
        paired_replicates=paired_replicates if use_mixed_model else False,
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        memento_n_cpus=memento_n_cpus,
        counts=perm_counts,
        min_counts_per_gene=min_counts_per_gene,
    )

    t_mask = shuffled_labels == target_group
    r_mask = shuffled_labels == reference_group
    delta_velocity, _, _gamma_ref, _ = _compute_perm_velocity_delta(
        velocity_source=velocity_source,
        uns_layer=uns_layer,
        spl_layer=spl_layer,
        t_mask=t_mask,
        r_mask=r_mask,
        prior_weight=prior_weight,
        gamma_method=gamma_method,
        eb_prior=eb_prior,
    )

    total_us_for_filter = np.asarray(total_us_for_filter)
    valid_expr = total_us_for_filter >= min_total_counts

    # Use the shared bias correction (DRY). bias_correction setting is respected
    # so that permuted scores are comparable to the real run.
    gene_length = adata_subset.var["gene_length"].values
    intron_number = adata_subset.var["intron_number"].values

    residual, _bias_info = _fit_huber_bias_correction(
        delta_velocity,
        gene_length,
        intron_number,
        total_us_for_filter,
        valid_feat,
        valid_expr,
        X_features,
        bias_correction=bias_correction,
    )

    s1 = _soft_scale(perm_de_df["logFC"].values, lambda_fc)
    s2 = _soft_scale(residual, lambda_res)
    s3 = _soft_scale(-np.log10(perm_de_df["p_adj"].values + 1e-300), lambda_pval)

    total_w = weight_fc + weight_unspliced + weight_pval
    perm_score = (weight_fc * s1 + weight_unspliced * s2 + weight_pval * s3) / total_w * 100.0
    return perm_score, residual


def run_permutation_test(
    *,
    n_perm: int,
    effective_n_jobs: int,
    random_seed: int,
    obs_labels: np.ndarray,
    target_group: str,
    reference_group: str,
    adata: Any,
    X_features: np.ndarray | None,
    valid_feat: np.ndarray,
    velocity_layer_for_perm_uns: Any,
    velocity_layer_for_perm_spl: Any,
    total_us_raw: np.ndarray,
    min_total_counts: int,
    weight_fc: float,
    weight_unspliced: float,
    weight_pval: float,
    lambda_fc: float,
    lambda_res: float,
    lambda_pval: float,
    is_pseudobulk: bool,
    perm_pb_backend: str,
    perm_de_method: str,
    prior_weight: float,
    gamma_method: str,
    de_preprocess: str,
    strict_pydeseq2_counts: bool,
    real_score: np.ndarray,
    real_residual: np.ndarray,
    eb_prior: dict[str, Any] | None = None,
    velocity_source: str = "heuristic_global_ratio",
    bias_correction: str = "huber_length_intron",
    # Memento forwarding for advanced consistent permutation
    use_memento_de: bool = False,
    memento_capture_rate: float = 0.07,
    memento_num_boot: int = 5000,
    memento_n_cpus: int = -1,
    # MixedLM forwarding for consistent active_score null (perm_de_backend='same')
    use_mixed_model: bool = False,
    sample_col: str | None = None,
    mixed_model_pval: str = "wald",
    paired_replicates: bool = False,
    valid_expr: np.ndarray | None = None,
    min_counts_per_gene: int = 10,
) -> tuple:
    """Run parallel permutation and return score/residual p-values and FDR arrays.

    Returns
    -------
    active_score_pval, active_score_fdr,
    unspliced_excess_pval, unspliced_excess_fdr,
    use_fdr, disabled_reason
    """
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
                gamma_method,
                de_preprocess,
                strict_pydeseq2_counts,
                eb_prior,
                velocity_source=velocity_source,
                bias_correction=bias_correction,
                use_memento_de=use_memento_de,
                memento_capture_rate=memento_capture_rate,
                memento_num_boot=memento_num_boot,
                memento_n_cpus=memento_n_cpus,
                use_mixed_model=use_mixed_model,
                sample_col=sample_col,
                mixed_model_pval=mixed_model_pval,
                paired_replicates=paired_replicates,
                min_counts_per_gene=min_counts_per_gene,
            )
            for i in range(n_perm)
        )

    # Filter out any replicates that failed to produce a distinct shuffle (returned NaNs).
    # This prevents using original labels for the null (which would bias p-values upward).
    valid_scores = []
    valid_residuals = []
    for r in perm_results:
        s, res = r
        if s is not None and res is not None and not np.any(np.isnan(np.asarray(s))):
            valid_scores.append(s)
            valid_residuals.append(res)
    n_success = len(valid_scores)
    if n_success == 0:
        logger.warning(
            "All %d permutation replicates failed to generate distinct shuffles.", n_perm
        )
        # Fall back to non-informative p-values (no power to reject).
        n_genes = adata.n_vars
        return (
            np.ones(n_genes),
            np.ones(n_genes),
            np.ones(n_genes),
            np.ones(n_genes),
            False,
            "permutation_shuffle_failed",
        )
    perm_scores_matrix = np.vstack(valid_scores)
    perm_residual_matrix = np.vstack(valid_residuals)

    exceed_count = np.sum(perm_scores_matrix >= real_score.reshape(1, -1), axis=0)
    active_score_pval = (1.0 + exceed_count) / (n_success + 1.0)

    # One-sided test for positive unspliced excess (matches active-gene direction filter).
    exceed_res = np.sum(
        perm_residual_matrix >= np.asarray(real_residual, dtype=float).reshape(1, -1), axis=0
    )
    unspliced_excess_pval = (1.0 + exceed_res) / (n_success + 1.0)

    active_score_fdr = np.ones(adata.n_vars)
    unspliced_excess_fdr = np.ones(adata.n_vars)
    if valid_expr is None:
        valid_expr = adata.var.get("valid_expr", np.ones(adata.n_vars, dtype=bool))
    valid_expr = np.asarray(valid_expr, dtype=bool)
    if valid_expr.sum() > 0:
        active_score_fdr[valid_expr] = multipletests(
            active_score_pval[valid_expr], method="fdr_bh"
        )[1]
        unspliced_excess_fdr[valid_expr] = multipletests(
            unspliced_excess_pval[valid_expr], method="fdr_bh"
        )[1]

    # FDR decision and disabled_reason:
    # - FDR is computed on valid_expr (always, for the pvals we have).
    # - For very small permutation spaces (n_perm < 100, common in pseudobulk with few samples),
    #   we mark use_fdr=False and provide reason so callers can avoid using FDR for significance
    #   (p-values become coarse; BH adjustment less reliable). This eliminates vestigial logic in tl.py.
    use_fdr = n_success >= 100
    disabled_reason = None if use_fdr else "small_permutation_space"

    return (
        active_score_pval,
        active_score_fdr,
        unspliced_excess_pval,
        unspliced_excess_fdr,
        use_fdr,
        disabled_reason,
    )
