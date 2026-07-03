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
from joblib import Parallel, delayed
from statsmodels.stats.multitest import multipletests

from ._de import _run_de_wrapper

# local import to avoid circulars at module load
from ._utils import (
    _apply_de_preprocess,
    _fit_huber_bias_correction,
    _soft_scale,
)
from ._velocity import _compute_velocity_delta

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


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

    For ``scvelo_moments*`` sources the layers are pre-smoothed Mu/Ms from the main
    analysis (label permutation keeps smoothed values fixed). For heuristic sources
    they are library-size-normalized unspliced/spliced layers.
    """
    if velocity_source.startswith("scvelo_moments"):
        # Same groupwise ratio estimator as _compute_moments_velocity_delta (on Mu/Ms).
        return _compute_velocity_delta(
            uns_layer,
            spl_layer,
            t_mask,
            r_mask,
            prior_weight,
            gamma_method=gamma_method,
            eb_prior=eb_prior,
        )
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
    velocity_source: str = "heuristic_global_ratio",
) -> tuple[np.ndarray, np.ndarray]:
    """One permutation replicate.

    Returns (active_score_vector, unspliced_excess_residual_vector) for that shuffle.

    bias_correction is forwarded to the shared bias correction routine so that
    permutation scores are computed under the same correction setting the user chose
    for the real data (default = on).
    """
    rng = np.random.default_rng(seed)
    for _ in range(50):
        shuffled_labels = rng.permutation(original_labels)
        if not np.array_equal(shuffled_labels, original_labels):
            break
    else:
        logger.warning(
            "Failed to generate a different permutation after 50 attempts. Skipping this replicate."
        )
        n = len(original_labels)
        return np.full(n, np.nan), np.full(n, np.nan)

    ad_temp = adata_subset.copy()

    _apply_de_preprocess(
        ad_temp,
        de_preprocess,
        skip_auto=is_pseudobulk and pb_backend == "pydeseq2",
    )

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
        use_memento_de=use_memento_de,
        memento_capture_rate=memento_capture_rate,
        memento_num_boot=memento_num_boot,
        memento_n_cpus=memento_n_cpus,
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
    valid_expr: np.ndarray | None = None,
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
