"""
scATrans internal velocity delta computation (heuristic + advanced scVelo moments track).

Reference gamma (U/S ratio) estimation supports multiple robust methods, including
"empirical_bayes" which implements hierarchical gamma shrinkage: per-gene
log-ratios in the reference group are shrunk toward a robust data-driven prior
(estimated via trimmed median + MAD on log-ratios). This borrows strength across
genes and improves stability especially for small reference groups.

Extracted from original tl.py.
"""

from __future__ import annotations

import logging
from typing import Any

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

from ._utils import _get_group_mean, _matrix_shape

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

_EPS = 1e-8
_MIN_GENES_FOR_EB_PRIOR = 10
_DEFAULT_TRIM_FRACTION = 0.05


def _robust_mad_scale(x: np.ndarray) -> float:
    """MAD-based robust scale estimate (≈ std for Gaussian)."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 2:
        return 0.5
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return max(1e-3, float(mad * 1.4826))


def _estimate_eb_prior_from_reference(
    U_r: np.ndarray,
    S_r: np.ndarray,
    *,
    eps: float = _EPS,
    trim_fraction: float = _DEFAULT_TRIM_FRACTION,
    count_pseudocount: float = 1.0,
    min_genes_for_prior: int = _MIN_GENES_FOR_EB_PRIOR,
    prior_weight: float = 5.0,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """
    Estimate log-ratio empirical-Bayes prior hyperparameters from reference-group means.

    Uses median + MAD on trimmed log-ratios; falls back to global median when too few genes.
    """
    U_r = np.asarray(U_r, dtype=float).ravel()
    S_r = np.asarray(S_r, dtype=float).ravel()
    r_g = np.log((U_r + eps) / (S_r + eps))
    expressed = (U_r + S_r) > 0
    r_valid = r_g[expressed & np.isfinite(r_g)]

    fallback_triggered = False
    n_genes_used = int(r_valid.size)

    if n_genes_used < min_genes_for_prior:
        # Too few expressed genes for any per-gene ratio prior -> global sum ratio fallback.
        fallback_triggered = True
        global_ratio = (float(np.sum(U_r)) + eps) / (float(np.sum(S_r)) + eps)
        prior_mean_log = float(np.log(global_ratio))
        tau = max(0.25, _robust_mad_scale(r_valid) if r_valid.size else 0.5)
        method_detail = "empirical_bayes_fallback_global_ratio"
        n_genes_used = max(n_genes_used, int(expressed.sum()))
    else:
        lo, hi = np.quantile(r_valid, [trim_fraction, 1.0 - trim_fraction])
        trim_mask = (r_valid >= lo) & (r_valid <= hi)
        n_trimmed = int(trim_mask.sum())
        if n_trimmed >= min_genes_for_prior:
            r_trim = r_valid[trim_mask]
            method_detail = "empirical_bayes_trimmed_median_mad"
        else:
            # Quantile trim would leave too few genes; use all valid ratios instead.
            r_trim = r_valid
            method_detail = "empirical_bayes_median_mad"
        prior_mean_log = float(np.median(r_trim))
        tau = _robust_mad_scale(r_trim)
        n_genes_used = int(r_trim.size)

    # prior_weight scales observation-precision pseudocount (not direct shrinkage like
    # heuristic_shrink). Use a low floor so the default tunable range (0.5–5.0) is not
    # flattened by the legacy count_pseudocount=1.0 ceiling.
    _count_pseudo_floor = 0.05
    count_pseudo = max(_count_pseudo_floor, prior_weight * 0.2) * float(count_pseudocount)

    eb_prior = {
        "prior_mean_log": prior_mean_log,
        "tau_squared": float(tau**2),
        "count_pseudocount": float(count_pseudo),
        "eps": float(eps),
        "trim_fraction": float(trim_fraction),
    }
    meta: dict[str, Any] = {
        "gamma_prior_mean": prior_mean_log,
        "gamma_prior_scale": float(tau),
        "gamma_prior_tau_squared": float(tau**2),
        "n_genes_used_for_prior": n_genes_used,
        "gamma_method_detailed": method_detail,
        "fallback_triggered": fallback_triggered,
        "count_pseudocount": float(count_pseudo),
        "eb_prior": eb_prior,
    }
    return eb_prior, meta


def _apply_empirical_bayes_gamma(
    U_r: np.ndarray,
    S_r: np.ndarray,
    eb_prior: dict[str, Any],
    *,
    n_ref: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Shrink per-gene reference log-ratios toward a fixed prior.

    n_ref: number of cells in the reference group. Required for correct delta-method
           variance of the mean log-ratio: Var(log(Ū/S̄)) ≈ 1/(n_r * Ū) + 1/(n_r * S̄).
           Previously the n_r factor was missing, causing over-shrinkage.

    Returns (gamma_ref, shrinkage_weight, posterior_log_sd).
    """
    eps = float(eb_prior.get("eps", _EPS))
    prior_mean = float(eb_prior["prior_mean_log"])
    tau2 = float(eb_prior["tau_squared"])
    c = float(eb_prior.get("count_pseudocount", 1.0))

    U_r = np.asarray(U_r, dtype=float).ravel()
    S_r = np.asarray(S_r, dtype=float).ravel()

    n_ref = max(1.0, float(n_ref))
    r_g = np.log((U_r + eps) / (S_r + eps))
    sigma2 = 1.0 / (n_ref * U_r + c) + 1.0 / (n_ref * S_r + c)
    sigma2 = np.maximum(sigma2, 1e-12)

    w_g = tau2 / (tau2 + sigma2)
    w_g = np.clip(w_g, 0.0, 1.0)
    r_post = w_g * r_g + (1.0 - w_g) * prior_mean
    gamma_ref = np.exp(r_post)

    post_var = tau2 * sigma2 / (tau2 + sigma2)
    posterior_log_sd = np.sqrt(np.maximum(post_var, 0.0))

    return (
        np.nan_to_num(gamma_ref),
        np.nan_to_num(w_g),
        np.nan_to_num(posterior_log_sd),
    )


def _shrinkage_summary(weights: np.ndarray) -> dict[str, float]:
    w = np.asarray(weights, dtype=float)
    w = w[np.isfinite(w)]
    if w.size == 0:
        return {
            "mean": np.nan,
            "q10": np.nan,
            "q25": np.nan,
            "q50": np.nan,
            "q75": np.nan,
            "q90": np.nan,
        }
    qs = np.quantile(w, [0.1, 0.25, 0.5, 0.75, 0.9])
    return {
        "mean": float(np.mean(w)),
        "q10": float(qs[0]),
        "q25": float(qs[1]),
        "q50": float(qs[2]),
        "q75": float(qs[3]),
        "q90": float(qs[4]),
    }


def _gamma_stats(
    gamma_ref: np.ndarray, posterior_log_sd: np.ndarray | None = None
) -> dict[str, Any]:
    g = np.asarray(gamma_ref, dtype=float)
    g_fin = g[np.isfinite(g)]
    stats: dict[str, Any] = {
        "median": float(np.median(g_fin)) if g_fin.size else np.nan,
        "mean": float(np.mean(g_fin)) if g_fin.size else np.nan,
        "min": float(np.min(g_fin)) if g_fin.size else np.nan,
        "max": float(np.max(g_fin)) if g_fin.size else np.nan,
        "n_finite": int(np.isfinite(g).sum()),
    }
    if posterior_log_sd is not None:
        sd = np.asarray(posterior_log_sd, dtype=float)
        sd_fin = sd[np.isfinite(sd)]
        stats["posterior_log_sd_median"] = float(np.median(sd_fin)) if sd_fin.size else np.nan
        stats["posterior_log_sd_mean"] = float(np.mean(sd_fin)) if sd_fin.size else np.nan
    return stats


def _compute_velocity_delta(
    uns_layer: Any,
    spl_layer: Any,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float = 5.0,
    gamma_method: str = "heuristic_shrink",
    eb_prior: dict[str, Any] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Groupwise velocity delta (U_t - gamma_ref * S_t) on arbitrary expression layers.

    Used by the heuristic track (normalized unspliced/spliced) and, after scVelo
    moments smoothing, by the advanced track (Mu/Ms).

    gamma_method:
        - "heuristic_shrink": original global ratio + additive shrinkage using prior_weight
        - "robust_median": variant of heuristic_shrink that uses the *median of per-gene*
          U/S ratios (instead of global sum ratio) as the anchor for additive pseudo-count
          shrinkage. This is *not* a Bayesian/empirical-Bayes method; documented as a
          robust heuristic alternative for small reference groups.
        - "empirical_bayes": hierarchical gamma estimation via robust log-ratio
          empirical Bayes shrinkage (see module docstring). Recommended for small ref.
        - "raw" or other: minimal/no shrinkage (per-gene raw ratio)

    eb_prior:
        When ``gamma_method="empirical_bayes"``, pass a fixed prior dict (from the main
        analysis) to reuse the same hyperparameters during permutation shuffles.

    Returns
    -------
    delta_velocity, total_us, gamma_ref, gamma_info
        gamma_info carries diagnostics (empty for legacy methods except basic stats).
    """
    U_t = _get_group_mean(uns_layer, t_mask)
    S_t = _get_group_mean(spl_layer, t_mask)
    U_r = _get_group_mean(uns_layer, r_mask)
    S_r = _get_group_mean(spl_layer, r_mask)

    eps = _EPS
    gamma_info: dict[str, Any] = {"gamma_method": gamma_method}

    n_ref = float(np.sum(r_mask)) if r_mask is not None else 1.0

    if gamma_method == "empirical_bayes":
        if eb_prior is None:
            eb_prior, prior_meta = _estimate_eb_prior_from_reference(
                U_r, S_r, prior_weight=prior_weight
            )
            gamma_info.update(prior_meta)
        else:
            gamma_info.update(
                {
                    "used_fixed_prior": True,
                    "gamma_method_detailed": "empirical_bayes_fixed_prior",
                    "fallback_triggered": False,
                    "eb_prior": eb_prior,
                    "gamma_prior_mean": eb_prior.get("prior_mean_log"),
                    "gamma_prior_tau_squared": eb_prior.get("tau_squared"),
                }
            )
        gamma_ref, w_g, post_sd = _apply_empirical_bayes_gamma(U_r, S_r, eb_prior, n_ref=n_ref)
        gamma_info["shrinkage_summary"] = _shrinkage_summary(w_g)
        gamma_info["shrinkage_weights"] = w_g
        gamma_info["posterior_log_sd"] = post_sd
        gamma_info["effective_gamma_stats"] = _gamma_stats(gamma_ref, post_sd)
        gamma_info["eb_prior"] = eb_prior
    elif gamma_method == "robust_median":
        # Exclude zero-expression genes from the median anchor: (eps/eps)≈1 would
        # otherwise dominate on sparse scRNA-seq and pull base_gamma toward 1.
        raw_ratios = (U_r + eps) / (S_r + eps)
        expressed = (U_r + S_r) > 0
        ratios_use = raw_ratios[expressed] if np.any(expressed) else raw_ratios
        base_gamma = (
            float(np.median(ratios_use))
            if ratios_use.size > 0
            else (float(np.sum(U_r)) + eps) / (float(np.sum(S_r)) + eps)
        )
        beta = prior_weight
        alpha = base_gamma * beta
        gamma_ref = (U_r + alpha) / (S_r + beta)
        gamma_info["effective_gamma_stats"] = _gamma_stats(gamma_ref)
        gamma_info["gamma_method_detailed"] = "robust_median_additive_shrink"
        gamma_info["n_genes_used_for_median_anchor"] = int(ratios_use.size)
    elif gamma_method in ("heuristic_shrink", None, ""):
        global_gamma = (np.sum(U_r) + eps) / (np.sum(S_r) + eps)
        beta = prior_weight
        alpha = global_gamma * beta
        gamma_ref = (U_r + alpha) / (S_r + beta)
        gamma_info["effective_gamma_stats"] = _gamma_stats(gamma_ref)
        gamma_info["gamma_method_detailed"] = "heuristic_shrink_additive"
    else:
        # raw: per-gene U/S; zero-expression ref genes get the global sum ratio
        # (not eps/eps≈1, which mis-scales excess when S_t > 0).
        expressed = (U_r + S_r) > 0
        gamma_ref = np.empty_like(U_r, dtype=float)
        gamma_ref[expressed] = (U_r[expressed] + eps) / (S_r[expressed] + eps)
        global_gamma = (float(np.sum(U_r)) + eps) / (float(np.sum(S_r)) + eps)
        gamma_ref[~expressed] = global_gamma
        gamma_info["effective_gamma_stats"] = _gamma_stats(gamma_ref)
        gamma_info["gamma_method_detailed"] = "raw_ratio"
        gamma_info["n_genes_raw_ratio"] = int(np.sum(expressed))
        gamma_info["n_genes_global_ratio_fallback"] = int(np.sum(~expressed))
        gamma_info["global_ratio_fallback"] = float(global_gamma)

    delta_velocity = U_t - (gamma_ref * S_t)

    total_uns = np.asarray(uns_layer.sum(axis=0)).ravel()
    total_spl = np.asarray(spl_layer.sum(axis=0)).ravel()
    total_us = total_uns + total_spl
    return (
        np.nan_to_num(delta_velocity),
        np.nan_to_num(total_us),
        np.nan_to_num(gamma_ref),
        gamma_info,
    )


def _compute_moments_velocity_delta(
    adata_comp: ad.AnnData,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float = 5.0,
    gamma_method: str = "heuristic_shrink",
    eb_prior: dict[str, Any] | None = None,
    n_neighbors: int = 30,
    n_pcs: int = 30,
    use_precomputed: bool = False,
    recompute_neighbors: bool = True,
    random_state: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Advanced track: use scVelo moments (Mu/Ms) for local smoothing before delta.

    Supports gamma_method for reference gamma estimation robustness.
    """
    try:
        import scvelo as scv
    except ImportError as e:
        raise ImportError(
            "mode='advanced' requires the 'scvelo' package. "
            "Please install it with: pip install scvelo or 'scatrans[advanced]'"
        ) from e

    info: dict[str, Any] = {}
    had_moments_before = "Mu" in adata_comp.layers and "Ms" in adata_comp.layers

    if use_precomputed and "Mu" in adata_comp.layers and "Ms" in adata_comp.layers:
        logger.info("Using precomputed Mu/Ms layers.")
        used_precomputed = True
        neighbors_source = "precomputed"
        n_neighbors_eff = None
        n_pcs_eff = None

        if _matrix_shape(adata_comp.layers["Mu"]) != tuple(adata_comp.shape):
            raise ValueError("Layer 'Mu' shape does not match adata shape.")
        if _matrix_shape(adata_comp.layers["Ms"]) != tuple(adata_comp.shape):
            raise ValueError("Layer 'Ms' shape does not match adata shape.")
    else:
        used_precomputed = False
        if "spliced" not in adata_comp.layers or "unspliced" not in adata_comp.layers:
            raise ValueError(
                "Advanced mode requires 'spliced' and 'unspliced' layers "
                "(after any automatic kb_python 'mature'/'nascent' remapping)."
            )

        n_obs = adata_comp.n_obs
        if n_obs < 5:
            raise ValueError(f"Advanced mode requires at least 5 observations (got {n_obs}).")

        if n_obs < 30:
            logger.warning(
                "Advanced mode has only %d observations; scVelo moments may be unstable.", n_obs
            )

        n_neighbors_eff = min(n_neighbors, max(2, n_obs - 1))
        n_pcs_eff = min(n_pcs, max(1, n_obs - 1), max(1, adata_comp.n_vars - 1))

        if n_pcs_eff < 2:
            raise ValueError("Too few observations/genes for reliable PCA in advanced mode.")

        if recompute_neighbors or "neighbors" not in adata_comp.uns:
            sc.pp.neighbors(
                adata_comp,
                n_neighbors=n_neighbors_eff,
                n_pcs=n_pcs_eff,
                random_state=random_state,
            )
            neighbors_source = "computed_by_scatrans"
        else:
            neighbors_source = "preexisting_neighbors"

        scv.pp.moments(adata_comp, n_pcs=n_pcs_eff, n_neighbors=n_neighbors_eff)

        if _matrix_shape(adata_comp.layers["Mu"]) != tuple(adata_comp.shape):
            raise ValueError("Layer 'Mu' shape does not match adata shape after moments.")
        if _matrix_shape(adata_comp.layers["Ms"]) != tuple(adata_comp.shape):
            raise ValueError("Layer 'Ms' shape does not match adata shape after moments.")

        info.update(
            {
                "n_neighbors_effective": n_neighbors_eff,
                "n_pcs_effective": n_pcs_eff,
            }
        )

    delta_velocity, total_us, gamma_ref, gamma_info = _compute_velocity_delta(
        adata_comp.layers["Mu"],
        adata_comp.layers["Ms"],
        t_mask,
        r_mask,
        prior_weight,
        gamma_method=gamma_method,
        eb_prior=eb_prior,
    )

    info.update(
        {
            "used_precomputed_moments": used_precomputed,
            "neighbors_source": neighbors_source,
            "had_moments_before": had_moments_before,
            "n_neighbors_requested": n_neighbors,
            "n_pcs_requested": n_pcs,
            "recompute_neighbors": recompute_neighbors,
            "random_state": random_state,
            "Mu_sparse": sparse.issparse(adata_comp.layers["Mu"]),
            "Ms_sparse": sparse.issparse(adata_comp.layers["Ms"]),
            "moments_shape": tuple(adata_comp.shape),
            "gamma_info": gamma_info,
        }
    )

    return delta_velocity, total_us, gamma_ref, info
