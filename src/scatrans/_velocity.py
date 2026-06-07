"""
scATrans internal velocity delta computation (heuristic + advanced scVelo moments track).

Extracted from original tl.py.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Tuple

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse

from ._utils import _get_group_mean

logger = logging.getLogger(__name__)


def _compute_velocity_delta(
    uns_layer: Any,
    spl_layer: Any,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float = 5.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """Classic (global ratio) velocity delta used by the heuristic track."""
    U_t = _get_group_mean(uns_layer, t_mask)
    S_t = _get_group_mean(spl_layer, t_mask)
    U_r = _get_group_mean(uns_layer, r_mask)
    S_r = _get_group_mean(spl_layer, r_mask)

    global_gamma = (np.sum(U_r) + 1e-8) / (np.sum(S_r) + 1e-8)
    beta = prior_weight
    alpha = global_gamma * beta
    gamma_ref = (U_r + alpha) / (S_r + beta)
    delta_velocity = U_t - (gamma_ref * S_t)

    total_uns = np.asarray(uns_layer.sum(axis=0)).ravel()
    total_spl = np.asarray(spl_layer.sum(axis=0)).ravel()
    total_us = total_uns + total_spl
    return np.nan_to_num(delta_velocity), np.nan_to_num(total_us)


def _compute_moments_velocity_delta(
    adata_comp: ad.AnnData,
    t_mask: np.ndarray,
    r_mask: np.ndarray,
    prior_weight: float = 5.0,
    n_neighbors: int = 30,
    n_pcs: int = 30,
    use_precomputed: bool = False,
    recompute_neighbors: bool = True,
    random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Advanced track: use scVelo moments (Mu/Ms) for local smoothing before delta."""
    try:
        import scvelo as scv
    except ImportError as e:
        raise ImportError(
            "mode='advanced' requires the 'scvelo' package. "
            "Please install it with: pip install scvelo or 'scatrans[advanced]'"
        ) from e

    info: Dict[str, Any] = {}
    had_moments_before = "Mu" in adata_comp.layers and "Ms" in adata_comp.layers

    if use_precomputed and "Mu" in adata_comp.layers and "Ms" in adata_comp.layers:
        logger.info("Using precomputed Mu/Ms layers.")
        used_precomputed = True
        neighbors_source = "precomputed"
        n_neighbors_eff = None
        n_pcs_eff = None

        if adata_comp.layers["Mu"].shape != adata_comp.shape:
            raise ValueError("Layer 'Mu' shape does not match adata shape.")
        if adata_comp.layers["Ms"].shape != adata_comp.shape:
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

        if adata_comp.layers["Mu"].shape != adata_comp.shape:
            raise ValueError("Layer 'Mu' shape does not match adata shape after moments.")
        if adata_comp.layers["Ms"].shape != adata_comp.shape:
            raise ValueError("Layer 'Ms' shape does not match adata shape after moments.")

        info.update(
            {
                "n_neighbors_effective": n_neighbors_eff,
                "n_pcs_effective": n_pcs_eff,
            }
        )

    delta_velocity, total_us = _compute_velocity_delta(
        adata_comp.layers["Mu"],
        adata_comp.layers["Ms"],
        t_mask,
        r_mask,
        prior_weight,
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
        }
    )

    return delta_velocity, total_us, info
