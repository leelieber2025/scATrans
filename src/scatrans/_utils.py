"""
scATrans internal utilities (not part of public API).

Small, pure or near-pure helper functions extracted from the original tl.py
to keep the core active_score readable and to enable reuse (esp. bias correction).
"""

from __future__ import annotations

import logging
from math import comb  # re-exported for permutation use
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.linear_model import HuberRegressor

logger = logging.getLogger(__name__)


# Re-export for modules that need it without importing math directly
__all__ = [
    "comb",
    "_is_integer_counts_like",
    "_warn_if_not_integer_counts_matrix",
    "_warn_if_low_counts_matrix",
    "_safe_add_matrices",
    "_normalize_velocity_layers_by_size_factor",
    "_get_group_mean",
    "_get_exponential_scale_lambda",
    "_soft_scale",
    "_pseudobulk_with_layers",
    "_fit_huber_bias_correction",
]


def _is_integer_counts_like(X: Any, max_check: int = 100000) -> bool:
    if sparse.issparse(X):
        data = X.data
        if data.size == 0:
            return True
        if not np.all(np.isfinite(data)):
            return False
        vals = data
    else:
        arr = np.asarray(X)
        if not np.all(np.isfinite(arr)):
            return False
        vals = arr.ravel()

    if vals.size == 0:
        return True

    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)

    return np.all(vals >= 0) and np.allclose(vals, np.round(vals))


def _warn_if_not_integer_counts_matrix(X: Any, max_check: int = 100000) -> None:
    if not _is_integer_counts_like(X, max_check=max_check):
        logger.warning(
            "Data passed to PyDESeq2 may not be raw non-negative integer counts. "
            "Please ensure the input contains unnormalized counts."
        )


def _warn_if_low_counts_matrix(X: Any, max_check: int = 100000) -> None:
    vals = X.data if sparse.issparse(X) else np.asarray(X).ravel()

    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return

    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)

    if vals.max() < 30:
        logger.warning(
            "Maximum count passed to PyDESeq2 is <30. This may be valid for small datasets, "
            "but please verify that the matrix contains raw counts, not normalized/log-transformed values."
        )


def _safe_add_matrices(a: Any, b: Any) -> Any:
    if sparse.issparse(a) or sparse.issparse(b):
        return sparse.csr_matrix(a) + sparse.csr_matrix(b)
    return np.asarray(a) + np.asarray(b)


def _normalize_velocity_layers_by_size_factor(
    uns_layer: Any, spl_layer: Any, target_sum: float | None = None
) -> tuple[Any, Any, np.ndarray, np.ndarray]:
    total_layer = _safe_add_matrices(uns_layer, spl_layer)
    row_totals = np.asarray(total_layer.sum(axis=1)).ravel()
    positive = row_totals > 0

    if positive.sum() == 0:
        return uns_layer, spl_layer, row_totals, np.ones_like(row_totals, dtype=float)

    if target_sum is None:
        target_sum = np.median(row_totals[positive])

    factors = target_sum / np.maximum(row_totals, 1e-8)

    if sparse.issparse(uns_layer) or sparse.issparse(spl_layer):
        uns_layer = sparse.csr_matrix(uns_layer)
        spl_layer = sparse.csr_matrix(spl_layer)
        D = sparse.diags(factors)
        return D @ uns_layer, D @ spl_layer, row_totals, factors

    return (
        np.asarray(uns_layer) * factors[:, None],
        np.asarray(spl_layer) * factors[:, None],
        row_totals,
        factors,
    )


def _get_group_mean(matrix: Any, mask: np.ndarray) -> np.ndarray:
    if np.sum(mask) == 0:
        raise ValueError("Cannot compute group mean for an empty group.")
    sub = matrix[mask]
    if sparse.issparse(sub):
        return np.asarray(sub.mean(axis=0)).ravel()
    return np.asarray(sub.mean(axis=0)).ravel()


def _get_exponential_scale_lambda(x: np.ndarray) -> float:
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x_pos = np.clip(x, 0.0, None)
    nonzero_x = x_pos[x_pos > 0]
    if len(nonzero_x) < 2:
        return 1e-8
    med = np.median(nonzero_x)
    return med / np.log(2.0) if med > 0 else 1e-8


def _soft_scale(x: np.ndarray, lam: float) -> np.ndarray:
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x_pos = np.clip(x, 0.0, None)
    if lam <= 1e-8:
        return np.zeros_like(x)
    return 1.0 - np.exp(-x_pos / lam)


def _pseudobulk_with_layers(
    adata: ad.AnnData,
    sample_col: str,
    groupby: str,
    layers: Iterable[str] = ("spliced", "unspliced"),
    x_layer: str | None = None,
    use_total_for_x: bool = False,
    min_cells: int = 10,
    min_counts: int = 1000,
) -> ad.AnnData:
    """Aggregate to pseudobulk while preserving the requested layers."""
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col='{sample_col}' not found.")
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby='{groupby}' not found.")
    for layer in layers:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers")

    if use_total_for_x:
        X_source = _safe_add_matrices(adata.layers["spliced"], adata.layers["unspliced"])
        x_source_name = "spliced + unspliced"
    else:
        if x_layer is not None and x_layer not in adata.layers:
            raise ValueError(f"x_layer '{x_layer}' not found in adata.layers")
        X_source = adata.X if x_layer is None else adata.layers[x_layer]
        x_source_name = "adata.X" if x_layer is None else f"layer '{x_layer}'"

    group_df = adata.obs[[sample_col, groupby]].copy()
    group_df[sample_col] = group_df[sample_col].astype(str)
    group_df[groupby] = group_df[groupby].astype(str)
    pb_key = group_df[sample_col] + "||" + group_df[groupby]
    unique_keys = pd.Index(pb_key.unique())

    X_rows, obs_rows = [], []
    layer_rows: dict[str, list] = {layer: [] for layer in layers}

    for key in unique_keys:
        mask = pb_key.values == key
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        x_sum = np.nan_to_num(np.asarray(X_source[mask].sum(axis=0)).ravel())
        if float(x_sum.sum()) < min_counts:
            continue

        sample_id, group_value = key.split("||", 1)
        X_rows.append(sparse.csr_matrix(x_sum.reshape(1, -1)))
        obs_rows.append(
            {
                sample_col: sample_id,
                groupby: group_value,
                "n_cells": n_cells,
                "total_counts": float(x_sum.sum()),
                "pb_x_source": x_source_name,
            }
        )
        for layer in layers:
            l_sum = np.nan_to_num(np.asarray(adata.layers[layer][mask].sum(axis=0)).ravel())
            layer_rows[layer].append(sparse.csr_matrix(l_sum.reshape(1, -1)))

    if not X_rows:
        raise ValueError("No samples remained after pseudobulk filtering.")

    adata_pb = ad.AnnData(
        X=sparse.vstack(X_rows).tocsr(),
        obs=pd.DataFrame(obs_rows),
        var=adata.var.copy(),
    )
    adata_pb.obs.index = (
        adata_pb.obs[sample_col].astype(str) + "_" + adata_pb.obs[groupby].astype(str)
    )
    for layer in layers:
        adata_pb.layers[layer] = sparse.vstack(layer_rows[layer]).tocsr()
    adata_pb.obs_names_make_unique()
    return adata_pb


def _is_bias_correction_enabled(val: Any) -> bool:
    """Return True unless the user explicitly disabled bias correction."""
    if val is None or val is False:
        return False
    if isinstance(val, str):
        v = val.strip().lower()
        if v in ("none", "off", "no", "false", "disable", ""):
            return False
    return True


def _fit_huber_bias_correction(
    delta_velocity: np.ndarray,
    gene_length: np.ndarray,
    intron_number: np.ndarray,
    total_us_for_weights: np.ndarray,
    valid_feat: np.ndarray,
    valid_expr: np.ndarray,
    X_features: np.ndarray | None,
    min_fit_obs: int = 30,
    huber_epsilon: float = 1.35,
    huber_max_iter: int = 500,
    bias_correction: str = "huber_length_intron",
) -> tuple[np.ndarray, dict[str, Any]]:
    """
    Shared Huber regression bias correction (or median fallback).

    Used by both the main analysis path and the permutation tasks so the
    correction logic stays in one place (DRY).

    bias_correction controls behavior:
      - "huber_length_intron" (default), "huber", "yes", "on": perform the
        length+intron Huber correction (with median fallback if regression
        cannot be fit).
      - "none", "off", False, None: disable correction entirely; residual is
        the raw delta_velocity (no subtraction of fit or median). This keeps
        the basic analysis clean for users who do not want the correction.

    Returns (residual, bias_info_dict) where bias_info contains:
      - "bias_corrected": bool
      - "method": the effective method used ("huber_length_intron" or "none")
      - "n_genes_used_for_fit": int
      - "fallback_to_median": bool
      - "coef_gene_length", "coef_intron_number" (if regression succeeded)
      - "intercept" (if available)
    """
    residual = np.zeros_like(delta_velocity, dtype=float)
    method = str(bias_correction) if bias_correction is not None else "huber_length_intron"
    bias_info: dict[str, Any] = {
        "bias_corrected": False,
        "method": method,
        "n_genes_used_for_fit": 0,
        "fallback_to_median": False,
        "coef_gene_length": np.nan,
        "coef_intron_number": np.nan,
        "intercept": np.nan,
    }

    if not _is_bias_correction_enabled(bias_correction):
        # No correction at all: residual == raw delta (clipped for invalid expr)
        residual = np.array(delta_velocity, dtype=float, copy=True)
        residual[~valid_expr] = 0.0
        bias_info["bias_corrected"] = False
        bias_info["fallback_to_median"] = False
        bias_info["n_genes_used_for_fit"] = 0
        return residual, bias_info

    fit_mask = valid_feat & valid_expr
    n_fit = int(fit_mask.sum())
    bias_info["n_genes_used_for_fit"] = n_fit

    regression_succeeded = False
    if X_features is not None and n_fit >= min_fit_obs:
        try:
            X_fit = np.column_stack(
                [
                    np.log1p(gene_length[fit_mask]),
                    np.log1p(intron_number[fit_mask]),
                ]
            )
            weights = np.clip(
                total_us_for_weights[fit_mask],
                a_min=None,
                a_max=np.percentile(total_us_for_weights[fit_mask], 95),
            )
            with warnings.catch_warnings():
                import warnings as _w

                _w.simplefilter("ignore")
                model = HuberRegressor(epsilon=huber_epsilon, max_iter=huber_max_iter).fit(
                    X_fit, delta_velocity[fit_mask], sample_weight=weights
                )
            pred = model.predict(X_features)
            residual[valid_feat] = delta_velocity[valid_feat] - pred
            regression_succeeded = True
            bias_info["bias_corrected"] = True
            if hasattr(model, "coef_") and len(model.coef_) >= 2:
                bias_info["coef_gene_length"] = float(model.coef_[0])
                bias_info["coef_intron_number"] = float(model.coef_[1])
            if hasattr(model, "intercept_"):
                bias_info["intercept"] = float(model.intercept_)
        except Exception as e:
            logger.warning("Bias correction failed. Falling back to median. Reason: %s", e)

    if not regression_succeeded and valid_expr.sum() > 0:
        residual[valid_expr] = delta_velocity[valid_expr] - np.nanmedian(delta_velocity[valid_expr])
        bias_info["fallback_to_median"] = True
        bias_info["bias_corrected"] = True  # median correction still applied

    residual[~valid_expr] = 0.0
    return residual, bias_info


# warnings is used inside _fit_huber_bias_correction
import warnings  # noqa: E402  (executed at import time)
