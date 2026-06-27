"""
scATrans internal utilities (not part of public API).

Small, pure or near-pure helper functions extracted from the original tl.py
to keep the core active_score readable and to enable reuse (esp. bias correction).
"""

from __future__ import annotations

import logging
import uuid
import warnings
from math import comb  # re-exported for permutation use
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import HuberRegressor

logger = logging.getLogger(__name__)


# Re-export for modules that need it without importing math directly
# Primary result column names (public API); legacy velocity_* aliases are kept in sync.
UNSPLICED_EXCESS_DELTA_COL = "unspliced_excess_delta"
UNSPLICED_EXCESS_RESIDUAL_COL = "unspliced_excess_residual"
UNSPLICED_EXCESS_PVAL_COL = "unspliced_excess_pval"
UNSPLICED_EXCESS_FDR_COL = "unspliced_excess_fdr"
LEGACY_VELOCITY_DELTA_COL = "velocity_delta_raw"
LEGACY_VELOCITY_RESIDUAL_COL = "velocity_residual"


def _resolve_results_column(
    df: pd.DataFrame, primary: str, legacy: str, *, required: bool = True
) -> str:
    """Return *primary* if present, else *legacy*; raise if neither and required."""
    if primary in df.columns:
        return primary
    if legacy in df.columns:
        return legacy
    if required:
        raise KeyError(f"Expected column '{primary}' (or legacy '{legacy}') in results DataFrame.")
    return primary


def _write_unspliced_excess_columns(
    var_df: pd.DataFrame,
    *,
    delta: np.ndarray,
    residual: np.ndarray,
) -> None:
    """Write primary unspliced-excess columns and deprecated velocity aliases."""
    var_df[UNSPLICED_EXCESS_DELTA_COL] = delta
    var_df[UNSPLICED_EXCESS_RESIDUAL_COL] = residual
    var_df[LEGACY_VELOCITY_DELTA_COL] = delta
    var_df[LEGACY_VELOCITY_RESIDUAL_COL] = residual


__all__ = [
    "comb",
    "UNSPLICED_EXCESS_DELTA_COL",
    "UNSPLICED_EXCESS_RESIDUAL_COL",
    "UNSPLICED_EXCESS_PVAL_COL",
    "UNSPLICED_EXCESS_FDR_COL",
    "LEGACY_VELOCITY_DELTA_COL",
    "LEGACY_VELOCITY_RESIDUAL_COL",
    "_resolve_results_column",
    "_write_unspliced_excess_columns",
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
    "_resolve_aligned_raw_counts",
    "_prepare_log_normalized_expression",
]


def _is_integer_counts_like(X: Any, max_check: int = 100000, atol: float = 1e-6) -> bool:
    """Return True if the matrix contains non-negative values that are integer-valued
    (within tolerance). This is tolerant of float64 summed counts that are exactly
    (or very nearly) integers, which commonly occurs after pseudobulk aggregation.
    """
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

    # Tolerant check: allows tiny floating point noise from summation / cast
    rounded = np.round(vals)
    return np.all(vals >= 0) and np.allclose(vals, rounded, atol=atol, rtol=1e-5)


def _warn_if_not_integer_counts_matrix(X: Any, max_check: int = 100000) -> None:
    if not _is_integer_counts_like(X, max_check=max_check):
        logger.warning(
            "Data passed to PyDESeq2 may not be raw non-negative integer counts. "
            "Please ensure the input contains unnormalized counts."
        )


def _resolve_aligned_raw_counts(
    adata: ad.AnnData,
    *,
    layer: str = "counts",
    require_integer: bool = True,
) -> Any | None:
    """Return a count matrix aligned to ``adata.n_vars``, or None if unsafe to use.

    Refuses matrices whose second dimension does not match the current gene count.
    When ``raw_gene_list`` in ``.uns`` differs in length (typical after HVG subsetting),
    a counts layer that matches current ``var_names`` is still accepted for DE backends.
    """
    candidates: list[tuple[str, Any]] = []
    if layer in adata.layers:
        candidates.append((f"layers['{layer}']", adata.layers[layer]))
    raw = getattr(adata, "raw", None)
    if (
        raw is not None
        and raw.shape[1] == adata.n_vars
        and hasattr(raw, "var_names")
        and np.array_equal(raw.var_names, adata.var_names)
    ):
        candidates.append(("adata.raw", raw.X))

    for source_name, mat in candidates:
        n_cols = mat.shape[1] if hasattr(mat, "shape") else 0
        if n_cols != adata.n_vars:
            logger.warning(
                "Counts from %s have %d columns but adata has %d genes; skipping for count-based DE.",
                source_name,
                n_cols,
                adata.n_vars,
            )
            continue
        if require_integer and not _is_integer_counts_like(mat):
            logger.warning(
                "Counts from %s do not look like raw integer counts; skipping for count-based DE.",
                source_name,
            )
            continue

        raw_gene_list = adata.uns.get("scatrans", {}).get("raw_gene_list")
        if raw_gene_list is not None:
            stored = np.asarray(raw_gene_list)
            current = adata.var_names.to_numpy()
            if len(stored) == adata.n_vars and not np.array_equal(stored, current):
                logger.warning(
                    "Counts from %s match n_vars but stored raw_gene_list order differs from "
                    "adata.var_names. Refusing misaligned counts for count-based DE. "
                    "Re-run store_raw_counts() on the current object.",
                    source_name,
                )
                continue
            if len(stored) != adata.n_vars:
                logger.info(
                    "raw_gene_list (%d genes) differs from current n_vars (%d). "
                    "Using %s aligned to current genes for DE; enrichment universe still uses "
                    "the preserved full raw_gene_list.",
                    len(stored),
                    adata.n_vars,
                    source_name,
                )
        return mat

    return None


def _prepare_log_normalized_expression(ad_expr: ad.AnnData) -> np.ndarray:
    """Dense log1p library-size normalized matrix for mixed models (no double log1p)."""
    ad_work = ad_expr.copy()
    if "log1p" in ad_work.uns:
        X = ad_work.X.toarray() if sparse.issparse(ad_work.X) else np.asarray(ad_work.X)
        return np.asarray(X, dtype=float)

    if _is_integer_counts_like(ad_work.X):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.normalize_total(ad_work, target_sum=1e4)
            sc.pp.log1p(ad_work)
        X = ad_work.X.toarray() if sparse.issparse(ad_work.X) else np.asarray(ad_work.X)
        return np.asarray(X, dtype=float)

    X = ad_work.X.toarray() if sparse.issparse(ad_work.X) else np.asarray(ad_work.X)
    X = np.asarray(X, dtype=float)
    finite = X[np.isfinite(X)]
    if finite.size:
        mx = np.nanmax(finite)
        has_neg = np.any(finite < 0)
        # Heuristic for already log-transformed data (e.g. log1p, scran, SCT residuals):
        # - negatives or max <=20 : treat as already transformed (closes the 5<mx<=20 gap)
        # - mx >20 on non-negative: likely raw-ish counts, apply log1p + warn
        if has_neg or mx <= 20:
            return X
        if mx > 20:
            logger.warning(
                "Mixed model input is neither marked log1p nor integer counts; applying log1p "
                "for stability (max value=%.1f). If this is already log-scale, set de_preprocess='none' "
                "or pre-log the input.",
                mx,
            )
            return np.log1p(np.clip(X, 0, None))
    return X


def _warn_if_low_counts_matrix(X: Any, max_check: int = 100000) -> None:
    vals = X.data if sparse.issparse(X) else np.asarray(X).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return

    # Use rounded view for the max check (handles float64 summed pseudobulk)
    rounded_max = np.max(np.round(vals))

    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)

    if rounded_max < 30:
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
    layers: Iterable[str] = (),
    x_layer: str | None = None,
    use_total_for_x: bool = False,
    min_cells: int = 10,
    min_counts: int = 1000,
) -> ad.AnnData:
    """Aggregate to pseudobulk while preserving the requested layers.

    layers: which .layers to aggregate and carry through (e.g. velocity layers).
            Default empty so pure-DE callers do not accidentally require spliced/unspliced.
    use_total_for_x=True requires spliced+unspliced to exist (independent of layers list).
    """
    if sample_col not in adata.obs.columns:
        raise ValueError(f"sample_col='{sample_col}' not found.")
    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby='{groupby}' not found.")
    for layer in layers:
        if layer not in adata.layers:
            raise ValueError(f"Layer '{layer}' not found in adata.layers")

    if use_total_for_x:
        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            raise ValueError(
                "use_total_for_x=True requires both 'spliced' and 'unspliced' layers to be present."
            )
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

    # Use a per-run UUID separator (printable, no embedded NUL) that is vanishingly unlikely
    # to appear in real sample/group names. Avoids pandas str-concat truncation with \0 and
    # eliminates "||" injection / mis-split risk from earlier fragile separator.
    _pb_sep = f"__scAT_PB_{uuid.uuid4().hex}__"
    pb_key = group_df[sample_col] + _pb_sep + group_df[groupby]
    unique_keys = pd.Index(pb_key.unique())

    X_rows, obs_rows = [], []
    layer_rows: dict[str, list] = {layer: [] for layer in layers}

    for key in unique_keys:
        mask = pb_key.values == key
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        x_sum = np.nan_to_num(np.asarray(X_source[mask].sum(axis=0)).ravel())
        # Clean to integer-valued floats for count-like data (pseudobulk sums).
        # Using round keeps the numeric value exact while float dtype is fine for AnnData/sparse.
        x_sum = np.round(x_sum).astype(np.float64, copy=False)

        if float(x_sum.sum()) < min_counts:
            continue

        sample_id, group_value = key.split(_pb_sep, 1)
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
            # Velocity layers (spliced/unspliced) can stay as summed float; only round for cleanliness
            l_sum = np.round(l_sum).astype(np.float64, copy=False)
            layer_rows[layer].append(sparse.csr_matrix(l_sum.reshape(1, -1)))

    if not X_rows:
        raise ValueError("No samples remained after pseudobulk filtering.")

    with warnings.catch_warnings():
        warnings.simplefilter(
            "ignore", category=UserWarning
        )  # pandas "Transforming to str index" during AnnData/obs construction is benign
        adata_pb = ad.AnnData(
            X=sparse.vstack(X_rows).tocsr(),
            obs=pd.DataFrame(obs_rows),
            var=adata.var.copy(),
        )
        # Use AnnData's obs_names setter (preferred)
        adata_pb.obs_names = (
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
                warnings.simplefilter("ignore")
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
            bias_info["fallback_reason"] = (
                f"huber_regression_failed: {type(e).__name__}: {str(e)[:200]}"
            )

    if not regression_succeeded and valid_expr.sum() > 0:
        residual[valid_expr] = delta_velocity[valid_expr] - np.nanmedian(delta_velocity[valid_expr])
        bias_info["fallback_to_median"] = True
        bias_info["bias_corrected"] = True  # median correction still applied

    residual[~valid_expr] = 0.0
    return residual, bias_info


# warnings imported at top of file (used inside _fit_huber_bias_correction and pseudobulk creation)
