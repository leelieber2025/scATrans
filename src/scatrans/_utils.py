"""
scATrans internal utilities (not part of public API).

Small, pure or near-pure helper functions extracted from the original tl.py
to keep the core active_score readable and to enable reuse (esp. bias correction).
"""

from __future__ import annotations

import logging
import uuid
import warnings
from collections.abc import Iterable
from math import comb  # re-exported for permutation use
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.linear_model import HuberRegressor

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# ---------------------------------------------------------------------------
# Matrix helpers
# ---------------------------------------------------------------------------
# AnnData types ``.X`` / ``.layers[...]`` as a very wide union (ndarray | sparse |
# dask | None | …). Calling ``.sum()`` / ``.toarray()`` / ``.copy()`` / ``.shape``
# directly triggers dozens of mypy [union-attr] / incomplete-type errors.
# Route all of that through these helpers so call sites stay readable.


def _require_matrix(X: Any, *, name: str = "matrix") -> Any:
    """Return *X* if non-None; raise :class:`ValueError` otherwise."""
    if X is None:
        raise ValueError(f"{name} is None")
    return X


def _dense_expression_matrix(X: Any) -> np.ndarray:
    """Dense ``float`` ndarray from AnnData ``.X``, a layer, sparse, or dense array."""
    X = _require_matrix(X, name="expression matrix")
    if sparse.issparse(X):
        return np.asarray(X.toarray(), dtype=float)
    # Some backends expose densify via toarray/todense without being scipy sparse.
    toarray = getattr(X, "toarray", None)
    if callable(toarray):
        try:
            return np.asarray(toarray(), dtype=float)
        except Exception:
            pass
    todense = getattr(X, "todense", None)
    if callable(todense):
        try:
            return np.asarray(todense(), dtype=float)
        except Exception:
            pass
    return np.asarray(X, dtype=float)


def _matrix_copy(X: Any) -> Any:
    """Copy a layer/X matrix without union-type attribute access at the call site."""
    X = _require_matrix(X)
    copy_fn = getattr(X, "copy", None)
    if callable(copy_fn):
        return copy_fn()
    return np.array(X, copy=True)


def _matrix_shape(X: Any) -> tuple[int, ...]:
    """Return ``X.shape`` as a plain int tuple (narrows incomplete layer types)."""
    X = _require_matrix(X)
    shape = getattr(X, "shape", None)
    if shape is None:
        raise TypeError(f"object of type {type(X)!r} has no shape")
    return tuple(int(s) for s in shape)


def _matrix_sum_axis0(X: Any) -> np.ndarray:
    """Column sums as a 1-d float array (works for sparse and dense)."""
    X = _require_matrix(X)
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=0)).ravel().astype(float, copy=False)
    arr = np.asarray(X, dtype=float)
    if arr.ndim == 1:
        return arr
    return arr.sum(axis=0).ravel()


def _matrix_sum_axis1(X: Any) -> np.ndarray:
    """Row sums as a 1-d float array (works for sparse and dense)."""
    X = _require_matrix(X)
    if sparse.issparse(X):
        return np.asarray(X.sum(axis=1)).ravel().astype(float, copy=False)
    arr = np.asarray(X, dtype=float)
    if arr.ndim == 1:
        return arr
    return arr.sum(axis=1).ravel()


def _matrix_row_subset_sum_axis0(X: Any, row_mask: Any) -> np.ndarray:
    """Column sums over a boolean row mask (pseudobulk aggregation helper)."""
    X = _require_matrix(X)
    sub = X[row_mask]
    return _matrix_sum_axis0(sub)


def _as_var_dataframe(adata: ad.AnnData) -> pd.DataFrame:
    """Narrow ``adata.var`` to :class:`~pandas.DataFrame` for ``.loc`` / sorting.

    Some anndata type stubs type ``.var`` as a broad union (e.g. including
    Dataset2D) that does not expose ``.loc``; runtime is always a DataFrame.
    """
    return cast(pd.DataFrame, adata.var)


def _normalize_group_label(val: Any) -> str | None:
    """Normalize a single group label for stable string matching.

    Handles NaN, bool, int/float (``1.0`` → ``"1"``), and **string** forms from
    CSV/Excel (``"1.0"`` → ``"1"``, ``"2.0"`` → ``"2"``) so callers can pass
    ``target_group=1`` or ``"1"`` against obs values stored as ``"1.0"``.
    """
    if val is None:
        return None
    try:
        if pd.isna(val):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(val, (bool, np.bool_)):
        return str(val)
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    if isinstance(val, (float, np.floating)):
        fv = float(val)
        if np.isnan(fv):
            return None
        if fv.is_integer():
            return str(int(fv))
        return str(fv)
    s = str(val).strip()
    if s.lower() in ("nan", "<na>", "none", ""):
        return None
    # Stringified numerics from CSV/Excel (e.g. "1.0", "2.00")
    try:
        fv = float(s)
        if np.isfinite(fv) and fv.is_integer():
            return str(int(fv))
    except (TypeError, ValueError):
        pass
    return s


def _normalize_label_array(labels: Any) -> np.ndarray:
    """Map labels through :func:`_normalize_group_label` (vectorized via pandas)."""
    ser = pd.Series(np.asarray(labels, dtype=object).ravel())
    return ser.map(_normalize_group_label).to_numpy()


def _as_contrast_categorical(
    labels: Any,
    reference_group: str,
    target_group: str,
) -> pd.Categorical:
    """Build a two-level Categorical with normalized labels (safe for float ``1.0`` vs ``"1"``)."""
    vals = _normalize_label_array(labels)
    # categories are expected already normalized by callers (public API validates)
    return pd.Categorical(vals, categories=[reference_group, target_group])


def _subset_obs_mask(
    obs_col: pd.Series,
    subset_values: Any,
) -> pd.Series:
    """Boolean mask for ``subset_col`` matching with group-label normalization.

    ``subset_values`` may be a single label or a sequence. Numeric / ``"1.0"``
    forms match the same way as ``groupby`` contrasts.
    """
    if isinstance(subset_values, (str, int, float, bool, np.integer, np.floating, np.bool_)):
        raw_list = [subset_values]
    else:
        try:
            raw_list = list(subset_values)
        except TypeError as exc:
            raise TypeError(
                f"subset_values must be a label or sequence of labels, got {type(subset_values)!r}"
            ) from exc
    wanted = {_normalize_group_label(v) for v in raw_list}
    wanted.discard(None)
    if not wanted:
        return pd.Series(False, index=obs_col.index)
    normed = obs_col.map(_normalize_group_label)
    return normed.isin(wanted)


def _merge_scatrans_uns(
    existing: dict[str, Any],
    meta: dict[str, Any],
    *,
    sticky_keys: tuple[str, ...] = ("raw_gene_list", "raw_gene_list_full", "history"),
) -> dict[str, Any]:
    """Merge run metadata into ``adata.uns['scatrans']`` without sticky stale fields.

    - Keys in ``sticky_keys`` are preserved from ``existing`` unless ``meta``
      provides a non-None replacement (e.g. updated ``history``).
    - All other keys from a previous run are **dropped**, then ``meta`` is applied.
    - ``meta`` values that are ``None`` **remove** the key (explicit "feature off")
      so a second analysis cannot inherit e.g. ``sample_col`` from a prior mixed-model run.
    """
    out: dict[str, Any] = {}
    for k in sticky_keys:
        if k in existing and k not in meta or k in existing and meta.get(k) is None:
            out[k] = existing[k]
    for k, v in meta.items():
        if v is None:
            out.pop(k, None)
        else:
            out[k] = v
    # history: prefer meta if present (caller usually already merged)
    if "history" in meta and meta["history"] is not None:
        out["history"] = meta["history"]
    elif "history" in existing:
        out["history"] = existing["history"]
    return out


def _validate_group_contrast(
    obs_col: pd.Series,
    *,
    groupby: str,
    target_group: str,
    reference_group: str,
) -> tuple[str, str, pd.Series]:
    """Validate target/reference exist; return normalized labels and per-cell series."""
    target_norm = _normalize_group_label(target_group)
    reference_norm = _normalize_group_label(reference_group)
    if target_norm is None or reference_norm is None:
        raise ValueError(
            f"target_group and reference_group must be valid labels for adata.obs['{groupby}']."
        )
    if target_norm.lower() == "nan" or reference_norm.lower() == "nan":
        raise ValueError(
            "target_group and reference_group cannot be the string 'nan' "
            f"(check missing values in adata.obs['{groupby}'])."
        )
    if target_norm == reference_norm:
        raise ValueError("target_group and reference_group must be different.")

    norm_groups = obs_col.map(_normalize_group_label)
    n_missing = int(norm_groups.isna().sum())
    if n_missing:
        logger.warning(
            "%d cells have missing %s labels and will be excluded from the contrast.",
            n_missing,
            groupby,
        )

    unique_valid = set(norm_groups.dropna().unique())
    if target_norm not in unique_valid:
        raise ValueError(
            f"target_group '{target_group}' not found in adata.obs['{groupby}']. "
            f"Available (non-missing): {sorted(unique_valid)[:20]}"
        )
    if reference_norm not in unique_valid:
        raise ValueError(
            f"reference_group '{reference_group}' not found in adata.obs['{groupby}']. "
            f"Available (non-missing): {sorted(unique_valid)[:20]}"
        )
    return target_norm, reference_norm, norm_groups


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
    "_warn_if_negative_layer_values",
    "_warn_if_not_integer_counts_matrix",
    "_warn_if_low_counts_matrix",
    "_safe_add_matrices",
    "_normalize_velocity_layers_by_size_factor",
    "_get_group_mean",
    "_get_exponential_scale_lambda",
    "_soft_scale",
    "_score_direction_effect",
    "_lambda_pval_for_active_score",
    "_composite_active_score_terms",
    "_pseudobulk_with_layers",
    "_fit_huber_bias_correction",
    "_resolve_aligned_raw_counts",
    "_x_looks_log_normalized",
    "_clear_log_preprocess_metadata",
    "_reconcile_log1p_marker",
    "_apply_de_preprocess",
    "_prepare_log_normalized_expression",
]


def _warn_if_negative_layer_values(layer: Any, layer_name: str, *, max_check: int = 100000) -> None:
    """Warn when a count layer contains negative values (physically invalid for RNA counts)."""
    if sparse.issparse(layer):
        vals = np.asarray(layer.data, dtype=float)
    else:
        vals = np.asarray(layer, dtype=float).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    if vals.size > max_check:
        rng = np.random.default_rng(0)
        vals = rng.choice(vals, size=max_check, replace=False)
    if np.any(vals < 0):
        n_neg = int(np.sum(vals < 0))
        logger.warning(
            "Layer '%s' contains %d negative value(s) (min=%.4g). "
            "RNA count layers should be non-negative; results may be unreliable. "
            "Check for data corruption or incorrect layer assignment.",
            layer_name,
            n_neg,
            float(np.min(vals)),
        )


def _is_integer_counts_like(X: Any, max_check: int = 100000, atol: float = 1e-6) -> bool:
    """Return True if the matrix contains non-negative values that are integer-valued
    (within tolerance). This is tolerant of float64 summed counts that are exactly
    (or very nearly) integers, which commonly occurs after pseudobulk aggregation.

    For large matrices (> max_check elements), a fixed-seed (0) random subsample is
    used for performance. This is a QC heuristic; it does not affect scientific results
    of the main analysis. Subsampling is deterministic for reproducibility.
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
        # Stride + random subsample: stride covers the full matrix deterministically;
        # random supplement catches sparse non-integer contamination missed by stride alone.
        stride = max(1, vals.size // (max_check // 2))
        stride_vals = vals[::stride]
        n_random = max_check - stride_vals.size
        if n_random > 0:
            random_vals = rng.choice(vals, size=min(n_random, vals.size), replace=False)
            vals = np.concatenate([stride_vals, random_vals])
        else:
            vals = stride_vals[:max_check]

    # Tolerant check: allows tiny floating point noise from summation / cast
    rounded = np.round(vals)
    return bool(np.all(vals >= 0) and np.allclose(vals, rounded, atol=atol, rtol=1e-5))


def _x_looks_zscore_scaled(finite: np.ndarray) -> bool:
    """Return True when values look z-score scaled (e.g. after ``sc.pp.scale``)."""
    if finite.size == 0 or not np.any(finite < 0):
        return False
    mx = float(np.nanmax(finite))
    mn = float(np.nanmin(finite))
    mean = float(np.nanmean(finite))
    std = float(np.nanstd(finite))
    if mx > 25.0:
        return False
    return abs(mean) < 1.0 and 0.25 < std < 5.0 and mn < -0.1


def _x_carries_library_size_signal(X: Any, *, min_cv: float = 0.12, max_cells: int = 5000) -> bool:
    """Return True when per-cell library sizes show meaningful sequencing-depth variation.

    Raw counts (including kallisto/salmon decimal UMIs) retain depth signal in row sums;
    normalize_total+log1p data does not. Used to avoid false 'already log-normalized'
    detection on small-magnitude decimal count matrices.
    """
    lib = _matrix_sum_axis1(X)
    lib = lib[np.isfinite(lib)]
    if lib.size < 2:
        return False
    if lib.size > max_cells:
        rng = np.random.default_rng(0)
        lib = rng.choice(lib, size=max_cells, replace=False)
    mean = float(np.mean(lib))
    if mean <= 0:
        return False
    cv = float(np.std(lib) / mean)
    return cv >= min_cv


def _x_gene_dispersion_looks_raw(
    X: Any, *, slope_threshold: float = 1.0, min_genes: int = 20, max_cells: int = 5000
) -> bool:
    """Return True when the cross-gene mean-variance relationship looks like raw counts.

    Raw RNA-seq counts (Poisson/NB) have variance that grows roughly linearly-to-quadratically
    with the mean across genes (log-log slope >= ~1). normalize_total+log1p compresses this
    relationship to a much flatter slope. Unlike per-cell library-size CV
    (``_x_carries_library_size_signal``), this is a cross-gene statistic and is not confounded
    by real biological heterogeneity between cells/cell types, which can inflate per-cell
    library-size CV even in properly log-normalized data (e.g. after ``anndata.concat`` drops
    the ``uns['log1p']`` marker). Returns False (conservatively "not raw-looking") when there
    are too few informative genes to estimate the relationship reliably.
    """
    if sparse.issparse(X):
        n_cells = X.shape[0]
        if n_cells > max_cells:
            rng = np.random.default_rng(0)
            idx = rng.choice(n_cells, size=max_cells, replace=False)
            X = X[idx]
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    else:
        arr = np.asarray(X, dtype=float)
        if arr.shape[0] > max_cells:
            rng = np.random.default_rng(0)
            idx = rng.choice(arr.shape[0], size=max_cells, replace=False)
            arr = arr[idx]
        mean = arr.mean(axis=0)
        mean_sq = (arr**2).mean(axis=0)
    var = np.clip(mean_sq - mean**2, 0, None)

    mask = np.isfinite(mean) & np.isfinite(var) & (mean > 0.05) & (var > 0)
    if int(mask.sum()) < min_genes:
        return False
    log_mean = np.log(mean[mask])
    log_var = np.log(var[mask])
    slope = float(np.polyfit(log_mean, log_var, 1)[0])
    return slope >= slope_threshold


def _x_looks_log_normalized(
    X: Any, *, max_check: int = 100000, has_log1p_marker: bool = False
) -> bool:
    """Return True when *X* is unlikely to be raw integer counts needing normalize+log1p."""
    if _is_integer_counts_like(X, max_check=max_check):
        return False
    arr = _dense_expression_matrix(X)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return False
    if _x_looks_zscore_scaled(finite):
        return False
    mx = float(np.nanmax(finite))
    has_neg = bool(np.any(finite < 0))
    if has_neg:
        return True
    # Trust an explicit uns['log1p'] marker for non-integer data even when max>20
    # (high-depth / bulk-like log matrices). Only reject when the matrix still
    # looks like raw library-size-dominated counts (stale marker after restore).
    if has_log1p_marker:
        # Trust marker unless matrix still looks raw + library-size dominated (stale marker).
        return not (
            mx > 20.0 and _x_carries_library_size_signal(X) and _x_gene_dispersion_looks_raw(X)
        )
    if mx > 20.0:
        return False
    # Per-cell library-size CV alone is unreliable: real biological heterogeneity (different
    # cell types/states with different expression breadth) can produce high CV even in data
    # that is already correctly log-normalized, most commonly after anndata.concat() drops the
    # uns['log1p'] marker. Require corroboration from the cross-gene mean-variance dispersion
    # check (robust to cell-level heterogeneity) before concluding "still needs normalization".
    return not (_x_carries_library_size_signal(X) and _x_gene_dispersion_looks_raw(X))


def _clear_log_preprocess_metadata(adata: ad.AnnData) -> None:
    """Drop scanpy log-transform markers after restoring raw counts into .X."""
    adata.uns.pop("log1p", None)


def _restore_log_from_scaled_x(adata: ad.AnnData) -> bool:
    """Reverse ``sc.pp.scale`` on ``.X`` when mean/std are stored in ``adata.var``."""
    arr = _dense_expression_matrix(adata.X)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0 or not _x_looks_zscore_scaled(finite):
        return False

    if "mean" not in adata.var.columns or "std" not in adata.var.columns:
        return False

    mean = pd.to_numeric(adata.var["mean"], errors="coerce").to_numpy(dtype=float)
    std = pd.to_numeric(adata.var["std"], errors="coerce").to_numpy(dtype=float)
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(std)):
        return False

    X = _dense_expression_matrix(adata.X)

    adata.X = X * std + mean
    logger.warning(
        ".X appeared z-score scaled (typical after sc.pp.scale). Restored log-normalized "
        "expression from adata.var['mean']/'std' for DE. Keep scaled data in a separate "
        "AnnData copy for PCA/clustering."
    )
    adata.uns.setdefault("log1p", {"base": None})
    return True


def _reconcile_log1p_marker(adata: ad.AnnData) -> bool:
    """Align ``uns['log1p']`` with the current ``.X`` scale.

    Mutates ``adata.uns`` when the marker is stale (common after ``restore_raw_counts``
    or manual reassignment of ``.X`` without clearing metadata).

    Returns True when callers should treat ``.X`` as already log-normalized and skip
    normalize_total + log1p.
    """
    has_marker = "log1p" in adata.uns
    arr = _dense_expression_matrix(adata.X)
    finite = arr[np.isfinite(arr)]
    x_is_scaled = _x_looks_zscore_scaled(finite) if finite.size else False
    x_is_log = _x_looks_log_normalized(adata.X, has_log1p_marker=has_marker)

    if has_marker and x_is_scaled:
        _clear_log_preprocess_metadata(adata)
        logger.warning(
            "Removed stale uns['log1p'] metadata: .X appears z-score scaled (e.g. after "
            "sc.pp.scale) while the log1p marker was still set. Downstream DE preprocessing "
            "will re-apply normalize_total + log1p when de_preprocess='auto'."
        )
        return False

    if has_marker and not x_is_log:
        _clear_log_preprocess_metadata(adata)
        logger.warning(
            "Removed stale uns['log1p'] metadata: .X appears to be raw or non-log "
            "counts while the log1p marker was still set. Downstream DE preprocessing "
            "will re-apply normalize_total + log1p when de_preprocess='auto'."
        )
        return False

    if x_is_log:
        if not has_marker:
            logger.warning(
                "DE preprocess: .X appears log-normalized without uns['log1p'] metadata; "
                "skipping re-normalization. If this is raw decimal counts (e.g. kallisto/salmon "
                "UMIs), pass de_preprocess='normalize_log1p' explicitly."
            )
        return True

    return False


def _apply_de_preprocess(
    adata: ad.AnnData,
    de_preprocess: str,
    *,
    skip_auto: bool = False,
) -> None:
    """Apply normalize_total + log1p according to ``de_preprocess`` mode.

    When ``skip_auto=True`` (count-based PyDESeq2 pseudobulk path), both
    ``auto`` *and* explicit ``normalize_log1p`` are skipped so counts are not
    log-transformed before DESeq2. ``normalize_log1p`` logs a warning in that case.
    """
    if de_preprocess in ("auto", "normalize_log1p") and not skip_auto:
        arr = _dense_expression_matrix(adata.X)
        finite = arr[np.isfinite(arr)]
        if finite.size and _x_looks_zscore_scaled(finite) and not _restore_log_from_scaled_x(adata):
            raise ValueError(
                ".X appears z-score scaled (e.g. after sc.pp.scale) but cannot be restored "
                "for DE (missing adata.var['mean']/'std']). Use an AnnData copy from before "
                "scaling, pass de_preprocess='none' with a suitable matrix, or store "
                "log-normalized expression in adata.raw or a layer."
            )

    if de_preprocess == "normalize_log1p":
        if skip_auto:
            logger.warning(
                "de_preprocess='normalize_log1p' is ignored on the count-based PyDESeq2 "
                "path (would log-transform integer counts). Leaving .X untransformed; "
                "pass de_preprocess='none' or 'auto' to silence this message."
            )
            return
        logger.info("DE preprocessing: applying normalize_total + log1p (explicit).")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.uns["log1p"] = {"base": None}
    elif de_preprocess == "auto" and not skip_auto:
        if _reconcile_log1p_marker(adata):
            logger.debug("DE preprocessing: 'auto' — .X already log-normalized, skipping.")
        else:
            logger.info(
                "DE preprocessing: 'auto' detected non-log .X; applying normalize_total + log1p."
            )
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            adata.uns["log1p"] = {"base": None}
    elif de_preprocess == "none":
        logger.debug("DE preprocessing: 'none' requested — no normalization applied.")
    elif de_preprocess == "auto" and skip_auto:
        logger.debug("DE preprocessing: 'auto' skipped for count-based pseudobulk backend.")


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
        # After pseudobulk aggregation, layers['counts'] is always rounded and
        # looks integer. Prefer the pre-aggregation verdict when available.
        if (
            require_integer
            and layer == "counts"
            and source_name == f"layers['{layer}']"
            and "pb_counts_is_count_like" in adata.uns
        ):
            if not bool(adata.uns["pb_counts_is_count_like"]):
                logger.warning(
                    "Counts from %s were aggregated from a non-count-like source "
                    "(pb_counts_is_count_like=False); skipping for count-based DE.",
                    source_name,
                )
                continue
        elif require_integer and not _is_integer_counts_like(mat):
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
    """Dense log1p library-size normalized matrix for mixed models (no double log1p).

    This is an internal helper for the LMM (mixedlm) DE path only.

    Logic (best-effort):
    1. Reconcile ``uns['log1p']`` with ``.X`` scale; if already log-normalized, return as-is.
    2. Else if matrix looks like integer counts, run normalize_total + log1p.
    3. Else if ``_x_looks_log_normalized`` is True, return as-is.
    4. Else apply log1p with a warning (large non-integer values).

    For rank_genes_groups / pydeseq2 etc. use the documented ``de_preprocess`` parameter,
    which shares the same reconciliation helpers via :func:`_apply_de_preprocess`.
    """
    ad_work = ad_expr.copy()
    if _reconcile_log1p_marker(ad_work):
        return _dense_expression_matrix(ad_work.X)

    if _is_integer_counts_like(ad_work.X):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.normalize_total(ad_work, target_sum=1e4)
            sc.pp.log1p(ad_work)
        return _dense_expression_matrix(ad_work.X)

    X = _dense_expression_matrix(ad_work.X)
    if _x_looks_log_normalized(X):
        return X

    finite = X[np.isfinite(X)]
    if finite.size:
        mx = float(np.nanmax(finite))
        logger.warning(
            "Mixed model input is neither log-normalized nor integer counts; applying log1p "
            "for stability (max value=%.1f). If this is already log-scale, set de_preprocess='none' "
            "or pre-log the input.",
            mx,
        )
        return np.log1p(np.clip(X, 0, None))
    return X


def _warn_if_low_counts_matrix(X: Any, max_check: int = 100000) -> None:
    """Warn if max rounded count looks suspiciously low for raw counts input to DE.

    For very large matrices, inspects only a deterministic fixed-seed subsample
    (see _is_integer_counts_like for rationale).
    """
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
    row_totals = np.asarray(total_layer.sum(axis=1)).ravel().astype(float)
    positive = row_totals > 0

    if positive.sum() == 0:
        return uns_layer, spl_layer, row_totals, np.ones_like(row_totals, dtype=float)

    if target_sum is None:
        target_sum = float(np.median(row_totals[positive]))

    # Zero-total rows keep factor=1 (stay zero). Do not use target/1e-8, which
    # produces ~1e9 factors that amplify float noise if any near-zero totals appear.
    factors = np.ones_like(row_totals, dtype=float)
    factors[positive] = target_sum / row_totals[positive]

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
    return np.asarray(sub.mean(axis=0)).ravel()


def _get_exponential_scale_lambda(x: np.ndarray) -> float:
    """Data-adaptive soft-scale length: median of positive values / ln(2).

    **Not transportable across runs.** Changing which genes (or how many) enter
    the vector changes λ and therefore all soft-scaled scores even if a gene's
    own raw statistic is fixed. ``active_score`` 0–100 ranks are therefore
    *within-analysis relative* only.
    """
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


def _score_direction_effect(
    logFC: np.ndarray,
    *,
    mixedlm_coef: np.ndarray | None = None,
) -> np.ndarray:
    """Effect used to gate the significance leg of ``active_score``.

    Prefer MixedLM fixed-effect coefficient when present (``p_adj`` tests that
    coefficient, not the sample-aware log2FC). Non-finite coefs fall back to logFC.
    """
    logFC_arr = np.asarray(logFC, dtype=float)
    if mixedlm_coef is None:
        return logFC_arr
    coef = np.asarray(mixedlm_coef, dtype=float)
    if coef.shape != logFC_arr.shape:
        raise ValueError(
            f"mixedlm_coef shape {coef.shape} does not match logFC shape {logFC_arr.shape}"
        )
    return np.where(np.isfinite(coef), coef, logFC_arr)


def _lambda_pval_for_active_score(
    p_adj: np.ndarray,
    direction_effect: np.ndarray,
    *,
    floor: float = 1.0,
) -> float:
    """Exponential scale for the p-value leg, estimated on direction-positive genes only.

    Estimating λ on all genes lets strongly downregulated (tiny p_adj) inflate the
    scale and shrink s3 for true upregulated genes. Fall back to the full vector
    when fewer than two positive-direction genes have finite -log10(p).
    """
    p = np.asarray(p_adj, dtype=float)
    effect = np.asarray(direction_effect, dtype=float)
    neglog = -np.log10(p + 1e-300)
    up = np.isfinite(effect) & (effect > 0.0) & np.isfinite(neglog)
    if int(up.sum()) >= 2:
        return max(_get_exponential_scale_lambda(neglog[up]), floor)
    return max(_get_exponential_scale_lambda(neglog), floor)


def _composite_active_score_terms(
    logFC: np.ndarray,
    residual: np.ndarray,
    p_adj: np.ndarray,
    lambda_fc: float,
    lambda_res: float,
    lambda_pval: float,
    *,
    direction_effect: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Soft-scale the three active_score legs with consistent direction gating.

    - s1 (logFC): ``_soft_scale`` clips negatives → only upregulation contributes.
    - s2 (unspliced residual): one-sided positive excess (same clip); independent of DE
      direction (nascent excess can exist without mature-RNA upregulation).
    - s3 (significance): ``-log10(p_adj)`` is directionless; multiply by
      ``(direction_effect > 0)`` so strongly downregulated genes cannot earn mid-range
      composite scores from p-values alone. ``direction_effect`` defaults to logFC;
      MixedLM paths should pass the fixed-effect coefficient (what ``p_adj`` tests).
      Observed and permutation paths must use this helper so the null matches.
    """
    logFC_arr = np.asarray(logFC, dtype=float)
    effect = logFC_arr if direction_effect is None else np.asarray(direction_effect, dtype=float)
    s1 = _soft_scale(logFC_arr, lambda_fc)
    s2 = _soft_scale(residual, lambda_res)
    s3 = _soft_scale(-np.log10(np.asarray(p_adj, dtype=float) + 1e-300), lambda_pval)
    s3 = np.where(effect > 0.0, s3, 0.0)
    return s1, s2, s3


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

    # Check count-likeness on the pre-aggregation, pre-rounding source matrix.
    # Summing+rounding below always produces integer-valued output, so checking
    # after aggregation would always pass and defeat strict_pydeseq2_counts.
    # Apply the same pre-agg check to every requested layer (esp. "counts"), not
    # only to .X — PyDESeq2 often consumes layers['counts'] after aggregation.
    x_source_is_count_like = _is_integer_counts_like(X_source)
    layer_is_count_like: dict[str, bool] = {
        layer: bool(_is_integer_counts_like(adata.layers[layer])) for layer in layers
    }

    valid_mask = adata.obs[sample_col].notna() & adata.obs[groupby].notna()
    if not valid_mask.any():
        raise ValueError(
            f"No observations with both '{sample_col}' and '{groupby}' labels for pseudobulk."
        )
    n_invalid = int((~valid_mask).sum())
    if n_invalid:
        logger.warning(
            "Excluding %d cells with missing %s or %s labels from pseudobulk aggregation.",
            n_invalid,
            sample_col,
            groupby,
        )
        adata = adata[valid_mask].copy()

    group_df = adata.obs[[sample_col, groupby]].copy()
    sample_labels = group_df[sample_col].map(_normalize_group_label).astype(str)
    group_labels = group_df[groupby].map(_normalize_group_label).astype(str)
    if (
        sample_labels.isin(["nan", "None", ""]).any()
        or group_labels.isin(["nan", "None", ""]).any()
    ):
        raise ValueError("Invalid sample/group labels remain after normalization for pseudobulk.")

    # Use a per-run UUID separator (printable, no embedded NUL) that is vanishingly unlikely
    # to appear in real sample/group names. Avoids pandas str-concat truncation with \0 and
    # eliminates "||" injection / mis-split risk from earlier fragile separator.
    _pb_sep = f"__scAT_PB_{uuid.uuid4().hex}__"
    # Force plain str (AnnData obs often stores Categorical; Categorical + str raises TypeError).
    pb_key = sample_labels + _pb_sep + group_labels
    unique_keys = pd.Index(pb_key.unique())

    X_rows, obs_rows = [], []
    layer_rows: dict[str, list] = {layer: [] for layer in layers}

    for key in unique_keys:
        mask = pb_key.values == key
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            continue
        x_sum = np.nan_to_num(_matrix_row_subset_sum_axis0(X_source, mask))
        # Clean to integer-valued floats for count-like data (pseudobulk sums).
        # Using round keeps the numeric value exact while float dtype is fine for AnnData/sparse.
        x_sum = np.round(x_sum).astype(np.float64, copy=False)

        if float(x_sum.sum()) < min_counts:
            continue

        sample_id, group_value = str(key).split(_pb_sep, 1)
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
            l_sum = np.nan_to_num(_matrix_row_subset_sum_axis0(adata.layers[layer], mask))
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
        adata_pb.obs[sample_col] = adata_pb.obs[sample_col].astype(str)
        adata_pb.obs[groupby] = adata_pb.obs[groupby].astype(str)
        adata_pb.obs_names = adata_pb.obs[sample_col] + "_" + adata_pb.obs[groupby]
        for layer in layers:
            adata_pb.layers[layer] = sparse.vstack(layer_rows[layer]).tocsr()
        adata_pb.obs_names_make_unique()
        # Carry the pre-aggregation count-likeness verdict so downstream consumers
        # (e.g. PyDESeq2 strict_pydeseq2_counts check) don't re-check the already-rounded X.
        adata_pb.uns["pb_x_is_count_like"] = bool(x_source_is_count_like)
        adata_pb.uns["pb_x_source_desc"] = x_source_name
        # Same trap for layers: post-aggregation np.round makes every layer look
        # integer-valued. Record pre-agg verdicts for each aggregated layer.
        adata_pb.uns["pb_layer_is_count_like"] = dict(layer_is_count_like)
        if "counts" in layer_is_count_like:
            # Convenience alias used by PyDESeq2 counts= / layers['counts'] paths.
            adata_pb.uns["pb_counts_is_count_like"] = bool(layer_is_count_like["counts"])
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
      - "bias_corrected": bool — True only when Huber length+intron regression succeeded
      - "method": the effective method used ("huber_length_intron" or "none")
      - "n_genes_used_for_fit": int
      - "fallback_to_median": bool — median centering used instead of Huber regression
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
            # Genes with expression but missing length/intron features cannot use
            # the multi-covariate prediction. Use the same median-centering as the
            # global fallback so they are not silently left at residual=0 (which
            # looks like "no excess" rather than "not annotated").
            missing_feat = valid_expr & ~valid_feat
            n_missing_feat = int(missing_feat.sum())
            bias_info["n_genes_residual_missing_features"] = n_missing_feat
            if n_missing_feat > 0:
                med = float(np.nanmedian(delta_velocity[valid_expr]))
                residual[missing_feat] = delta_velocity[missing_feat] - med
                bias_info["missing_features_residual"] = "median_centered_delta"
        except (ValueError, TypeError, np.linalg.LinAlgError, ArithmeticError) as e:
            logger.warning("Bias correction failed. Falling back to median. Reason: %s", e)
            bias_info["fallback_reason"] = (
                f"huber_regression_failed: {type(e).__name__}: {str(e)[:200]}"
            )

    if not regression_succeeded and valid_expr.sum() > 0:
        residual[valid_expr] = delta_velocity[valid_expr] - np.nanmedian(delta_velocity[valid_expr])
        bias_info["fallback_to_median"] = True
        # Median centering is a safe fallback but is not length/intron Huber correction.
        bias_info["bias_corrected"] = False
        bias_info["n_genes_residual_missing_features"] = int((valid_expr & ~valid_feat).sum())

    residual[~valid_expr] = 0.0
    return residual, bias_info


# warnings imported at top of file (used inside _fit_huber_bias_correction and pseudobulk creation)
