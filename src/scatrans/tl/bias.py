"""Abundance-/length-normalized unspliced-excess residual.

The shipped ``unspliced_excess_residual`` regresses the excess on intron number
and gene length, but that does not remove the **abundance / nuclear-retention**
confounder: nuclear-retained lncRNAs (e.g. *MALAT1*) and very long genes
dominate the top of the ranking even when they are not induced (their
``unspliced_excess_delta`` is huge in absolute terms). This module adds an
**additive** post-hoc residual (it does not modify ``active_score``):

* ``method="abundance"`` (default) — scale-free excess:
  ``delta / (total_us_counts + quantile(total_us_counts, floor_quantile))``.
  The abundance floor (default 75th percentile) deflates extreme-abundance
  outliers like MALAT1 while adding little noise to low-abundance genes.
* ``method="abundance_length"`` — the above, then a **gentle** robust
  residualization on ``log1p(gene_length)`` and ``log1p(intron_number)`` to also
  suppress length-dominated artifacts (e.g. very long *Ptpr* genes). Keep the
  length term mild (larger ``length_alpha``) so genuinely long induced genes are
  not over-penalized.

The result is written to ``unspliced_excess_residual_abnorm``.

.. note::
   This improves the **interpretability** of the residual ranking (removes
   length/abundance artifacts from the top); empirically it does not change the
   residual's *reliability* on steady-state velocity snapshots, where the
   nascent proxy is limited for kinetic reasons. See also
   :mod:`scatrans.tl.adaptive`.
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RESID_COL = "unspliced_excess_residual_abnorm"
_METHODS = ("abundance", "abundance_length")


def _impute_log1p(series: pd.Series) -> np.ndarray:
    v = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    finite = np.isfinite(v)
    v = np.where(finite, v, np.median(v[finite])) if finite.any() else np.zeros_like(v)
    return np.log1p(v)


def _standardize(x: np.ndarray) -> np.ndarray:
    return (x - np.median(x)) / (x.std() + 1e-9)


def add_abundance_normalized_residual(
    all_results: pd.DataFrame,
    *,
    method: str = "abundance",
    floor_quantile: float = 0.75,
    length_alpha: float = 1.0,
    inplace: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add ``unspliced_excess_residual_abnorm`` to an ``active_score`` table.

    Parameters
    ----------
    all_results
        Table from :func:`active_score` / :func:`run_default_pipeline` (needs
        ``unspliced_excess_delta``, ``total_us_counts``, ``valid_expr``; plus
        ``gene_length`` and ``intron_number`` for ``method="abundance_length"``).
    method
        ``"abundance"`` (default) or ``"abundance_length"``.
    floor_quantile
        Abundance floor as a quantile of ``total_us_counts`` (default 0.75).
    length_alpha
        L2 strength of the Huber length regression (larger = milder length
        correction). Only used for ``method="abundance_length"``.
    inplace
        Modify and return the input frame instead of a copy.

    Returns
    -------
    (all_results, diagnostics)
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {_METHODS}, got {method!r}")
    need = ["unspliced_excess_delta", "total_us_counts", "valid_expr"]
    if method == "abundance_length":
        need += ["gene_length", "intron_number"]
    missing = [c for c in need if c not in all_results.columns]
    if missing:
        raise KeyError(f"all_results missing required columns: {missing}")
    if not 0.0 <= floor_quantile <= 1.0:
        raise ValueError("floor_quantile must be in [0, 1]")

    ar = all_results if inplace else all_results.copy()
    expr = ar[ar["valid_expr"] == True]  # noqa: E712
    y = pd.to_numeric(expr["unspliced_excess_delta"], errors="coerce").to_numpy(float)
    ab = pd.to_numeric(expr["total_us_counts"], errors="coerce").to_numpy(float)
    ab = np.where(np.isfinite(ab), ab, 0.0)

    floor = float(np.quantile(ab, floor_quantile))
    v = y / (ab + floor + 1e-9)  # abundance-normalized excess

    if method == "abundance_length":
        X = np.column_stack(
            [_impute_log1p(expr["gene_length"]), _impute_log1p(expr["intron_number"])]
        )
        X = (X - X.mean(0)) / (X.std(0) + 1e-9)
        try:
            from sklearn.linear_model import HuberRegressor

            pred = (
                HuberRegressor(epsilon=1.35, alpha=length_alpha, max_iter=500).fit(X, v).predict(X)
            )
        except Exception:
            Xi = np.column_stack([np.ones(len(v)), X])
            beta, *_ = np.linalg.lstsq(Xi, v, rcond=None)
            pred = Xi @ beta
        v = v - pred

    ar[RESID_COL] = np.nan
    ar.loc[expr.index, RESID_COL] = _standardize(v)

    diagnostics = {
        "method": method,
        "floor_quantile": floor_quantile,
        "abundance_floor": floor,
        "n_expressed": int(len(expr)),
        "length_alpha": length_alpha if method == "abundance_length" else None,
    }
    logger.info(
        "abundance-normalized residual (%s): floor=%.1f on %d genes -> %s",
        method,
        floor,
        len(expr),
        RESID_COL,
    )
    return ar, diagnostics
