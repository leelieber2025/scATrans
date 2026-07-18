"""Reliability-adaptive weighting of the unspliced-excess (nascent) leg.

The composite :func:`active_score` blends a differential-expression leg with the
unspliced-excess (nascent-transcription) leg. Empirically the nascent leg's
usefulness is **regime dependent**: on metabolic-labeling data or early
stimulation time points it tracks newly transcribed RNA well, but on
steady-state / late single-cell snapshots the reference-gamma excess can
*anti-correlate* with net induction, in which case a fixed weight drags the
composite below plain DE.

This module adds an **additive** wrapper (it does not modify ``active_score``)
that:

1. estimates a proxy **reliability** from the data itself, and
2. produces an ``adaptive_score`` whose nascent leg is weighted by that
   reliability — shrunk to 0 when the proxy is uninformative/anti-correlated,
   and up-weighted (>1) when it is highly reliable.

Reliability is the AUC of ``unspliced_excess_residual`` recovering the obvious
DE-induced genes (``logFC >= 1`` and ``p_adj < 0.05``).

.. note::
   The DE anchor rewards a nascent proxy that *corroborates* DE; it therefore
   cannot credit the proxy for ranking genes that DE misses. This is an
   accepted trade-off (its job here is robustness, not beyond-DE discovery).
   The anchor is isolated in :func:`_de_induced_anchor` so it can be swapped
   later. ``adaptive_score`` is a heuristic rank, not a calibrated FDR.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import Any, Union

import numpy as np
import pandas as pd
from scipy.stats import rankdata

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

K_DEFAULT = 4.0
W_MAX_DEFAULT = 2.0
_REQUIRED_COLS = ("active_score", "unspliced_excess_residual", "logFC", "p_adj", "valid_expr")

# An anchor is the "truth-ish" induced-gene set the proxy's reliability is scored
# against. Either the string "de" (the built-in DE anchor), a callable
# ``expr -> array/Series of 0/1``, or a precomputed boolean/int array/Series.
AnchorSpec = Union[str, Callable[[pd.DataFrame], Any], "pd.Series", np.ndarray, list]


def _auc(score: np.ndarray, label: np.ndarray) -> float:
    """Rank-based ROC-AUC (Mann-Whitney U); NaN if a class is empty."""
    score = np.asarray(score, dtype=float)
    label = np.asarray(label, dtype=int)
    ok = np.isfinite(score)
    score, label = score[ok], label[ok]
    n1 = int(label.sum())
    n0 = len(label) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(score)
    return float((r[label == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _de_induced_anchor(expr: pd.DataFrame) -> np.ndarray:
    """Boolean anchor of 'obviously induced' genes (swap point for the anchor)."""
    return ((expr["logFC"] >= 1.0) & (expr["p_adj"] < 0.05)).to_numpy(dtype=int)


def labeling_anchor(
    column: str = "new_log2fc", threshold: float = 1.0
) -> Callable[[pd.DataFrame], np.ndarray]:
    """Build a labeling-truth anchor callable ``expr -> {0,1}`` from a column.

    Use with metabolic-labeling data (e.g. scNT-seq / sci-fate) where a per-gene
    newly-transcribed-RNA log2FC is available (default column ``new_log2fc``):
    ``add_adaptive_score(df, anchor=labeling_anchor())``. Unlike the DE anchor,
    this credits the proxy for tracking genuinely newly-transcribed genes that
    plain DE misses — the empirically correct reliability signal on labeling data
    (see the scNT-seq gating result). The named column must be present on the
    ``valid_expr`` rows; missing/NaN values count as not-induced.
    """

    def _anchor(expr: pd.DataFrame) -> np.ndarray:
        if column not in expr.columns:
            raise KeyError(
                f"labeling_anchor: column {column!r} not in results; "
                "supply the labeling truth (e.g. merge new-RNA log2FC) first."
            )
        vals = pd.to_numeric(expr[column], errors="coerce").to_numpy()
        return (np.nan_to_num(vals, nan=-np.inf) >= threshold).astype(int)

    _anchor.__name__ = f"labeling_anchor[{column}>={threshold}]"
    return _anchor


def _resolve_anchor(anchor: AnchorSpec, expr: pd.DataFrame) -> tuple[np.ndarray, str]:
    """Resolve an anchor spec to a 0/1 array aligned to ``expr`` rows + a label."""
    if isinstance(anchor, str):
        if anchor != "de":
            raise ValueError(f"unknown string anchor {anchor!r}; use 'de' or pass a callable/array")
        return _de_induced_anchor(expr), "de"
    if callable(anchor):
        vec = anchor(expr)
        if isinstance(vec, pd.Series):
            vec = vec.reindex(expr.index).to_numpy()
        vec = np.nan_to_num(np.asarray(vec, dtype=float), nan=0.0).astype(int)
        if len(vec) != len(expr):
            raise ValueError(f"anchor callable returned {len(vec)} values, expected {len(expr)}")
        return vec, getattr(anchor, "__name__", "callable")
    # array-like / Series aligned to the valid_expr rows
    if isinstance(anchor, pd.Series):
        vec = anchor.reindex(expr.index).to_numpy()
    else:
        vec = np.asarray(anchor)
        if len(vec) != len(expr):
            raise ValueError(f"anchor array has {len(vec)} values, expected {len(expr)}")
    vec = np.nan_to_num(np.asarray(vec, dtype=float), nan=0.0).astype(int)
    return vec, "array"


def adaptive_weight(
    reliability: float, k: float = K_DEFAULT, w_max: float = W_MAX_DEFAULT
) -> float:
    """Map a reliability AUC to a nascent-leg weight in ``[0, w_max]``.

    ``w = clip(k * (reliability - 0.5), 0, w_max)`` — 0.5 (chance) -> 0,
    higher reliability -> larger weight (may exceed 1 so a strong proxy leads).
    """
    if reliability is None or reliability != reliability:  # NaN
        reliability = 0.5
    return float(np.clip(k * (reliability - 0.5), 0.0, w_max))


def add_adaptive_score(
    all_results: pd.DataFrame,
    *,
    k: float = K_DEFAULT,
    w_max: float = W_MAX_DEFAULT,
    anchor: AnchorSpec = "de",
    inplace: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add ``adaptive_score`` (+ ``adaptive_score_pct``) to an ``active_score`` table.

    Parameters
    ----------
    all_results
        The ``all_results`` DataFrame from :func:`active_score` /
        :func:`run_default_pipeline` (needs columns ``active_score``,
        ``unspliced_excess_residual``, ``logFC``, ``p_adj``, ``valid_expr``).
    k, w_max
        Slope and cap of the reliability→weight map (see :func:`adaptive_weight`).
    anchor
        The induced-gene set the proxy's reliability is scored against.
        ``"de"`` (default) uses the built-in DE anchor (``logFC>=1 & p_adj<0.05``);
        pass :func:`labeling_anchor` (or any callable ``expr -> 0/1`` / array /
        Series aligned to the ``valid_expr`` rows) to anchor on metabolic-labeling
        truth instead. On labeling data the DE anchor under-estimates reliability
        (it is dominated by fast IEGs with depleted unspliced signal); a labeling
        anchor is the empirically correct choice there.
    inplace
        Modify and return the input frame instead of a copy.

    Returns
    -------
    (all_results, diagnostics)
        ``diagnostics`` has ``reliability_auc``, ``w_proxy``, ``anchor``,
        ``n_anchor_induced``, ``n_anchor_de_induced`` (back-compat alias),
        ``n_expressed``, ``k``, ``w_max``, ``verdict``.
    """
    missing = [c for c in _REQUIRED_COLS if c not in all_results.columns]
    if missing:
        raise KeyError(f"all_results missing required columns: {missing}")

    ar = all_results if inplace else all_results.copy()
    expr = ar[ar["valid_expr"] == True]  # noqa: E712
    anchor_vec, anchor_label = _resolve_anchor(anchor, expr)
    reliability = _auc(expr["unspliced_excess_residual"].to_numpy(), anchor_vec)
    w_proxy = adaptive_weight(reliability, k=k, w_max=w_max)

    r_fc = rankdata(expr["logFC"].to_numpy())
    r_pv = rankdata(-np.log10(expr["p_adj"].clip(lower=1e-300).to_numpy()))
    r_px = rankdata(expr["unspliced_excess_residual"].to_numpy())
    raw = (r_fc + r_pv + w_proxy * r_px) / (2.0 + w_proxy)
    span = raw.max() - raw.min()
    scaled = 100.0 * (raw - raw.min()) / (span + 1e-12)

    ar["adaptive_score"] = np.nan
    ar.loc[expr.index, "adaptive_score"] = scaled
    ar["adaptive_score_pct"] = ar["adaptive_score"].rank(pct=True) * 100.0

    if w_proxy == 0.0:
        verdict = "nascent leg DISABLED (uninformative/anti-correlated); adaptive == DE"
    elif w_proxy < 1.0:
        verdict = f"nascent leg down-weighted (w={w_proxy:.2f})"
    else:
        verdict = f"nascent leg up-weighted / leading (w={w_proxy:.2f})"
    n_anchor = int(anchor_vec.sum())
    diagnostics = {
        "reliability_auc": reliability,
        "w_proxy": w_proxy,
        "anchor": anchor_label,
        "n_anchor_induced": n_anchor,
        "n_anchor_de_induced": n_anchor,  # back-compat alias
        "n_expressed": int(len(expr)),
        "k": k,
        "w_max": w_max,
        "verdict": verdict,
    }
    logger.info(
        "adaptive_score: reliability AUC=%.3f (anchor=%s, n=%d) -> w_proxy=%.2f | %s",
        reliability if reliability == reliability else float("nan"),
        anchor_label,
        n_anchor,
        w_proxy,
        verdict,
    )
    return ar, diagnostics


def adaptive_active_score(
    adata: Any = None,
    *,
    all_results: pd.DataFrame | None = None,
    groupby: str = "condition",
    target_group: str | None = None,
    reference_group: str | None = None,
    organism: str = "mouse",
    k: float = K_DEFAULT,
    w_max: float = W_MAX_DEFAULT,
    anchor: AnchorSpec = "de",
    add_gene_features: bool = False,
    **pipeline_kwargs: Any,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Score with :func:`run_default_pipeline`, then add the adaptive score.

    Supply either a precomputed ``all_results`` (fast, re-uses a prior run) or an
    ``adata`` to score first. Returns ``(all_results_with_adaptive, diagnostics)``.
    """
    if all_results is None:
        if adata is None:
            raise ValueError("provide either `adata` or `all_results`")
        if target_group is None or reference_group is None:
            raise ValueError(
                "target_group and reference_group are required when scoring from `adata`"
            )
        import scatrans as scat  # lazy: avoid import cycle

        if add_gene_features:
            adata = scat.add_gene_features(adata, organism=organism)
        result = scat.run_default_pipeline(
            adata,
            groupby=groupby,
            target_group=target_group,
            reference_group=reference_group,
            organism=organism,
            run_go_enrichment=False,
            **pipeline_kwargs,
        )
        all_results = result.all_results.copy()

    return add_adaptive_score(all_results, k=k, w_max=w_max, anchor=anchor, inplace=True)
