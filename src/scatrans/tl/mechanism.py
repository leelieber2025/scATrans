"""scatrans.tl.mechanism — transcription-vs-stabilization mechanism annotation.

Additive helpers for the "DE selects, proxy annotates" workflow (see
``run_default_pipeline(select_by="de")``). All operate on an ``active_score`` /
pipeline results table and never gate gene-list membership:

  * :func:`annotate_mechanism_class` — per-gene static mechanism label
    (transcription-driven / stabilization-driven / ambiguous) from the nascent
    residual. Per-gene accuracy is modest by design; pair with
    :func:`program_mechanism`.
  * :func:`threshold_sensitivity` — DE-selected list size and overlap versus a
    reference cut across padj / logFC grids.
  * :func:`program_mechanism` — threshold-free program-level inference: pool
    per-gene transcription support over each gene set and test against background
    (competitive Mann–Whitney + BH FDR).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from .._utils import (
    LEGACY_VELOCITY_RESIDUAL_COL,
    UNSPLICED_EXCESS_RESIDUAL_COL,
)
from .bias import RESID_COL as ABNORM_RESID_COL

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

SUPPORT_COL = "transcription_support"
CLASS_COL = "mechanism_class"
CONF_COL = "mechanism_confidence"

# Preference order for the residual used as the transcription-support signal:
# bias-corrected (demotes MALAT1/long-gene artifacts) > raw > legacy.
_RESID_PREFERENCE = (ABNORM_RESID_COL, UNSPLICED_EXCESS_RESIDUAL_COL, LEGACY_VELOCITY_RESIDUAL_COL)


def _pick_residual_col(df: pd.DataFrame, residual_col: str | None) -> str:
    if residual_col is not None:
        if residual_col not in df.columns:
            raise KeyError(f"residual_col={residual_col!r} not in results columns")
        return residual_col
    for c in _RESID_PREFERENCE:
        if c in df.columns:
            return c
    raise KeyError(
        "no unspliced-excess residual column found; expected one of "
        f"{_RESID_PREFERENCE}. Run active_score (optionally bias_method=...)."
    )


def _robust_z(x: np.ndarray) -> np.ndarray:
    """Median/MAD standardization (robust to the heavy residual tails)."""
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    scale = 1.4826 * mad if mad > 0 else (np.nanstd(x) + 1e-9)
    return (x - med) / scale


def annotate_mechanism_class(
    results: pd.DataFrame,
    *,
    residual_col: str | None = None,
    logfc_col: str = "logFC",
    class_threshold: float = 0.5,
    reliability: float = 1.0,
    inplace: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Add a STATIC transcription-vs-stabilization mechanism label per gene.

    Annotation only — does NOT change gene-list membership. The nascent
    unspliced-excess residual is robust-standardized into ``transcription_support``
    (a z-score); up-regulated genes are then labeled:

    - ``transcription-driven``   : support >= +``class_threshold``
    - ``stabilization-driven``   : support <= -``class_threshold``
    - ``ambiguous``              : in between
    - ``unclassified_down``      : logFC < 0 (validated contrast is up-regulation)
    - ``unknown``                : residual missing

    Parameters
    ----------
    results
        ``active_score`` / ``run_default_pipeline`` table (needs a residual column;
        bias-corrected ``unspliced_excess_residual_abnorm`` is preferred when present).
    residual_col
        Override the residual column; default auto-picks the best available.
    logfc_col
        Direction column. If absent, all genes are treated as up-regulated.
    class_threshold
        Support magnitude (robust-z units) separating a class call from ambiguous.
    reliability
        Dataset-level proxy reliability in [0, 1] (e.g. from the regime diagnosis);
        scales ``mechanism_confidence`` so steady-state runs report low confidence.
    inplace
        Modify and return the input frame instead of a copy.

    Returns
    -------
    (results, diagnostics)
    """
    if class_threshold < 0:
        raise ValueError("class_threshold must be >= 0")
    if not 0.0 <= reliability <= 1.0:
        raise ValueError("reliability must be in [0, 1]")

    rc = _pick_residual_col(results, residual_col)
    df = results if inplace else results.copy()
    resid = pd.to_numeric(df[rc], errors="coerce").to_numpy(float)
    support = _robust_z(resid)
    df[SUPPORT_COL] = support

    # Missing/NaN logFC is neither up nor down (stays "unknown"); only a finite
    # negative logFC is a down-regulated gene. When no logFC column exists, treat
    # all genes as up-regulated (the validated contrast).
    if logfc_col in df.columns:
        lfc = pd.to_numeric(df[logfc_col], errors="coerce").to_numpy(float)
        up = lfc > 0
        down = lfc < 0
    else:
        up = np.ones(len(df), bool)
        down = np.zeros(len(df), bool)
    finite = np.isfinite(support)
    label = np.full(len(df), "unknown", dtype=object)
    label[finite & down] = "unclassified_down"
    up_ok = finite & up
    label[up_ok & (support >= class_threshold)] = "transcription-driven"
    label[up_ok & (support <= -class_threshold)] = "stabilization-driven"
    label[up_ok & (np.abs(support) < class_threshold)] = "ambiguous"
    df[CLASS_COL] = label

    # confidence: how far past the threshold, squashed to (0,1], times reliability
    conf = np.clip(np.abs(support) / (class_threshold + 1.0), 0.0, 1.0) * reliability
    conf[~up_ok] = np.nan
    df[CONF_COL] = conf

    counts = pd.Series(label).value_counts().to_dict()
    diagnostics = {
        "residual_col": rc,
        "class_threshold": class_threshold,
        "reliability": reliability,
        "class_counts": counts,
        "n_classified": int(up_ok.sum()),
    }
    logger.info(
        "mechanism annotation (support=%s, thr=%.2f, reliability=%.2f): %s",
        rc,
        class_threshold,
        reliability,
        counts,
    )
    return df, diagnostics


def threshold_sensitivity(
    results: pd.DataFrame,
    *,
    padj_grid: Sequence[float] = (0.01, 0.05, 0.1),
    logfc_grid: Sequence[float] = (0.58, 1.0, 1.5),
    direction: str = "up",
    reference: tuple[float, float] = (0.05, 1.0),
) -> pd.DataFrame:
    """Tabulate DE-selected list size vs (padj, logFC) cutoffs.

    Re-runs the DE selection (``filter_active_genes(select_by="de")``) over the
    grid and reports, for each cell, the number selected and the Jaccard overlap
    with the ``reference`` (padj, logFC) cut — so a report can show robustness to
    the threshold instead of defending one number.

    Returns a long DataFrame: columns ``padj_cutoff``, ``logfc_cutoff``,
    ``n_selected``, ``jaccard_vs_reference``, ``is_reference``.
    """
    from .filter import filter_active_genes  # local import: avoid import cycle

    def _select(padj: float, logfc: float) -> pd.Index:
        sub = filter_active_genes(
            results,
            select_by="de",
            padj_cutoff=padj,
            logfc_cutoff=logfc,
            logfc_direction=direction,
        )
        return sub.index

    ref_padj, ref_logfc = reference
    ref_idx = _select(ref_padj, ref_logfc)
    rows = []
    for padj in padj_grid:
        for logfc in logfc_grid:
            idx = _select(padj, logfc)
            inter = len(idx.intersection(ref_idx))
            union = len(idx.union(ref_idx))
            rows.append(
                {
                    "padj_cutoff": padj,
                    "logfc_cutoff": logfc,
                    "n_selected": int(len(idx)),
                    "jaccard_vs_reference": (inter / union) if union else 1.0,
                    "is_reference": bool(
                        np.isclose(padj, ref_padj) and np.isclose(logfc, ref_logfc)
                    ),
                }
            )
    return pd.DataFrame(rows)


def program_mechanism(
    results: pd.DataFrame,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    support_col: str | None = None,
    residual_col: str | None = None,
    min_genes: int = 5,
    restrict_index: Sequence[str] | None = None,
    alpha: float = 0.05,
    exclude_other_programs: bool = True,
) -> pd.DataFrame:
    """THRESHOLD-FREE program-level transcription-vs-stabilization inference.

    For each gene set, pool the per-gene ``transcription_support`` and test it
    against the background (all other tested genes) with a competitive
    Mann–Whitney U. No per-gene hard cutoff — program-level pooling is preferred
    over single-gene mechanism claims.

    Parameters
    ----------
    results
        Results table. If it lacks ``transcription_support`` it is computed on the
        fly from the residual column (robust-z), same as
        :func:`annotate_mechanism_class`.
    gene_sets
        Mapping ``{program_name: [gene, ...]}`` (genes matched against
        ``results.index``).
    support_col
        Column holding the per-gene transcription support (default: auto).
    residual_col
        Residual column used to derive support when ``support_col`` is absent.
    min_genes
        Skip programs with fewer than this many tested genes.
    restrict_index
        Optional subset of genes (e.g. up-regulated / DE-tested) to use as the
        analysis universe for both foreground and background.
    alpha
        FDR level flagged in the ``significant`` column (BH across programs).
    exclude_other_programs
        When True (default), each program is compared against a "generic"
        background = tested genes in NONE of the provided gene sets. This avoids a
        strong program inflating another's background (which can make a null
        program read "stabilization-driven"). Set False for the plain competitive
        background (all tested genes outside the current set).

    Returns
    -------
    DataFrame per program: ``n_genes``, ``mean_support``, ``bg_mean_support``,
    ``direction`` (transcription-driven / stabilization-driven / ns), ``U``,
    ``p_value``, ``fdr``, ``significant``. Sorted by ``p_value``.
    """
    df = results
    if support_col is not None:
        if support_col not in df.columns:
            raise KeyError(f"support_col={support_col!r} not in results columns")
        support = pd.to_numeric(df[support_col], errors="coerce").to_numpy(float)
    elif SUPPORT_COL in df.columns:
        support = pd.to_numeric(df[SUPPORT_COL], errors="coerce").to_numpy(float)
    else:
        rc = _pick_residual_col(df, residual_col)
        support = _robust_z(pd.to_numeric(df[rc], errors="coerce").to_numpy(float))

    s = pd.Series(support, index=df.index)
    if restrict_index is not None:
        s = s.reindex(pd.Index(restrict_index).intersection(df.index))
    s = s[np.isfinite(s.to_numpy(float))]
    universe = s.index

    # union of all provided programs (for the generic background)
    all_set = universe.intersection(
        pd.Index([str(g) for genes in gene_sets.values() for g in genes])
    )

    rows = []
    for name, genes in gene_sets.items():
        in_set = universe.intersection(pd.Index([str(g) for g in genes]))
        n = len(in_set)
        if n < min_genes:
            continue
        fg = s.loc[in_set].to_numpy(float)
        # background: generic genes (in no program) when exclude_other_programs,
        # else all tested genes outside the current set.
        bg_idx = (
            universe.difference(all_set) if exclude_other_programs else universe.difference(in_set)
        )
        bg = s.loc[bg_idx].to_numpy(float)
        if len(bg) == 0 and exclude_other_programs:
            # The generic background is empty because every tested gene falls in
            # some program (common when the user passes a complete partition of the
            # DE list). Fall back to the plain competitive background (all tested
            # genes outside this set) rather than silently dropping the program.
            logger.warning(
                "program_mechanism: generic background empty for %r (all tested genes "
                "are inside some program); falling back to competitive background.",
                name,
            )
            bg = s.loc[universe.difference(in_set)].to_numpy(float)
        if len(bg) == 0:
            continue
        U, p = stats.mannwhitneyu(fg, bg, alternative="two-sided")
        fg_mean, bg_mean = float(np.mean(fg)), float(np.mean(bg))
        rows.append(
            {
                "program": name,
                "n_genes": int(n),
                "mean_support": fg_mean,
                "bg_mean_support": bg_mean,
                "U": float(U),
                "p_value": float(p),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out.assign(direction=[], fdr=[], significant=[])

    # BH FDR across programs
    p = out["p_value"].to_numpy(float)
    order = np.argsort(p)
    ranked = p[order]
    m = len(p)
    fdr_sorted = ranked * m / (np.arange(1, m + 1))
    fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
    fdr = np.empty(m)
    fdr[order] = np.clip(fdr_sorted, 0, 1)
    out["fdr"] = fdr
    out["significant"] = out["fdr"] < alpha
    # Direction requires strict inequality; equal foreground/background mean support
    # (a tie) is not evidence for either mechanism → "ns".
    diff = out["mean_support"] - out["bg_mean_support"]
    out["direction"] = np.where(
        ~out["significant"] | (diff == 0),
        "ns",
        np.where(diff > 0, "transcription-driven", "stabilization-driven"),
    )
    return out.sort_values("p_value").reset_index(drop=True)
