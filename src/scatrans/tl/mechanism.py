"""scatrans.tl.mechanism — transcription-vs-stabilization mechanism annotation.

Additive helpers for the "DE selects, proxy annotates" workflow (see
``run_default_pipeline(select_by="de")``). All operate on an ``active_score`` /
pipeline results table and never gate gene-list membership:

  * :func:`annotate_mechanism_class` — per-gene static mechanism label
    (transcription-driven / stabilization-driven / ambiguous) from the nascent
    residual. Per-gene accuracy is modest by design; pair with
    :func:`program_mechanism` / :func:`program_mechanism_induction_matched`.
  * :func:`threshold_sensitivity` — DE-selected list size and overlap versus a
    reference cut across padj / logFC grids.
  * :func:`program_mechanism` — threshold-free program-level inference: pool
    per-gene transcription support over each gene set and test against background
    (competitive Mann–Whitney + BH FDR).
  * :func:`program_mechanism_induction_matched` — induction-controlled program
    tests (OLS + nearest-logFC neighbors) that avoid the high-induction snapshot
    confound that misleads naive per-gene class enrichment.
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
INDUCTION_CONFOUND_COL = "induction_confounded"

# Minimum number of induced genes needed to estimate the |logFC| percentile the
# induction-confound flag depends on; below this the flag is skipped (tiny tables,
# e.g. unit fixtures, are left untouched).
_MIN_INDUCED_FOR_FLAG = 50

# Below this reliability, hard per-gene mechanism labels are suppressed by default
# (regime high/low unspliced → proxy untrustworthy; confidence already ~0).
_DEFAULT_MIN_RELIABILITY_FOR_HARD_LABELS = 0.05

# Named parameter bundles for :func:`annotate_mechanism_class`. ``high_precision``
# trades recall for a lower transcription-driven false-positive rate: on scEU-seq
# RPE1 (matched-abundance gamma terciles) raising the threshold 0.5 -> 1.0 halves
# the transcription-driven FPR (0.187 -> 0.086) at the same precision. Callers can
# still override the field explicitly.
#
# NOTE: an earlier design also flipped on an asymmetric confidence discount for
# transcription-driven calls, on the premise that the proxy resolves stabilization
# more cleanly (per-gene mean robust-z stab -0.31 vs synth +0.05). A direct
# confusion-matrix check falsified that premise: at the operating point,
# transcription-driven precision (0.664) is HIGHER than stabilization-driven (0.591)
# on RPE1 — the symmetric +/-threshold sits off-centre from the negatively-shifted
# class pair, so the transcription cut is effectively stricter/purer. Down-weighting
# transcription confidence would penalise the MORE precise call, so it was dropped.
# (The residual precision asymmetry, if one wants to address it, calls for an
# asymmetric THRESHOLD, not a confidence discount.)
_PRESETS: dict[str, dict[str, Any]] = {
    "high_precision": {"class_threshold": 1.0},
}

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
    class_threshold: float | None = None,
    preset: str | None = None,
    reliability: float = 1.0,
    suppress_hard_labels_when_unreliable: bool = True,
    min_reliability_for_hard_labels: float = _DEFAULT_MIN_RELIABILITY_FOR_HARD_LABELS,
    flag_induction_confound: bool = True,
    induction_confound_quantile: float = 0.90,
    induction_confound_penalty: str = "graded",
    induction_logfc_floor: float = 1.0,
    padj_col: str = "p_adj",
    induction_padj_cutoff: float = 0.05,
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

    Per-gene accuracy is modest by design (proxy matched-abundance AUC ~0.61); the
    ``"high_precision"`` preset lowers the false-positive rate, and program-level
    pooling (:func:`program_mechanism` / :func:`program_mechanism_induction_matched`)
    is the real precision lever. Prefer program calls over single-gene mechanism
    claims.

    When ``reliability`` is below ``min_reliability_for_hard_labels`` (default from
    regime extremes such as snRNA / extreme high-unspliced libraries), hard
    transcription/stabilization labels are reassigned to ``ambiguous`` by default
    so a near-zero-confidence proxy cannot mint decisive-looking classes.

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
        Defaults to the ``preset`` value, or ``0.5`` when no preset is given.
        An explicit value always wins over the preset.
    preset
        Named parameter bundle. ``"high_precision"`` sets ``class_threshold=1.0``,
        trading recall for a lower transcription-driven false-positive rate (halved
        on scEU-seq RPE1 at the same precision). An explicit ``class_threshold``
        overrides the preset.
    reliability
        Dataset-level proxy reliability in [0, 1] (e.g. from the regime diagnosis);
        scales ``mechanism_confidence`` so steady-state runs report low confidence.
    suppress_hard_labels_when_unreliable
        When True (default), if ``reliability < min_reliability_for_hard_labels``,
        rewrite per-gene ``transcription-driven`` / ``stabilization-driven`` labels
        to ``ambiguous`` (support and continuous confidence still computed; down /
        unknown unchanged). Set False to keep hard labels at extreme regimes.
    min_reliability_for_hard_labels
        Reliability floor for hard labels when suppression is enabled (default 0.05).
    flag_induction_confound
        Mark and down-weight per-gene stabilization calls in the high-induction
        regime, where the static residual is single-snapshot NON-identifiable.
        Strong induction inflates FALSE stabilization: at a snapshot the mature pool
        accumulates faster than a resting U/S reference predicts, pushing the residual
        negative even for transcription-driven genes (e.g. ISGs; documented as the
        Fig. 2 time-decay). The confound overlaps genuine high-induction stabilization
        (e.g. ARE cytokines) on every per-gene observable, so this is an HONEST
        reliability flag, not a corrector: it never relabels and never consults the
        nascent detection score (the axes stay decoupled). It adds a boolean
        ``induction_confounded`` column and multiplies ``mechanism_confidence`` by a
        penalty for flagged genes; program-level calls (:func:`program_mechanism`,
        which run on ``transcription_support``) are unaffected. Requires ``logfc_col``
        and at least ``50`` induced genes; otherwise it is skipped.
        WARNING — ``induction_confounded=True`` marks a *danger zone*, NOT a correction:
        it does **not** mean the gene is "not stabilization-driven" (it flags both
        mislabeled transcription-driven ISGs and genuine high-induction ARE, which are
        indistinguishable per-gene). Resolve mechanism at the induction-matched
        program level, not by reading the boolean. NOTE — this defaults to ``True``, so
        it is a soft behaviour change: ``mechanism_confidence`` (not the label or
        support) is lowered for flagged genes relative to older releases; set ``False``
        to restore the prior confidence values.
    induction_confound_quantile
        A stabilization call is flagged when its ``|logFC|`` is at/above this quantile
        **within the induced genes** (``|logFC| >= induction_logfc_floor``). Default
        ``0.90`` (the top-decile of induction magnitude).
    induction_confound_penalty
        Confidence-discount shape for flagged stabilization calls. ``"graded"``
        (default): ``1 - 0.7 * clip((pct - q)/(0.99 - q), 0, 1)`` for genes in the
        flagged zone (1.0 at the quantile down to a 0.3 floor at the 99th percentile),
        leaving lower-induction stabilization calls untouched. ``"smooth"``: a
        thresholdless quadratic ramp ``1 - 0.7 * clip((pct - 0.5)/0.5, 0, 1)**2``
        applied across all induced stabilization calls (mid-induction calls are
        discounted mildly, extreme ones strongly). Note: under ``"smooth"`` a
        mid-induction gene can have discounted confidence while ``induction_confounded``
        is ``False`` — the boolean always marks the top-decile high-risk zone
        (``induction_confound_quantile``), whereas the smooth penalty starts below it.
        The penalty constants (0.7 weight, 0.3 floor, quantile 0.90) are heuristic
        defaults chosen from a threshold/shape sweep on GSE226488 + Kang, not a
        calibrated model; tune them for a specific dataset if needed.
    induction_logfc_floor
        ``|logFC|`` bar defining the "induced" reference population for the percentile
        (default ``1.0``, the standard effect-size gate). Prevents weakly-changed genes
        from setting the induction scale.
    inplace
        Modify and return the input frame instead of a copy.

    Returns
    -------
    (results, diagnostics)
    """
    preset_cfg: dict[str, Any] = {}
    if preset is not None:
        if preset not in _PRESETS:
            raise ValueError(f"unknown preset {preset!r}; choose from {sorted(_PRESETS)} or None")
        preset_cfg = _PRESETS[preset]
    if class_threshold is None:
        class_threshold = preset_cfg.get("class_threshold", 0.5)

    if class_threshold < 0:
        raise ValueError("class_threshold must be >= 0")
    if not 0.0 <= reliability <= 1.0:
        raise ValueError("reliability must be in [0, 1]")
    if not 0.0 <= min_reliability_for_hard_labels <= 1.0:
        raise ValueError("min_reliability_for_hard_labels must be in [0, 1]")
    if induction_confound_penalty not in ("graded", "smooth"):
        raise ValueError("induction_confound_penalty must be 'graded' or 'smooth'")
    if not 0.0 < induction_confound_quantile < 1.0:
        raise ValueError("induction_confound_quantile must be in (0, 1)")

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

    # Item H: extreme regimes (reliability ~ 0) — do not mint hard classes.
    n_hard_suppressed = 0
    hard_suppressed = bool(
        suppress_hard_labels_when_unreliable and reliability < min_reliability_for_hard_labels
    )
    if hard_suppressed:
        hard = (label == "transcription-driven") | (label == "stabilization-driven")
        n_hard_suppressed = int(hard.sum())
        label[hard] = "ambiguous"
        logger.info(
            "mechanism: reliability=%.3f < %.3f — suppressed %d hard per-gene "
            "labels to ambiguous (support unchanged; set "
            "suppress_hard_labels_when_unreliable=False to keep hard labels).",
            reliability,
            min_reliability_for_hard_labels,
            n_hard_suppressed,
        )

    df[CLASS_COL] = label

    # confidence: how far past the threshold, squashed to (0,1], times reliability
    conf = np.clip(np.abs(support) / (class_threshold + 1.0), 0.0, 1.0) * reliability
    conf[~up_ok] = np.nan

    # Induction-confound flag: strong induction inflates FALSE stabilization at a
    # snapshot (mature accumulates -> residual down), so a high-induction stabilization
    # call is single-snapshot non-identifiable. Mark it and discount its confidence;
    # never relabel, never use nascent z (axes stay decoupled). See docstring.
    confounded = np.zeros(len(df), bool)
    n_confounded = 0
    if flag_induction_confound and logfc_col in df.columns:
        absfc = np.abs(lfc)
        is_stab = label == "stabilization-driven"
        # Reference population that sets the induction scale = significant induced
        # genes (the DE program). Restricting to padj-significant genes is essential:
        # ranking against ALL |logFC|>=floor genes lets noisy low-expression genes
        # with huge non-significant fold-changes inflate the tail and hide the truly
        # extreme (e.g. ISG15). A gene's percentile is its |logFC| ECDF against that
        # reference distribution, so genes outside the reference are still scored.
        ref = up_ok & np.isfinite(absfc) & (absfc >= induction_logfc_floor)
        if padj_col in df.columns:
            padj = pd.to_numeric(df[padj_col], errors="coerce").to_numpy(float)
            ref = ref & np.isfinite(padj) & (padj < induction_padj_cutoff)
        if int(ref.sum()) >= _MIN_INDUCED_FOR_FLAG:
            ref_vals = np.sort(absfc[ref])
            pct = np.full(len(df), np.nan)
            fin = np.isfinite(absfc)
            pct[fin] = np.searchsorted(ref_vals, absfc[fin], side="right") / len(ref_vals)
            q = induction_confound_quantile
            in_zone = is_stab & np.isfinite(pct) & (pct >= q)
            confounded = in_zone
            if induction_confound_penalty == "graded":
                r = np.clip((pct - q) / (0.99 - q), 0.0, 1.0)
                pen = np.where(in_zone, 1.0 - 0.7 * r, 1.0)
            else:  # "smooth": thresholdless quadratic ramp from the reference median
                ramp = np.clip((pct - 0.5) / 0.5, 0.0, 1.0)
                pen = np.where(is_stab & np.isfinite(pct), 1.0 - 0.7 * ramp**2, 1.0)
            pen = np.where(np.isfinite(pen), pen, 1.0)
            conf = conf * pen
            n_confounded = int(confounded.sum())
    df[INDUCTION_CONFOUND_COL] = confounded
    df[CONF_COL] = conf

    counts = pd.Series(label).value_counts().to_dict()
    diagnostics = {
        "residual_col": rc,
        "class_threshold": class_threshold,
        "preset": preset,
        "reliability": reliability,
        "suppress_hard_labels_when_unreliable": bool(suppress_hard_labels_when_unreliable),
        "min_reliability_for_hard_labels": float(min_reliability_for_hard_labels),
        "hard_labels_suppressed": hard_suppressed,
        "n_hard_labels_suppressed": n_hard_suppressed,
        "class_counts": counts,
        "n_classified": int(up_ok.sum()),
        "flag_induction_confound": bool(flag_induction_confound),
        "induction_confound_penalty": induction_confound_penalty,
        "induction_confound_quantile": induction_confound_quantile,
        "n_induction_confounded": n_confounded,
    }
    logger.info(
        "mechanism annotation (support=%s, thr=%.2f, preset=%s, reliability=%.2f): %s",
        rc,
        class_threshold,
        preset,
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

    .. warning::
       **Use mechanism-coherent, length-aware gene sets — not arbitrary functional
       pathways.** Pooling sharpens the weak per-gene mechanism signal only when the
       set's genes actually share a mechanism *and* are not confounded with gene
       length. The unspliced-excess support is length-entangled, and
       functionally-defined libraries (KEGG / GO) differ enormously in gene-length
       composition, so a naive pathway screen ranks pathways by **gene length, not
       mechanism**. Validated on scEU-seq RPE1 against independent pulse-chase
       ``gamma``: pathway ``mean_support`` correlated with the labeling truth
       (Spearman +0.39) but that vanished after controlling for gene length
       (partial Spearman −0.06, p≈0.5); the "most transcription-driven" pathways
       were simply the longest-gene ones (axon guidance, LTP — not even expressed in
       the epithelial line). Prefer **curated mechanism sets** (e.g. ARE-containing
       transcripts) and/or a **length-matched background**; treat KEGG/GO
       ``program_mechanism`` output as descriptive only. The against-background test
       is also underpowered for impure sets — it detects strong curated programs
       (e.g. ARE stabilization) but not weak, mixed ones.

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


def _winsorize_1d(x: np.ndarray, lo_pct: float, hi_pct: float) -> np.ndarray:
    """Clip finite values to the [lo_pct, hi_pct] percentiles (NaN preserved)."""
    out = np.asarray(x, dtype=float).copy()
    fin = np.isfinite(out)
    if int(fin.sum()) < 3:
        return out
    lo, hi = np.percentile(out[fin], [lo_pct, hi_pct])
    out[fin] = np.clip(out[fin], lo, hi)
    return out


def _ols_coef_se_p(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """OLS with classical SE and two-sided t p-values. X includes intercept column."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = max(int(len(y) - X.shape[1]), 1)
    sigma2 = float(resid @ resid) / dof
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        tval = np.where(se > 0, beta / se, np.nan)
    pval = 2.0 * stats.t.sf(np.abs(tval), dof)
    return beta, se, pval


def program_mechanism_induction_matched(
    results: pd.DataFrame,
    gene_sets: Mapping[str, Sequence[str]],
    *,
    support_col: str | None = None,
    residual_col: str | None = None,
    logfc_col: str = "logFC",
    padj_col: str = "p_adj",
    padj_cutoff: float = 0.05,
    logfc_min: float = 0.0,
    min_genes: int = 5,
    winsorize_support_pct: tuple[float, float] = (1.0, 99.0),
    nearest_k: int = 5,
    methods: Sequence[str] = ("regression", "nearest"),
) -> pd.DataFrame:
    """Induction-controlled program mechanism tests (preferred for claims).

    Naive competitive pooling (:func:`program_mechanism`) and especially
    per-gene ``mechanism_class`` ORA can be **induction-confounded** at a single
    snapshot: highly induced genes accumulate mature mRNA and read
    stabilization-driven regardless of true mechanism. This function tests each
    curated program **at matched induction strength** on the induced universe.

    Methods (select with ``methods``; both run by default):

    * ``"regression"`` — OLS ``support ~ 1 + logFC + is_program`` on induced genes;
      the ``is_program`` coefficient is the induction-adjusted effect (negative ⇒
      more stabilization-weighted than induction alone predicts).
    * ``"nearest"`` — for each program gene, match ``k`` non-program induced genes
      by closest ``logFC``; paired Wilcoxon on ``support_prog - median(support_ctrl)``.

    Parameters
    ----------
    results
        Table with ``transcription_support`` (or a residual to derive it) and DE
        columns (``logFC``, ``p_adj``).
    gene_sets
        ``{program_name: [gene, ...]}`` — prefer mechanism-coherent curated sets,
        not arbitrary KEGG/GO pathways.
    padj_cutoff, logfc_min
        Induced universe: ``p_adj < padj_cutoff`` and ``logFC > logfc_min`` with
        finite support (default: significant up genes).
    winsorize_support_pct
        Percentile clip on support before tests (default 1–99; matches the
        GSE226488 induction-matched validation script). Set ``(0, 100)`` to skip.
    nearest_k
        Controls per program gene for the nearest-logFC test.
    methods
        Subset of ``{"regression", "nearest"}``.

    Returns
    -------
    DataFrame one row per program (programs below ``min_genes`` in the induced
    universe are omitted). Key columns include ``n_genes``, ``mean_support``,
    ``regression_beta`` / ``regression_p`` (if requested), ``nearest_median_delta``
    / ``nearest_p`` (if requested), and a summary ``direction`` /
    ``preferred_p`` from the first available method with a finite p-value
    (regression preferred when both run).
    """
    allowed = {"regression", "nearest"}
    methods_set = {str(m).lower() for m in methods}
    if not methods_set or not methods_set.issubset(allowed):
        raise ValueError(f"methods must be a non-empty subset of {sorted(allowed)}")
    if nearest_k < 1:
        raise ValueError("nearest_k must be >= 1")
    if logfc_col not in results.columns:
        raise KeyError(f"logfc_col={logfc_col!r} required for induction-matched tests")

    if support_col is not None:
        if support_col not in results.columns:
            raise KeyError(f"support_col={support_col!r} not in results columns")
        support = pd.to_numeric(results[support_col], errors="coerce")
    elif SUPPORT_COL in results.columns:
        support = pd.to_numeric(results[SUPPORT_COL], errors="coerce")
    else:
        rc = _pick_residual_col(results, residual_col)
        support = pd.Series(
            _robust_z(pd.to_numeric(results[rc], errors="coerce").to_numpy(float)),
            index=results.index,
        )

    logfc = pd.to_numeric(results[logfc_col], errors="coerce")
    if padj_col in results.columns:
        padj = pd.to_numeric(results[padj_col], errors="coerce")
        padj_ok = np.isfinite(padj.to_numpy(float)) & (padj.to_numpy(float) < padj_cutoff)
    else:
        logger.warning(
            "program_mechanism_induction_matched: %r missing — induced universe "
            "uses logFC only (no padj gate).",
            padj_col,
        )
        padj_ok = np.ones(len(results), bool)

    ok = (
        np.isfinite(support.to_numpy(float))
        & np.isfinite(logfc.to_numpy(float))
        & (logfc.to_numpy(float) > logfc_min)
        & padj_ok
    )
    uni = results.index[ok]
    if len(uni) < min_genes + 2:
        return pd.DataFrame(
            columns=[
                "program",
                "n_genes",
                "n_universe",
                "mean_support",
                "direction",
                "preferred_p",
            ]
        )

    S = support.reindex(uni).to_numpy(float)
    L = logfc.reindex(uni).to_numpy(float)
    if winsorize_support_pct is not None:
        lo_p, hi_p = winsorize_support_pct
        if not (0.0 <= lo_p < hi_p <= 100.0):
            raise ValueError("winsorize_support_pct must satisfy 0 <= lo < hi <= 100")
        if lo_p > 0.0 or hi_p < 100.0:
            S = _winsorize_1d(S, lo_p, hi_p)

    rows: list[dict[str, Any]] = []
    for name, genes in gene_sets.items():
        in_set = uni.intersection(pd.Index([str(g) for g in genes]))
        n = len(in_set)
        if n < min_genes:
            continue
        is_prog = np.asarray(uni.isin(in_set), dtype=bool)
        fg = S[is_prog]
        row: dict[str, Any] = {
            "program": name,
            "n_genes": int(n),
            "n_universe": int(len(uni)),
            "mean_support": float(np.mean(fg)),
            "mean_logfc": float(np.mean(L[is_prog])),
        }

        reg_p = np.nan
        if "regression" in methods_set:
            # need both classes present
            if is_prog.sum() >= min_genes and (~is_prog).sum() >= 3:
                X = np.column_stack([np.ones(len(uni)), L, is_prog.astype(float)])
                beta, se, pval = _ols_coef_se_p(X, S)
                row["regression_beta"] = float(beta[2])
                row["regression_se"] = float(se[2])
                row["regression_p"] = float(pval[2])
                row["regression_logfc_beta"] = float(beta[1])
                reg_p = float(pval[2])
            else:
                row["regression_beta"] = np.nan
                row["regression_se"] = np.nan
                row["regression_p"] = np.nan
                row["regression_logfc_beta"] = np.nan

        nn_p = np.nan
        if "nearest" in methods_set:
            ctrl_mask = ~is_prog
            if int(ctrl_mask.sum()) >= nearest_k and int(is_prog.sum()) >= 1:
                L_ctrl = L[ctrl_mask]
                S_ctrl = S[ctrl_mask]
                deltas = []
                for lg, sg in zip(L[is_prog], S[is_prog], strict=True):
                    d = np.abs(L_ctrl - lg)
                    k_idx = np.argpartition(d, min(nearest_k - 1, len(d) - 1))[:nearest_k]
                    ctrl = float(np.median(S_ctrl[k_idx]))
                    deltas.append(sg - ctrl)
                deltas_a = np.asarray(deltas, dtype=float)
                row["nearest_median_delta"] = float(np.median(deltas_a))
                row["nearest_mean_delta"] = float(np.mean(deltas_a))
                # two-sided Wilcoxon on paired deltas (zero-null)
                if len(deltas_a) >= 5 and np.any(deltas_a != 0):
                    try:
                        w = stats.wilcoxon(deltas_a, alternative="two-sided")
                        row["nearest_p"] = float(w.pvalue)
                        nn_p = float(w.pvalue)
                    except ValueError:
                        row["nearest_p"] = np.nan
                else:
                    row["nearest_p"] = np.nan
            else:
                row["nearest_median_delta"] = np.nan
                row["nearest_mean_delta"] = np.nan
                row["nearest_p"] = np.nan

        # Prefer regression effect for direction; fall back to nearest delta.
        effect = row.get("regression_beta", np.nan)
        if not np.isfinite(effect):
            effect = row.get("nearest_median_delta", np.nan)
        pref_p = reg_p if np.isfinite(reg_p) else nn_p
        if not np.isfinite(effect) or not np.isfinite(pref_p) or pref_p >= 0.05:
            direction = "ns"
        elif effect > 0:
            direction = "transcription-driven"
        elif effect < 0:
            direction = "stabilization-driven"
        else:
            direction = "ns"
        row["direction"] = direction
        row["preferred_p"] = pref_p if np.isfinite(pref_p) else np.nan
        rows.append(row)

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # BH across programs on preferred_p
    p = out["preferred_p"].to_numpy(float)
    finite = np.isfinite(p)
    fdr = np.full(len(p), np.nan)
    if finite.any():
        pf = p[finite]
        order = np.argsort(pf)
        ranked = pf[order]
        m = len(pf)
        fdr_sorted = ranked * m / (np.arange(1, m + 1))
        fdr_sorted = np.minimum.accumulate(fdr_sorted[::-1])[::-1]
        fdr_f = np.empty(m)
        fdr_f[order] = np.clip(fdr_sorted, 0, 1)
        fdr[finite] = fdr_f
    out["fdr"] = fdr
    out["significant"] = out["fdr"] < 0.05
    # re-mark ns when FDR fails (keep direction only when significant)
    sig = out["significant"].fillna(False).to_numpy(bool)
    out.loc[~sig, "direction"] = "ns"
    return out.sort_values("preferred_p", na_position="last").reset_index(drop=True)
