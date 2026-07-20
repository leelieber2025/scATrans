from __future__ import annotations

import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def unspliced_global(adata, spliced_key="spliced", unspliced_key="unspliced", warn_threshold=0.5):
    """
    Calculate the global fraction of unspliced reads across all cells.
    Generates a warning if the fraction is abnormally high (e.g., > 50-60%),
    which could indicate nuclear enrichment or genomic DNA contamination.
    """
    if spliced_key not in adata.layers or unspliced_key not in adata.layers:
        raise ValueError(f"Layers '{spliced_key}' and/or '{unspliced_key}' not found in adata.")

    S = adata.layers[spliced_key]
    U = adata.layers[unspliced_key]

    # Both numpy arrays and scipy sparse matrices support .sum(); no need for the branch.
    sum_s = float(S.sum())
    sum_u = float(U.sum())

    total = sum_s + sum_u
    if total == 0:
        logger.warning("Total counts in spliced and unspliced layers are zero.")
        return 0.0

    unspliced_ratio = sum_u / total

    logger.info("Global Unspliced Fraction: %.2f%%", unspliced_ratio * 100)

    if unspliced_ratio > warn_threshold:
        logger.warning(
            "WARNING: The overall unspliced fraction (%.2f%%) is very high (> %.0f%%).",
            unspliced_ratio * 100,
            warn_threshold * 100,
        )
        logger.warning(
            "   This may indicate technical issues such as nuclear RNA enrichment or genomic DNA contamination."
        )

    return float(unspliced_ratio)


def _reliability_from_unspliced_fraction(
    f: float,
    *,
    low_floor: float = 0.02,
    low_ok: float = 0.10,
    high_ok: float = 0.45,
    high_ceil: float = 0.70,
) -> float:
    """Map the global unspliced fraction to a proxy-reliability scalar in [0, 1].

    U-shaped: the nascent proxy / gamma fit is trustworthy in a normal unspliced
    band and degrades at both extremes — too LOW (little nascent signal, noise
    dominated) and too HIGH (nuclear/gDNA contamination -> gamma mis-fit; e.g.
    HTEC 67.7% unspliced where the proxy failed, AUC 0.49). Full reliability on
    ``[low_ok, high_ok]``, linear ramps to 0 at ``low_floor`` / ``high_ceil``.
    """
    if not (low_floor < low_ok <= high_ok < high_ceil):
        raise ValueError("require low_floor < low_ok <= high_ok < high_ceil")
    if f <= low_floor or f >= high_ceil:
        return 0.0
    if low_ok <= f <= high_ok:
        return 1.0
    if f < low_ok:
        return (f - low_floor) / (low_ok - low_floor)
    return (high_ceil - f) / (high_ceil - high_ok)


def regime_diagnosis(
    adata,
    *,
    spliced_key: str = "spliced",
    unspliced_key: str = "unspliced",
    high_ok: float = 0.45,
    high_ceil: float = 0.70,
    low_ok: float = 0.10,
) -> dict:
    """Pre-flight proxy-reliability diagnosis from the global unspliced fraction.

    Returns a dataset-level verdict + a ``reliability`` scalar in [0, 1] suitable
    to pass to :func:`scatrans.tl.annotate_mechanism_class` (``reliability=``) so
    the per-gene mechanism confidence reflects data quality.

    NOTE — scope: this is the DATA-QUALITY / gamma-reliability half of the regime
    check (unspliced-fraction QC). It does NOT distinguish *dynamic vs
    steady-state* transcription, which is what actually governs whether the proxy
    beats DE. The candidate signal for that half — per-timepoint RNA-velocity
    ``velocity_length`` — was cross-validated on sci-fate and REJECTED (it tracked
    true proxy reliability on scNT at +0.80 but not across platforms, -0.50 on
    sci-fate: its absolute magnitude is not comparable across depth/chemistry), so
    no dynamic-vs-steady signal is provided; a validated label-free one is an open
    problem. Treat a high ``reliability`` here as "the proxy is not obviously
    corrupted", not as "the proxy is in its winning regime".

    Returns
    -------
    dict with ``unspliced_fraction``, ``reliability``, ``regime``
    ("ok" / "low_unspliced" / "high_unspliced"), ``basis`` and ``message``.
    """
    f = unspliced_global(
        adata, spliced_key=spliced_key, unspliced_key=unspliced_key, warn_threshold=high_ok
    )
    reliability = _reliability_from_unspliced_fraction(
        f, low_ok=low_ok, high_ok=high_ok, high_ceil=high_ceil
    )
    if f > high_ok:
        regime = "high_unspliced"
        msg = (
            f"unspliced fraction {f:.1%} is high (>= {high_ok:.0%}): possible nuclear/gDNA "
            "contamination -> gamma fit and the nascent proxy may be unreliable; "
            "mechanism annotations down-weighted."
        )
    elif f < low_ok:
        regime = "low_unspliced"
        msg = (
            f"unspliced fraction {f:.1%} is low (< {low_ok:.0%}): little nascent signal, "
            "the proxy is noise-dominated; mechanism annotations down-weighted."
        )
    else:
        regime = "ok"
        msg = f"unspliced fraction {f:.1%} is in the normal band; proxy not obviously corrupted."
    return {
        "unspliced_fraction": float(f),
        "reliability": float(reliability),
        "regime": regime,
        "basis": "unspliced_fraction",
        "message": msg,
    }


__all__ = ["unspliced_global", "regime_diagnosis"]


def __dir__():
    return sorted(__all__)
