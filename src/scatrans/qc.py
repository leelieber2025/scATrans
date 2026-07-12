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


__all__ = ["unspliced_global"]


def __dir__():
    return sorted(__all__)
