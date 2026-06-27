import logging

from scipy import sparse

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
        logger.warning("⚠️ Total counts in spliced and unspliced layers are zero.")
        return 0.0

    unspliced_ratio = sum_u / total

    logger.info(f"📊 Global Unspliced Fraction: {unspliced_ratio:.2%}")

    if unspliced_ratio > warn_threshold:
        logger.warning(
            f"⚠️ WARNING: The overall unspliced fraction ({unspliced_ratio:.2%}) is very high (> {warn_threshold:.0%})."
        )
        logger.warning(
            "   This may indicate technical issues such as nuclear RNA enrichment or genomic DNA contamination."
        )

    return float(unspliced_ratio)
