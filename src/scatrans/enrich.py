"""Enrichment analysis wrapper for scATrans.

Currently a thin convenience wrapper around gseapy (if installed).
Full integration (including background gene list from adata) is planned.
"""

import warnings
from typing import List, Optional, Union, Dict, Any


def run_enrichment(
    gene_list: List[str],
    gene_sets: Union[str, List[str]] = "GO_Biological_Process_2023",
    organism: str = "mouse",
    background: Optional[List[str]] = None,
    cutoff: float = 0.05,
    **kwargs
) -> Any:
    """
    Run functional enrichment analysis on a list of active genes.

    This is a convenience wrapper. For full control, use gseapy directly.

    Parameters
    ----------
    gene_list : list of str
        List of gene symbols (or IDs) that are significant active drivers.
    gene_sets : str or list
        Gene set database, e.g. "GO_Biological_Process_2023", "KEGG_2021_Human", etc.
        See gseapy documentation for available libraries.
    organism : {"mouse", "human"}
        Organism for the gene sets.
    background : list of str, optional
        Background gene universe (recommended: all genes expressed in your dataset).
        If None, gseapy will use its default.
    cutoff : float
        FDR cutoff to return significant terms.
    **kwargs
        Passed to gseapy.enrichr or prerank.

    Returns
    -------
    gseapy object or DataFrame with enrichment results.
    """
    try:
        import gseapy as gp
    except ImportError:
        raise ImportError(
            "run_enrichment requires the 'gseapy' package. "
            "Install with: pip install gseapy"
        ) from None

    if isinstance(gene_sets, str):
        gene_sets = [gene_sets]

    # Simple enrichr call (most common for GO/KEGG)
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism=organism,
        background=background,
        cutoff=cutoff,
        **kwargs
    )
    return enr


# Also expose a simple dotplot if user wants (but recommend using gseapy's plot)
def enrich_dotplot(enrichment_result, title: str = "Enrichment", save_path: Optional[str] = None, **kwargs):
    """Placeholder / thin wrapper. Use gseapy.plot.dotplot directly for best results."""
    try:
        import gseapy as gp
        fig = gp.plot.dotplot(enrichment_result.res2d, title=title, **kwargs)
        if save_path:
            fig.savefig(save_path, bbox_inches="tight", dpi=300)
        return fig
    except Exception as e:
        warnings.warn("enrich_dotplot is a thin wrapper. For full customization use gseapy directly.")
        raise NotImplementedError("Advanced dotplot styling not yet fully in scATrans.pl") from e
