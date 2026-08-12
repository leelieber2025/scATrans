"""scATrans enrichment (``enrich``) package.

ORA / GSEA / GO / KEGG helpers with bundled gene-set support.
Import from ``scatrans.enrich`` (or ``scatrans``) for the public surface.

Private helpers live in submodules (e.g. ``scatrans.enrich._data``) and are
not part of ``__all__``.
"""

from __future__ import annotations

import logging

from ._data import (
    BUNDLED_GENE_SET_PROVENANCE,
    GSEA_COLUMNS,
    ORA_COLUMNS,
    list_bundled_gene_sets,
)
from .compare import (
    compare_enrichment,
    concat_compare_results,
    extract_gene_lists,
)
from .gsea import run_gsea
from .ora import run_enrichment, run_go, run_kegg
from .report import expand_enrichment_genes, save_enrichment_report
from .simplify import simplify_enrichment

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

__all__ = [
    "run_enrichment",
    "run_kegg",
    "run_go",
    "run_gsea",
    "simplify_enrichment",
    "save_enrichment_report",
    "expand_enrichment_genes",
    "list_bundled_gene_sets",
    "compare_enrichment",
    "extract_gene_lists",
    "concat_compare_results",
    "BUNDLED_GENE_SET_PROVENANCE",
    "ORA_COLUMNS",
    "GSEA_COLUMNS",
]


def __dir__():
    dunders = {n for n in globals() if n.startswith("__") and n.endswith("__")}
    return sorted(set(__all__) | dunders)
