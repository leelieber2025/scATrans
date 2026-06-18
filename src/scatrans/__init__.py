"""
scATrans public API.

Recommended usage:
    import scatrans as scat
    scat.active_score(...)
    scat.add_gene_features(...)
    scat.pl.set_style()
    scat.run_enrichment(...)

Submodules `pl` and `qc` are intentionally exposed (scanpy-style convention).
Other internal modules are not part of the stable public surface.
"""

from . import pl, qc
from .enrich import (
    expand_enrichment_genes,
    list_bundled_gene_sets,
    run_enrichment,
    run_go,
    run_kegg,
    save_enrichment_report,
    simplify_enrichment,
)
from .generate_gene_features import main as generate_gene_features_main
from .pp_bias import add_gene_features, list_available_gene_features
from .tl import (
    active_score,
    diagnose_design,
    differential_expression,
    ensure_raw_counts,
    filter_active_genes,
    restore_raw_counts,
    store_raw_counts,
)

__all__ = [
    "active_score",
    "differential_expression",
    "diagnose_design",
    "filter_active_genes",
    "ensure_raw_counts",
    "restore_raw_counts",
    "store_raw_counts",
    "add_gene_features",
    "list_available_gene_features",
    "run_enrichment",
    "run_kegg",
    "run_go",
    "simplify_enrichment",
    "save_enrichment_report",
    "expand_enrichment_genes",
    "list_bundled_gene_sets",
    "pl",
    "qc",
    "generate_gene_features_main",
    "__version__",
]

# Version is provided dynamically when possible
try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.8.0"

# Optional: prevent some internal modules from appearing too prominently
# in casual inspection while still allowing advanced users to do
# `import scatrans.tl as tl` if they really need it.
# We do not delete them aggressively to preserve scanpy-like ergonomics.
