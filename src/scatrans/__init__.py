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
from .enrich import list_bundled_gene_sets, run_enrichment, run_kegg, simplify_enrichment
from .generate_gene_features import main as generate_gene_features_main
from .pp_bias import add_gene_features, list_available_gene_features
from .tl import active_score, diagnose_design, filter_active_genes

__all__ = [
    "active_score",
    "diagnose_design",
    "filter_active_genes",
    "add_gene_features",
    "list_available_gene_features",
    "run_enrichment",
    "run_kegg",
    "simplify_enrichment",
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
    __version__ = "0.7.0.dev0"

# Optional: prevent some internal modules from appearing too prominently
# in casual inspection while still allowing advanced users to do
# `import scatrans.tl as tl` if they really need it.
# We do not delete them aggressively to preserve scanpy-like ergonomics.
