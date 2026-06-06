from .tl import active_score
from .pp_bias import add_gene_features
from .enrich import run_enrichment, run_kegg, simplify_enrichment
from . import pl
from . import qc
from .generate_gene_features import main as generate_gene_features_main

__all__ = [
    "active_score",
    "add_gene_features",
    "run_enrichment",
    "run_kegg",
    "simplify_enrichment",
    "pl",
    "qc",
    "generate_gene_features_main",
]

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.7.0-dev"
