"""
scATrans public API.

Recommended usage:
    import scatrans as scat
    scat.active_score(...)
    scat.add_gene_features(...)
    scat.generate_gene_features_from_gtf(...)   # if you need custom tables from GTF
    scat.pl.set_style()
    scat.run_enrichment(...)

Submodules `pl` and `qc` are intentionally exposed (scanpy-style convention).
Other internal modules are not part of the stable public surface.
"""

from __future__ import annotations

from . import pl, qc
from ._version import __version__
from .enrich import (
    compare_enrichment,
    concat_compare_results,
    expand_enrichment_genes,
    extract_gene_lists,
    list_bundled_gene_sets,
    run_enrichment,
    run_go,
    run_gsea,
    run_kegg,
    save_enrichment_report,
    simplify_enrichment,
)
from .pp_bias import (
    add_gene_features,
    generate_gene_features_from_gtf,
    list_available_gene_features,
)
from .tl import (
    WORKFLOW_PRESETS,
    PartitionResult,
    PipelineResult,
    active_score,
    active_score_simple,
    adaptive_active_score,
    adaptive_weight,
    add_abundance_normalized_residual,
    add_adaptive_score,
    annotate_mechanism_class,
    diagnose_design,
    differential_expression,
    differential_expression_simple,
    ensure_raw_counts,
    filter_active_genes,
    labeling_anchor,
    partition_de_by_mechanism,
    program_mechanism,
    recommend_workflow,
    restore_raw_counts,
    run_default_pipeline,
    store_raw_counts,
    threshold_sensitivity,
)

__all__ = [
    "active_score",
    "adaptive_active_score",
    "add_adaptive_score",
    "adaptive_weight",
    "labeling_anchor",
    "add_abundance_normalized_residual",
    "annotate_mechanism_class",
    "threshold_sensitivity",
    "program_mechanism",
    "partition_de_by_mechanism",
    "PartitionResult",
    "active_score_simple",
    "differential_expression",
    "differential_expression_simple",
    "diagnose_design",
    "recommend_workflow",
    "WORKFLOW_PRESETS",
    "run_default_pipeline",
    "PipelineResult",
    "filter_active_genes",
    "ensure_raw_counts",
    "restore_raw_counts",
    "store_raw_counts",
    "add_gene_features",
    "generate_gene_features_from_gtf",
    "list_available_gene_features",
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
    "pl",
    "qc",
    "__version__",
]
