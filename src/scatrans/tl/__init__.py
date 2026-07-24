"""scATrans tools (``tl``) package.

Public entry points for active-transcription scoring, DE, filtering, and
guided workflows. Implementation is split across submodules; import from
``scatrans.tl`` (or ``scatrans``) for the stable surface.

Private helpers live in submodules (e.g. ``scatrans.tl._common``) and are
not part of ``__all__``.
"""

from __future__ import annotations

import logging

from ._common import (
    HEURISTIC_FILTER_DEFAULTS,
    MIXED_MODEL_MIN_SAMPLES_PER_GROUP,
    MIXED_MODEL_MIN_TOTAL_SAMPLES,
    PSEUDOBULK_FILTER_DEFAULTS,
)
from .active import active_score
from .adaptive import (
    adaptive_active_score,
    adaptive_weight,
    add_adaptive_score,
    labeling_anchor,
)
from .bias import add_abundance_normalized_residual
from .de import (
    differential_expression,
    ensure_raw_counts,
    restore_raw_counts,
    store_raw_counts,
)
from .design import WORKFLOW_PRESETS, diagnose_design, recommend_workflow
from .filter import filter_active_genes
from .mechanism import (
    annotate_mechanism_class,
    program_mechanism,
    program_mechanism_induction_matched,
    threshold_sensitivity,
)
from .nascent import nascent_activity_score
from .partition import PartitionResult, partition_de_by_mechanism
from .pipeline import (
    PipelineResult,
    active_score_simple,
    differential_expression_simple,
    run_default_pipeline,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

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
    "program_mechanism_induction_matched",
    "nascent_activity_score",
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
    "HEURISTIC_FILTER_DEFAULTS",
    "PSEUDOBULK_FILTER_DEFAULTS",
    "MIXED_MODEL_MIN_SAMPLES_PER_GROUP",
    "MIXED_MODEL_MIN_TOTAL_SAMPLES",
]


def __dir__():
    dunders = {n for n in globals() if n.startswith("__") and n.endswith("__")}
    return sorted(set(__all__) | dunders)
