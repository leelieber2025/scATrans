"""Internal helpers and shared constants for ``scatrans.tl``."""

from __future__ import annotations

import logging
import warnings
from typing import Any

import anndata as ad
import pandas as pd

# qc is imported lazily inside active_score to keep startup light, but exposed at package level

try:
    from .. import _version

    VERSION = _version.version
except (ImportError, AttributeError):
    VERSION = "0.0.0+unknown"

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

_NOT_PROVIDED = object()

# Shared defaults for filter_active_genes(preset="heuristic") and the built-in
# active_score() significant list (single-cell scale). Keep in sync — both entry
# points should agree under default parameters.
HEURISTIC_FILTER_DEFAULTS: dict[str, float | None] = {
    "active_score_cutoff": 55.0,
    "pval_cutoff": 0.05,
    "unspliced_excess_residual_cutoff": 1.0,
    "logfc_cutoff": 0.35,
    "active_score_fdr_cutoff": 0.25,
    "unspliced_excess_fdr_cutoff": 0.05,
    # Gamma bounds only when user sets them — reference U/S ratios often lie
    # outside (0.05, 1.0); applying defaults after show_effective_gamma=True
    # would silently drop most genes.
    "effective_gamma_min": None,
    "effective_gamma_max": None,
    "delta_variance_min": None,
}

# Pseudobulk scale: residuals / composite scores are much smaller after sample
# aggregation. Used by filter_active_genes(preset="pseudobulk") and the built-in
# significant list when is_pseudobulk=True.
PSEUDOBULK_FILTER_DEFAULTS: dict[str, float | None] = {
    "active_score_cutoff": 5.0,
    "pval_cutoff": 0.05,
    "unspliced_excess_residual_cutoff": 0.05,
    "logfc_cutoff": 0.2,
    "active_score_fdr_cutoff": 0.25,
    "unspliced_excess_fdr_cutoff": 0.05,
    "effective_gamma_min": None,
    "effective_gamma_max": None,
    "delta_variance_min": None,
}

# MixedLM needs enough biological replicates for identifiable random effects.
MIXED_MODEL_MIN_SAMPLES_PER_GROUP = 4
MIXED_MODEL_MIN_TOTAL_SAMPLES = 6

# Keep in sync with _permutation.run_permutation_test (n_success threshold for FDR).
_PERM_FDR_MIN_SUCCESS = 100


def _materialize_if_view(adata: ad.AnnData) -> ad.AnnData:
    """Return a writable AnnData when ``adata`` is an obs/var view (needed before in-place ``.obs`` writes)."""
    if getattr(adata, "is_view", False):
        return adata.copy()
    return adata


def _select_obs(
    adata: ad.AnnData,
    mask: pd.Series,
    *,
    copy_input: bool,
) -> ad.AnnData:
    """Subset to obs ``mask``; call ``AnnData.copy()`` only when ``copy_input=True``."""
    mask = mask.reindex(adata.obs_names, fill_value=False)
    if copy_input:
        return adata[mask].copy()
    if bool(mask.all()):
        return adata
    return adata[mask]


def _select_var(
    adata: ad.AnnData,
    mask: pd.Series,
    *,
    copy_input: bool,
) -> ad.AnnData:
    """Subset to var ``mask``; call ``AnnData.copy()`` only when ``copy_input=True``."""
    if copy_input:
        return adata[:, mask].copy()
    if bool(mask.all()):
        return adata
    return adata[:, mask]


def _resolve_velocity_layer_keys(adata: ad.AnnData) -> tuple[str, str] | None:
    """Return (spliced_key, unspliced_key) for QC, including kb_python mature/nascent."""
    layers = set(adata.layers.keys())
    if "spliced" in layers and "unspliced" in layers:
        return "spliced", "unspliced"
    if "mature" in layers and "nascent" in layers:
        return "mature", "nascent"
    return None


def _de_first_sort_keys(columns: Any) -> tuple[list[str], list[bool]]:
    """Sort spec for DE-first ordering: ``p_adj`` ascending, then ``logFC`` descending.

    Each key is paired with its own direction, so a missing ``p_adj`` cannot flip
    ``logFC`` to ascending. Falls back to the legacy composite ``active_score``
    (descending) only when neither DE column is present. Returns ``([], [])`` if
    none of the keys exist (caller should then leave the order untouched).
    """
    spec = [(c, asc) for c, asc in (("p_adj", True), ("logFC", False)) if c in columns]
    if not spec and "active_score" in columns:
        spec = [("active_score", False)]
    return [c for c, _ in spec], [asc for _, asc in spec]


def _validate_de_common_options(
    *,
    de_preprocess: str,
    pseudobulk_de_backend: str,
    n_jobs: int,
    use_permutation: bool,
    n_perm: int,
    use_mixed_model: bool,
    mixed_model_pval: str,
    paired_replicates: bool,
    use_memento_de: bool,
    memento_capture_rate: float,
    memento_num_boot: int,
    min_cells: int | None = None,
    min_counts: int | None = None,
) -> None:
    """Shared early validation for DE-related options used by both
    active_score() and differential_expression()."""
    if de_preprocess not in {"auto", "normalize_log1p", "none"}:
        raise ValueError("de_preprocess must be one of {'auto', 'normalize_log1p', 'none'}.")

    if pseudobulk_de_backend not in {"pydeseq2", "scanpy"}:
        raise ValueError("pseudobulk_de_backend must be 'pydeseq2' or 'scanpy'.")

    if use_mixed_model and mixed_model_pval not in ("wald", "lrt"):
        raise ValueError("mixed_model_pval must be 'wald' or 'lrt'.")

    if use_memento_de:
        if not (0 < memento_capture_rate < 1):
            raise ValueError(
                "memento_capture_rate must be in (0, 1). "
                "Typical values: ~0.07 for 10x v1, ~0.15 for v2."
            )

        if memento_num_boot < 100:
            raise ValueError(
                "memento_num_boot should be reasonably large (>=100) for stable estimates."
            )

    if min_cells is not None and min_cells < 1:
        raise ValueError("min_cells must be >= 1.")

    if min_counts is not None and min_counts < 0:
        raise ValueError("min_counts must be non-negative.")

    # Extra type guards to prevent bool/int confusion and weird values
    if not isinstance(use_permutation, bool):
        raise ValueError("use_permutation must be boolean.")

    if not isinstance(use_mixed_model, bool):
        raise ValueError("use_mixed_model must be boolean.")

    if not isinstance(paired_replicates, bool):
        raise ValueError("paired_replicates must be boolean.")

    if not isinstance(use_memento_de, bool):
        raise ValueError("use_memento_de must be boolean.")

    if use_memento_de and (
        not isinstance(memento_num_boot, int) or isinstance(memento_num_boot, bool)
    ):
        raise ValueError("memento_num_boot must be an integer.")

    if not isinstance(n_perm, int) or isinstance(n_perm, bool):
        raise ValueError("n_perm must be an integer.")

    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool):
        raise ValueError("n_jobs must be an integer.")
    if n_jobs < -1 or n_jobs == 0:
        raise ValueError(f"n_jobs must be -1 (all CPUs) or a positive integer; got {n_jobs}.")

    if min_cells is not None and (not isinstance(min_cells, int) or isinstance(min_cells, bool)):
        raise ValueError("min_cells must be an integer >= 1.")

    if min_counts is not None and (not isinstance(min_counts, int) or isinstance(min_counts, bool)):
        raise ValueError("min_counts must be a non-negative integer.")


def _coerce_memento_de_preprocess(use_memento_de: bool, de_preprocess: str) -> str:
    """If use_memento_de, force de_preprocess to 'none' (Memento needs raw counts)."""
    if use_memento_de and de_preprocess != "none":
        logger.info(
            "use_memento_de=True: forcing de_preprocess='none' "
            "(Memento method-of-moments works on raw counts)."
        )
        return "none"
    return de_preprocess


def _require_explicit_groups(
    target_group: str | None,
    reference_group: str | None,
    *,
    func_name: str,
) -> None:
    """Require explicit contrast labels; prevents silent GA/Ctrl mismatches."""
    if target_group is None or reference_group is None:
        raise ValueError(
            f"{func_name}() requires explicit target_group and reference_group "
            "matching adata.obs[groupby] values. "
            'Historical defaults "GA"/"Ctrl" were removed to prevent silent mismatches '
            "that yield empty subsets and NaN results. "
            "Use active_score_simple() / differential_expression_simple() for "
            "Disease/Control convenience defaults, or call recommend_workflow() first."
        )


def _resolve_deprecated_active_score_kwargs(kwargs: dict[str, Any]) -> tuple[float, bool]:
    """Pop deprecated active_score kwargs and emit DeprecationWarning when used."""
    active_fdr_cutoff = 0.05
    prioritize_velocity = False
    if "active_fdr_cutoff" in kwargs:
        active_fdr_cutoff = kwargs.pop("active_fdr_cutoff")
        warnings.warn(
            "active_fdr_cutoff is deprecated and no longer used for the built-in "
            "'significant' gene list. Use unspliced_excess_fdr_cutoff instead.",
            DeprecationWarning,
            stacklevel=3,
        )
        if not (0 < active_fdr_cutoff <= 1):
            raise ValueError("active_fdr_cutoff must be in (0, 1].")
    if "prioritize_velocity" in kwargs:
        prioritize_velocity = kwargs.pop("prioritize_velocity")
        warnings.warn(
            "prioritize_velocity is deprecated; use ranking_mode='nascent_excess' instead.",
            DeprecationWarning,
            stacklevel=3,
        )
    if kwargs:
        raise TypeError(
            f"active_score() got unexpected keyword argument(s): {', '.join(sorted(kwargs))}"
        )
    return active_fdr_cutoff, prioritize_velocity
