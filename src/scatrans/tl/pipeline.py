"""scatrans.tl.pipeline — internal package module."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any, ClassVar

import anndata as ad
import numpy as np
import pandas as pd

from .._utils import (
    _normalize_group_label,
)
from ._common import (
    _PERM_FDR_MIN_SUCCESS,
    VERSION,
)
from .active import active_score
from .de import differential_expression
from .filter import filter_active_genes

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _usable_gene_length_mask(series: Any) -> pd.Series:
    """True where gene_length is finite and > 0 (Huber valid_feat requirement)."""
    gl = pd.to_numeric(series, errors="coerce")
    return gl.notna() & (gl > 0)


def _fill_missing_gene_features_from_bundle(adata: Any, organism: str) -> None:
    """Fill only missing / non-positive gene_length (and missing intron) from bundled tables.

    Never overwrites existing usable (length > 0) values — safe for partial GTF
    feature tables that would otherwise block full auto-attach.
    """
    from ..pp_bias import add_gene_features

    # Snapshot usable lengths before full attach (add_gene_features rewrites columns).
    old_gl = (
        pd.to_numeric(adata.var["gene_length"], errors="coerce")
        if "gene_length" in adata.var.columns
        else pd.Series(np.nan, index=adata.var_names)
    )
    old_intr = (
        pd.to_numeric(adata.var["intron_number"], errors="coerce")
        if "intron_number" in adata.var.columns
        else pd.Series(np.nan, index=adata.var_names)
    )
    usable = _usable_gene_length_mask(old_gl)

    add_gene_features(adata, organism=organism)

    new_gl = pd.to_numeric(adata.var["gene_length"], errors="coerce")
    new_intr = pd.to_numeric(adata.var["intron_number"], errors="coerce")
    # Restore user/GTF positive lengths; fill elsewhere from bundle.
    merged_gl = new_gl.copy()
    merged_gl.loc[usable] = old_gl.loc[usable]
    adata.var["gene_length"] = merged_gl

    # Keep prior intron where present; fill NaN from bundle.
    merged_intr = new_intr.copy()
    prior_intr = old_intr.notna()
    merged_intr.loc[prior_intr] = old_intr.loc[prior_intr]
    adata.var["intron_number"] = merged_intr


def _maybe_add_gene_features(adata: Any, organism: str) -> Any:
    """Attach / complete bundled gene features for bias correction.

    Semantics:
    - If either column is absent, or no gene has *usable* length (finite, > 0)
      and a finite intron → full ``add_gene_features`` attach.
    - If usable length coverage is low (< 50% of genes) → fill missing/zero
      lengths from the bundled table **without overwriting** existing length > 0
      (avoids the "1% real, 99% NaN/0 never auto-completes" trap with partial GTF tables).
    - If coverage is adequate but some genes lack length → log only; active_score
      excludes them via ``valid_feat`` (``gene_length > 0``).
    """
    has_length = "gene_length" in adata.var.columns
    has_intron = "intron_number" in adata.var.columns
    if not has_length or not has_intron:
        from ..pp_bias import add_gene_features

        add_gene_features(adata, organism=organism)
        logger.info("Attached bundled gene features (organism=%s) for bias correction.", organism)
        return adata

    gl = pd.to_numeric(adata.var["gene_length"], errors="coerce")
    intr = pd.to_numeric(adata.var["intron_number"], errors="coerce")
    usable = _usable_gene_length_mask(gl)
    n_genes = max(int(adata.n_vars), 1)
    n_usable = int(usable.sum())
    frac = n_usable / n_genes

    if n_usable == 0 or not intr.notna().any():
        from ..pp_bias import add_gene_features

        add_gene_features(adata, organism=organism)
        logger.info(
            "Attached bundled gene features (organism=%s): no usable gene_length > 0 "
            "was present on the object.",
            organism,
        )
        return adata

    if frac < 0.5:
        logger.warning(
            "gene features: only %.1f%% of genes have usable gene_length > 0 "
            "(%d/%d). Filling missing/non-positive lengths from bundled features "
            "(organism=%s) without overwriting existing positive lengths. "
            "Partial GTF tables alone leave most genes out of Huber fit.",
            100.0 * frac,
            n_usable,
            n_genes,
            organism,
        )
        _fill_missing_gene_features_from_bundle(adata, organism)
        return adata

    n_bad = n_genes - n_usable
    if n_bad > 0:
        logger.info(
            "gene features: %d/%d genes lack usable gene_length > 0; "
            "they are excluded from Huber bias correction (valid_feat).",
            n_bad,
            n_genes,
        )
    return adata


def _resolve_simple_backend_kwargs(
    adata: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None,
) -> dict[str, Any]:
    """Pick pseudobulk vs single-cell defaults from replicate structure."""
    kwargs: dict[str, Any] = {
        "de_method": "wilcoxon",
        "use_pseudobulk": False,
        "sample_col": None,
        "pseudobulk_de_backend": "pydeseq2",
    }
    if sample_col and sample_col in adata.obs.columns:
        norm_groups = adata.obs[groupby].map(_normalize_group_label)
        t_mask = norm_groups == _normalize_group_label(target_group)
        r_mask = norm_groups == _normalize_group_label(reference_group)
        n_s_t = int(adata.obs.loc[t_mask, sample_col].nunique())
        n_s_r = int(adata.obs.loc[r_mask, sample_col].nunique())
        if min(n_s_t, n_s_r) >= 3:
            kwargs.update(
                {
                    "use_pseudobulk": True,
                    "sample_col": sample_col,
                    "de_method": "t-test_overestim_var",
                }
            )
            logger.info(
                "Simple path: detected >=3 samples per group — using pseudobulk + PyDESeq2."
            )
        else:
            logger.info(
                "Simple path: few samples per group (target=%d, reference=%d) — "
                "using single-cell Wilcoxon DE.",
                n_s_t,
                n_s_r,
            )
    return kwargs


def active_score_simple(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",  # convenience defaults (distinct from core active_score)
    sample_col: str | None = None,
    organism: str = "mouse",
    *,
    show_plot: bool = False,
    copy_input: bool = True,
    pydeseq2_min_counts: int = 10,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """
    Recommended entry point for new users (minimal parameters).

    Wraps :func:`active_score` with sensible defaults:
    - Uses "Disease"/"Control" as group defaults (unlike core active_score which defaults
      to the historical "GA"/"Ctrl").
    - heuristic mode, no permutation (inspect ``all_results`` + ``filter_active_genes``)
    - auto-attaches bundled gene features when missing
    - pseudobulk + PyDESeq2 when ``sample_col`` has >=3 replicates per group;
      otherwise single-cell Wilcoxon DE

    For full control (permutation, advanced mode, mixed models, etc.) use
    :func:`active_score` directly.

    copy_input : bool, default True
        If True, isolate a working copy before attaching gene features or DE.
        If False, still isolates before mutation so the caller's AnnData is not
        modified (same contract as :func:`active_score`).
    """
    # Always isolate before _maybe_add_gene_features (writes .var) so copy_input=False
    # does not pollute the caller's feature columns.
    adata = adata_input.copy()
    if not copy_input:
        logger.info(
            "active_score_simple(copy_input=False): still isolating a working copy "
            "before attaching gene features / DE so the caller's AnnData is unchanged."
        )
    _maybe_add_gene_features(adata, organism)
    backend = _resolve_simple_backend_kwargs(
        adata, groupby, target_group, reference_group, sample_col
    )
    return active_score(
        adata,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        mode="heuristic",
        use_permutation=False,
        show_plot=show_plot,
        copy_input=False,
        pydeseq2_min_counts=pydeseq2_min_counts,
        **backend,
    )


def differential_expression_simple(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",
    sample_col: str | None = None,
    *,
    copy_input: bool = True,
    pydeseq2_min_counts: int = 10,
) -> tuple[ad.AnnData, pd.DataFrame]:
    """
    Minimal-parameter differential expression (no velocity layers required).

    Same backend auto-selection as :func:`active_score_simple`.
    For Memento, mixed models, or custom preprocess use :func:`differential_expression`.
    """
    # Isolate so DE never mutates the caller's object (labels / preprocess).
    adata = adata_input.copy()
    if not copy_input:
        logger.info(
            "differential_expression_simple(copy_input=False): still isolating a "
            "working copy so the caller's AnnData is unchanged."
        )
    backend = _resolve_simple_backend_kwargs(
        adata, groupby, target_group, reference_group, sample_col
    )
    return differential_expression(
        adata,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        copy_input=False,
        pydeseq2_min_counts=pydeseq2_min_counts,
        **backend,
    )


_PIPELINE_READONLY_MSG = "PipelineResult is read-only; use result.to_dict() for a mutable copy"


class PipelineResult(dict):
    """Structured result of :func:`run_default_pipeline`.

    A **read-only** :class:`dict` subclass so legacy code keeps working:

    - ``isinstance(result, dict)`` is ``True``
    - ``result["candidates"]``, ``for k in result``, ``dict(result)``, ``{**result}``
    - attribute access: ``result.candidates`` (same keys)

    Mutation is blocked, including C-level dict paths that bypass a plain
    Python ``__setitem__`` override:

    - item / attribute assignment, ``update`` / ``pop`` / ``clear`` / ``setdefault``
    - in-place merge ``result |= other`` (``__ior__``)
    - ``result.copy()`` returns a **mutable** plain ``dict`` (same as
      :meth:`to_dict`); ``copy.copy(result)`` returns another read-only
      :class:`PipelineResult`
    - ``result | other`` returns a **new mutable** ``dict`` (does not mutate
      ``result`` and does not silently drop attribute access — callers that need
      ``.candidates`` should keep the original or rebuild a ``PipelineResult``)

    Nested values such as ``backend`` / ``meta`` dicts are still plain dicts.
    ``copy.deepcopy`` / pickle reconstruct via :meth:`__init__`.
    """

    _KEYS: ClassVar[tuple[str, ...]] = (
        "adata",
        "significant",
        "all_results",
        "candidates",
        "enrichment",
        "filter_preset",
        "backend",
        "meta",
    )

    __slots__ = ("_readonly",)

    def __init__(
        self,
        adata: Any,
        significant: pd.DataFrame,
        all_results: pd.DataFrame,
        candidates: pd.DataFrame,
        enrichment: pd.DataFrame | None,
        filter_preset: str,
        backend: dict[str, Any],
        meta: dict[str, Any] | None = None,
    ) -> None:
        object.__setattr__(self, "_readonly", False)
        super().__init__(
            adata=adata,
            significant=significant,
            all_results=all_results,
            candidates=candidates,
            enrichment=enrichment,
            filter_preset=filter_preset,
            backend=backend,
            meta={} if meta is None else meta,
        )
        object.__setattr__(self, "_readonly", True)

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as exc:
            raise AttributeError(f"{type(self).__name__!s} has no attribute {name!r}") from exc

    def __setattr__(self, name: str, value: Any) -> None:
        if name == "_readonly":
            object.__setattr__(self, name, value)
            return
        if getattr(self, "_readonly", False):
            raise TypeError(_PIPELINE_READONLY_MSG)
        object.__setattr__(self, name, value)

    def __setitem__(self, key: str, value: Any) -> None:  # type: ignore[override]
        if getattr(self, "_readonly", False):
            raise TypeError(_PIPELINE_READONLY_MSG)
        super().__setitem__(key, value)

    def __delitem__(self, key: str) -> None:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def clear(self) -> None:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def pop(self, *args: Any, **kwargs: Any) -> Any:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def popitem(self) -> Any:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def setdefault(self, *args: Any, **kwargs: Any) -> Any:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def update(self, *args: Any, **kwargs: Any) -> None:  # type: ignore[override]
        raise TypeError(_PIPELINE_READONLY_MSG)

    def __ior__(self, other: Any) -> PipelineResult:  # type: ignore[override,misc]
        # dict.__ior__ is C-level and would otherwise mutate in place silently.
        raise TypeError(_PIPELINE_READONLY_MSG)

    def __or__(self, other: Any) -> dict[str, Any]:  # type: ignore[override]
        """Return a **new mutable** ``dict`` (``self`` is unchanged)."""
        if not isinstance(other, Mapping):
            return NotImplemented
        return {**self, **dict(other)}

    def __ror__(self, other: Any) -> dict[str, Any]:  # type: ignore[override]
        if not isinstance(other, Mapping):
            return NotImplemented
        return {**dict(other), **self}

    def copy(self) -> dict[str, Any]:  # type: ignore[override]
        """Mutable shallow :class:`dict` copy (alias of :meth:`to_dict`).

        Note: this differs from ``copy.copy(self)``, which returns another
        read-only :class:`PipelineResult` via :meth:`__copy__`.
        """
        return self.to_dict()

    def __copy__(self) -> PipelineResult:
        """Shallow copy that stays a read-only :class:`PipelineResult`."""
        return PipelineResult(**{k: self[k] for k in self._KEYS})

    def __deepcopy__(self, memo: dict[int, Any]) -> PipelineResult:
        """Deep copy via ``__init__`` (avoids frozen ``__setitem__`` during rebuild)."""
        import copy as _copy

        if id(self) in memo:
            return memo[id(self)]
        kwargs = {k: _copy.deepcopy(self[k], memo) for k in self._KEYS}
        new = PipelineResult(**kwargs)
        memo[id(self)] = new
        return new

    def __reduce__(self) -> tuple[type, tuple[Any, ...]]:
        """Pickle via constructor args (all protocols / joblib).

        ``__reduce__`` is the sole serialization path: once defined, pickle
        never calls ``__getstate__`` / ``__setstate__``, so those hooks are
        intentionally absent (``copy.deepcopy`` uses ``__deepcopy__``).
        """
        return (type(self), tuple(self[k] for k in self._KEYS))

    def to_dict(self) -> dict[str, Any]:
        """Shallow mutable :class:`dict` copy of the result fields."""
        return dict(self)

    def summary(self) -> dict[str, Any]:
        """Compact counts for logging and quick inspection."""
        enrichment = self["enrichment"]
        n_enrich = None if enrichment is None else int(len(enrichment))
        return {
            "n_significant": int(len(self["significant"])),
            "n_candidates": int(len(self["candidates"])),
            "n_all_results": int(len(self["all_results"])),
            "n_enrichment_terms": n_enrich,
            "filter_preset": self["filter_preset"],
            "backend": dict(self["backend"]),
        }


def run_default_pipeline(
    adata_input: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",  # convenience defaults; always prefer explicit values
    sample_col: str | None = None,
    organism: str = "mouse",
    *,
    run_go_enrichment: bool = True,
    gene_sets: str = "GO_Biological_Process",
    filter_preset: str | None = None,
    show_plot: bool = False,
) -> PipelineResult:
    """
    End-to-end recommended workflow for first-time users.

    Steps: active scoring → ``filter_active_genes`` → optional GO enrichment.
    Uses "Disease"/"Control" convenience defaults for target/reference.

    The default for ``filter_preset`` is now **auto-detected** from the experimental
    design (via ``_resolve_simple_backend_kwargs``): "pseudobulk" when sample_col
    is provided with >=3 samples per group (which triggers pseudobulk inside
    active_score_simple), otherwise "heuristic". This keeps the thresholds
    consistent with the actual scale of active_score / unspliced_excess_residual
    (see WORKFLOW_PRESETS["pseudobulk_report"]).

    Returns a :class:`PipelineResult` (mapping-compatible) with fields:
      - ``adata``, ``significant``, ``all_results``, ``candidates``
      - ``enrichment`` (DataFrame or None)
      - ``filter_preset``, ``backend`` (kwargs used for DE)
      - ``meta``: always includes ``scatrans_version`` and ``organism``; when
        ``active_score`` wrote ``adata.uns["scatrans"]``, also surfaces its
        nested ``diagnostics`` block plus a few high-value run flags
        (``use_permutation``, ``gamma_method``, ``mode``,
        ``unspliced_global_fraction``, …). Full run metadata remains on
        ``result.adata.uns["scatrans"]``.
    """
    # Resolve once so we can pick a matching filter_preset (addresses mismatch
    # between auto-pseudobulk in active_score_simple and hardcoded "heuristic").
    backend = _resolve_simple_backend_kwargs(
        adata_input, groupby, target_group, reference_group, sample_col
    )
    if filter_preset is None:
        filter_preset = "pseudobulk" if backend.get("use_pseudobulk") else "heuristic"

    adata_res, significant, all_results = active_score_simple(
        adata_input,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        organism=organism,
        show_plot=show_plot,
    )
    candidates = filter_active_genes(all_results, preset=filter_preset)

    # active_score writes rich run metadata under .uns["scatrans"] (diagnostics,
    # gamma_method, use_permutation, …). Keep this separate from the result
    # ``meta`` dict so the two are not confused.
    scatrans_uns = adata_res.uns.get("scatrans", {})
    if not isinstance(scatrans_uns, dict):
        scatrans_uns = {}

    if len(significant) == 0 and not scatrans_uns.get("use_permutation"):
        logger.info(
            "Built-in significant is empty (use_permutation=False). "
            "Candidate genes come from filter_active_genes(preset=%r). "
            "For permutation-based lists, re-run with use_permutation=True only when "
            "the design supports >=%d successful shuffles.",
            filter_preset,
            _PERM_FDR_MIN_SUCCESS,
        )
    elif (
        len(significant) > 0
        and filter_preset == "significant"
        and not candidates.index.equals(significant.index)
    ):
        logger.warning(
            "filter_active_genes(preset='significant') did not exactly match the built-in "
            "significant list; prefer the significant return value directly."
        )

    enrichment = None
    if run_go_enrichment and len(candidates) > 0:
        # Lazy import: enrichment pulls gene-set I/O / optional gseapy paths;
        # not required when run_go_enrichment=False (no circular import risk).
        from ..enrich import run_enrichment

        enrichment = run_enrichment(
            candidates.index.tolist(),
            gene_sets=gene_sets,
            organism=organism,
            adata=adata_res,
            pval_cutoff=0.05,
        )

    # Package-level result meta: version + organism always; fold in diagnostics
    # and selected run flags from active_score when present.
    result_meta: dict[str, Any] = {
        "scatrans_version": VERSION,
        "organism": organism,
    }
    diag = scatrans_uns.get("diagnostics")
    if diag is not None:
        result_meta["diagnostics"] = diag
    for key in (
        "use_permutation",
        "gamma_method",
        "mode",
        "velocity_source",
        "unspliced_global_fraction",
        "de_method",
        "use_pseudobulk",
        "sample_col",
        "ranking_mode",
    ):
        if key in scatrans_uns:
            result_meta[key] = scatrans_uns[key]

    # We already resolved backend above (avoids calling the resolver a second time
    # just to "guess" what active_score_simple decided internally).
    return PipelineResult(
        adata=adata_res,
        significant=significant,
        all_results=all_results,
        candidates=candidates,
        enrichment=enrichment,
        filter_preset=filter_preset,
        backend=backend,
        meta=result_meta,
    )
