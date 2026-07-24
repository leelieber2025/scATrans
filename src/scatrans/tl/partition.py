"""scatrans.tl.partition — the DE→mechanism primary workflow.

scATrans does not replace differential expression for gene discovery. Instead:

1. a **standard DE** step selects the changed genes (package backends via
   :func:`~scatrans.tl.differential_expression`, or an external table /
   callable);
2. scATrans **partitions** those DE genes by mechanism — the nascent residual
   scores *transcription support*, labeling transcription-driven versus
   stabilization-driven changes (both can be real expression changes).

Design rules:

- per-gene labels are soft, low-confidence hints
  (``transcription_support``, ``mechanism_class``, ``mechanism_confidence``),
  with confidence scaled by a reliability pre-flight
  (:func:`~scatrans.qc.regime_diagnosis`);
- stronger claims should use program-level pooling
  (:func:`~scatrans.tl.program_mechanism`);
- the proxy never filters or removes DE hits — it only annotates and ranks.

Down-regulation is not yet mechanism-resolved (``unclassified_down``).
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from ..qc import regime_diagnosis
from ._common import VERSION
from .de import differential_expression
from .filter import filter_active_genes
from .mechanism import (
    CLASS_COL,
    INDUCTION_CONFOUND_COL,
    annotate_mechanism_class,
    program_mechanism,
    program_mechanism_induction_matched,
)
from .nascent import (
    REPRO_COL,
    nascent_activity_score,
)
from .pipeline import active_score_simple

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# "builtin" reuses the DE that active_score already computed on all_results.
_BUILTIN_DE = "builtin"


@dataclass
class PartitionResult:
    """Result of :func:`partition_de_by_mechanism`.

    Attributes
    ----------
    adata
        The scored AnnData (proxy columns in ``.var`` / run meta in ``.uns``).
    regime
        Reliability pre-flight (:func:`~scatrans.qc.regime_diagnosis`): the
        ``reliability`` scalar here scales every ``mechanism_confidence`` and
        tells you how far to trust the per-gene calls on this dataset.
    gene_table
        FULL scored table (all tested genes) with the mechanism annotation
        columns (``transcription_support`` / ``mechanism_class`` /
        ``mechanism_confidence``) added. Annotation is non-destructive: DE
        membership is NOT decided here.
    selected
        The DE-selected genes (rows of ``gene_table``) — the changed program,
        chosen by the DE step alone.
    programs
        Program-level transcription-vs-stabilization table
        (:func:`~scatrans.tl.program_mechanism`) when ``gene_sets`` was given,
        else ``None``. This is where the decisive calls live (prefer over
        per-gene ``mechanism_class``).
    enrichment
        Optional GO/pathway enrichment on ``selected`` (or ``None``). Do **not**
        enrich gene lists split by per-gene ``mechanism_class`` (induction trap).
    meta
        Run metadata: version, organism, DE source, thresholds, regime, and the
        mechanism-annotation diagnostics.
    programs_induction_matched
        Induction-controlled program table
        (:func:`~scatrans.tl.program_mechanism_induction_matched`) when requested,
        else ``None``. Prefer this for claims when induction strength varies.
        Placed after ``meta`` so existing positional ``PartitionResult``
        constructions keep working.
    """

    adata: Any
    regime: dict[str, Any]
    gene_table: pd.DataFrame
    selected: pd.DataFrame
    programs: pd.DataFrame | None = None
    enrichment: pd.DataFrame | None = None
    meta: dict[str, Any] = field(default_factory=dict)
    programs_induction_matched: pd.DataFrame | None = None

    def summary(self) -> dict[str, Any]:
        """Compact counts for logging / quick inspection.

        Program-level fields are listed first: per-gene ``mechanism_class`` counts
        are soft hints only (do not use them for ORA / hard claims).
        """
        soft_counts = {
            k: int((self.selected.get(CLASS_COL) == k).sum())
            for k in (
                "transcription-driven",
                "stabilization-driven",
                "ambiguous",
                "unclassified_down",
                "unknown",
            )
            if CLASS_COL in self.selected.columns
        }
        n_prog = None if self.programs is None else int(len(self.programs))
        n_prog_sig = None
        if self.programs is not None and len(self.programs) and "significant" in self.programs:
            n_prog_sig = int(self.programs["significant"].fillna(False).sum())
        n_im = (
            None
            if self.programs_induction_matched is None
            else int(len(self.programs_induction_matched))
        )
        n_im_sig = None
        if (
            self.programs_induction_matched is not None
            and len(self.programs_induction_matched)
            and "significant" in self.programs_induction_matched
        ):
            n_im_sig = int(self.programs_induction_matched["significant"].fillna(False).sum())
        return {
            # --- program-first (decisive unit) ---
            "n_programs": n_prog,
            "n_programs_significant": n_prog_sig,
            "n_programs_induction_matched": n_im,
            "n_programs_induction_matched_significant": n_im_sig,
            # --- design / reliability ---
            "n_selected": int(len(self.selected)),
            "n_gene_table": int(len(self.gene_table)),
            "regime": self.regime.get("regime"),
            "reliability": self.regime.get("reliability"),
            "hard_labels_suppressed": bool(
                (self.meta.get("mechanism") or {}).get("hard_labels_suppressed", False)
            ),
            # --- soft per-gene hints (NOT for ORA-by-class) ---
            "per_gene_labels_are_soft": True,
            "per_gene_class_counts_selected": soft_counts,
            # backward-compatible alias
            "class_counts_selected": {
                k: soft_counts[k]
                for k in ("transcription-driven", "stabilization-driven", "ambiguous")
                if k in soft_counts
            },
            "n_induction_confounded_selected": (
                int(self.selected[INDUCTION_CONFOUND_COL].sum())
                if INDUCTION_CONFOUND_COL in self.selected.columns
                else 0
            ),
            "note": (
                "Prefer program_mechanism / program_mechanism_induction_matched over "
                "per-gene mechanism_class. Do not run ORA on genes split by "
                "mechanism_class (induction-confounded at single snapshots)."
            ),
        }


def _resolve_de_membership(
    adata_res: Any,
    all_results: pd.DataFrame,
    de: Any,
    *,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None,
    padj_cutoff: float,
    logfc_cutoff: float,
    logfc_direction: str,
    de_logfc_col: str,
    de_padj_col: str,
    de_kwargs: Mapping[str, Any] | None,
) -> tuple[pd.Index, str, dict[str, Any], pd.DataFrame | None]:
    """Return (selected gene index, de_source label, de diagnostics, de_stats).

    ``de`` dispatch:
      - ``"builtin"``  : select on the DE columns already on ``all_results``.
      - ``str``/``dict``: run :func:`differential_expression` (a de_method name,
        or kwargs) — the package's multi-backend DE front-end.
      - ``DataFrame``  : a precomputed DE table (needs logFC + padj columns).
      - ``callable``   : ``de(adata_res)`` -> such a DataFrame.
    Membership is always ``padj < padj_cutoff`` AND ``logFC > logfc_cutoff``
    in ``logfc_direction`` — the DE step alone decides the list.

    ``de_stats`` is ``None`` for the builtin path (the DE columns already live on
    ``all_results``); for an external DE it is a frame with standardized
    ``logFC`` / ``p_adj`` aligned to ``all_results.index``, so the caller can
    write the SELECTING DE's stats onto the result (reported logFC/p_adj and the
    mechanism direction then match the gates that chose the genes).
    """
    de_diag: dict[str, Any] = {}

    if isinstance(de, str) and de.lower() == _BUILTIN_DE:
        sel = filter_active_genes(
            all_results,
            select_by="de",
            padj_cutoff=padj_cutoff,
            logfc_cutoff=logfc_cutoff,
            logfc_direction=logfc_direction,
        )
        return (
            sel.index,
            "builtin",
            {"de_method": all_results.attrs.get("de_method", "builtin")},
            None,
        )

    # obtain an external DE table
    if callable(de):
        de_df = de(adata_res)
        de_source = "callable"
    elif isinstance(de, pd.DataFrame):
        de_df = de
        de_source = "dataframe"
    elif isinstance(de, (str, Mapping)):
        kw = dict(de_kwargs or {})
        if isinstance(de, str):
            kw.setdefault("de_method", de)
        else:
            kw.update(de)
        de_df = differential_expression(
            adata_res,
            groupby=groupby,
            target_group=target_group,
            reference_group=reference_group,
            sample_col=sample_col,
            **kw,
        )
        de_source = (
            f"differential_expression({kw.get('de_method', kw.get('pseudobulk_de_backend', '?'))})"
        )
        de_diag["de_kwargs"] = kw
    else:
        raise TypeError(
            f"de must be 'builtin', a de_method str, a kwargs dict, a DataFrame, "
            f"or a callable; got {type(de).__name__}"
        )

    # differential_expression / user callables may return (adata, results) — take
    # the DataFrame member.
    if isinstance(de_df, tuple):
        de_df = next((x for x in de_df if isinstance(x, pd.DataFrame)), None)
    if not isinstance(de_df, pd.DataFrame):
        raise TypeError(f"DE source returned {type(de_df).__name__}, expected a DataFrame")
    for col in (de_logfc_col, de_padj_col):
        if col not in de_df.columns:
            raise KeyError(
                f"DE table missing column {col!r}; pass de_logfc_col / de_padj_col "
                f"to map your columns. Available: {list(de_df.columns)}"
            )

    # Align the external DE onto the scored universe and select on standardized
    # columns. Genes absent from the DE table are simply not selectable.
    sel_frame = pd.DataFrame(index=all_results.index)
    sel_frame["logFC"] = pd.to_numeric(
        de_df[de_logfc_col].reindex(all_results.index), errors="coerce"
    )
    sel_frame["p_adj"] = pd.to_numeric(
        de_df[de_padj_col].reindex(all_results.index), errors="coerce"
    )
    # Carry a raw p-value if the external table has one; otherwise NaN — so the
    # caller can clear the builtin DE's stale ``p_val`` rather than leave a value
    # from a different DE backend next to the external logFC/p_adj.
    _pval_col = next(
        (c for c in ("p_val", "pval", "p_value", "pvalue", "pvals") if c in de_df.columns),
        None,
    )
    sel_frame["p_val"] = (
        pd.to_numeric(de_df[_pval_col].reindex(all_results.index), errors="coerce")
        if _pval_col is not None
        else float("nan")
    )
    n_matched = int(sel_frame["logFC"].notna().sum())
    sel = filter_active_genes(
        sel_frame,
        select_by="de",
        padj_cutoff=padj_cutoff,
        logfc_cutoff=logfc_cutoff,
        logfc_direction=logfc_direction,
    )
    de_diag.update({"n_de_genes": int(len(de_df)), "n_matched_to_scored": n_matched})
    return sel.index, de_source, de_diag, sel_frame


def partition_de_by_mechanism(
    adata: Any,
    groupby: str = "condition",
    target_group: str = "Disease",
    reference_group: str = "Control",
    *,
    de: Any = _BUILTIN_DE,
    de_logfc_col: str = "logFC",
    de_padj_col: str = "p_adj",
    de_kwargs: Mapping[str, Any] | None = None,
    sample_col: str | None = None,
    organism: str = "mouse",
    padj_cutoff: float = 0.05,
    logfc_cutoff: float = 1.0,
    logfc_direction: str = "up",
    add_nascent_score: bool = False,
    class_threshold: float | None = None,
    mechanism_preset: str | None = None,
    flag_induction_confound: bool = True,
    induction_confound_penalty: str = "graded",
    suppress_hard_labels_when_unreliable: bool = True,
    min_reliability_for_hard_labels: float = 0.05,
    gene_sets: Mapping[str, Sequence[str]] | None = None,
    program_min_genes: int = 5,
    program_restrict_to_selected: bool = True,
    induction_matched: bool = False,
    induction_matched_methods: Sequence[str] = ("regression", "nearest"),
    run_go_enrichment: bool = False,
    go_gene_sets: str = "GO_Biological_Process",
    show_plot: bool = False,
) -> PartitionResult:
    """DE selects; scATrans partitions by mechanism (primary workflow).

    Runs, in order: a reliability pre-flight, pluggable DE selection, soft
    per-gene mechanism annotation, and — when ``gene_sets`` is given — a
    program-level transcription-vs-stabilization table. The proxy only
    annotates/ranks; it never removes DE hits.

    .. note::
       This **always runs one ``active_score`` pass** for the nascent unspliced-
       excess residual (the mechanism signal). A non-``"builtin"`` ``de=`` only
       replaces the DE **membership and stats** (logFC/p_adj/p_val on the result);
       it does not skip scoring, so bringing your own DE does not avoid that cost.

    Parameters
    ----------
    adata
        AnnData with spliced/unspliced (velocity) layers and a ``groupby`` obs
        column holding ``target_group`` / ``reference_group``.
    de
        DE source that SELECTS the gene list:

        - ``"builtin"`` (default): reuse the DE scATrans already computes.
        - a **de_method** name (``"wilcoxon"``, ``"t-test_overestim_var"``, …)
          or a **kwargs dict** for :func:`~scatrans.tl.differential_expression`
          (e.g. ``{"use_pseudobulk": True, "pseudobulk_de_backend": "pydeseq2"}``).
        - a **precomputed DataFrame** (indexed by gene; needs logFC + padj — map
          names with ``de_logfc_col`` / ``de_padj_col``).
        - a **callable** ``adata -> DataFrame`` of the same shape.
    de_logfc_col, de_padj_col
        Column names for logFC / adjusted-p in an external DE table.
    de_kwargs
        Extra kwargs merged into :func:`differential_expression` (when ``de`` is
        a method name or dict).
    padj_cutoff, logfc_cutoff, logfc_direction
        DE selection thresholds (membership = ``padj < padj_cutoff`` AND
        ``|logFC| > logfc_cutoff`` in the given direction — strict ``>``, matching
        ``filter_active_genes``). Report sensitivity with
        :func:`~scatrans.tl.threshold_sensitivity` rather than defending one.
    add_nascent_score
        When ``True``, append active-transcription **detection** columns from
        :func:`~scatrans.tl.nascent_activity_score` (``nascent_poisson_z``,
        ``dlog_*``, ``de_reproducible`` / ``de_repro_frac``). Detection is
        **decoupled** from mechanism: the Poisson-z is induction-coupled and is
        never used for transcription-vs-stabilization labels (those stay on the
        residual). Additive and fail-soft; inspect
        :func:`~scatrans.qc.regime_diagnosis` on low-capture data.
    class_threshold
        Soft-label boundary (robust-z units) for the per-gene 3-way call.
        Defaults to the ``mechanism_preset`` value, or ``0.5``.
    mechanism_preset
        Passed to :func:`~scatrans.tl.annotate_mechanism_class` as ``preset=``.
        ``"high_precision"`` raises the threshold to 1.0 (lower transcription-driven
        false-positive rate at the same precision).
    flag_induction_confound, induction_confound_penalty
        Passed to :func:`~scatrans.tl.annotate_mechanism_class`. Mark and down-weight
        per-gene stabilization calls in the high-induction (single-snapshot
        non-identifiable) regime; adds an ``induction_confounded`` column and discounts
        ``mechanism_confidence`` without relabeling or touching program-level calls.
        ``penalty="graded"`` (default) or ``"smooth"``.
    suppress_hard_labels_when_unreliable, min_reliability_for_hard_labels
        Passed to :func:`~scatrans.tl.annotate_mechanism_class`. When regime
        reliability is near zero, hard per-gene classes become ``ambiguous``
        (default on).
    gene_sets
        ``{program: [gene, ...]}`` — when given, adds the program-level
        mechanism table (restricted to selected genes by default).
    program_restrict_to_selected
        Pool the program test over the DE-selected genes only (default) vs all
        tested genes.
    induction_matched
        When True and ``gene_sets`` is given, also run
        :func:`~scatrans.tl.program_mechanism_induction_matched` (OLS + nearest
        logFC controls). Recommended for claims when induction varies widely.
    induction_matched_methods
        Methods for the induction-matched table (default ``("regression", "nearest")``).
    run_go_enrichment, go_gene_sets
        Optional GO/pathway ORA on the **DE-selected** genes (not by
        ``mechanism_class``). Enriching genes split by per-gene mechanism labels
        is induction-confounded — use program tables instead.

    Returns
    -------
    PartitionResult
    """
    # Pseudoreplication guard (P0): cell-level DE without sample_col inflates p-values.
    if sample_col is None:
        logger.warning(
            "partition_de_by_mechanism: sample_col is None — DE is cell-level and "
            "subject to pseudoreplication (inflated significance). Prefer "
            "pseudobulk / mixed models with sample_col when biological replicates "
            "exist. Effect sizes (logFC) are more robust than p-values under "
            "pseudoreplication; mechanism labels use support ranks, not DE p."
        )
    elif sample_col not in getattr(adata, "obs", {}):
        logger.warning(
            "partition_de_by_mechanism: sample_col=%r not found in adata.obs — "
            "pseudobulk/sample-aware DE cannot use it.",
            sample_col,
        )

    # 1. score (this is the single active_score pass; also gives the builtin DE)
    adata_res, _significant, all_results = active_score_simple(
        adata,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        organism=organism,
        show_plot=show_plot,
    )

    # 2. reliability pre-flight (MANDATORY — scales confidence, fail-soft). A failed
    #    pre-flight must NOT imply full confidence: default to a cautious 0.5.
    try:
        regime = regime_diagnosis(adata_res)
    except Exception as exc:  # noqa: BLE001 — diagnostic is optional but expected
        logger.warning("regime_diagnosis failed; using cautious reliability=0.5: %s", exc)
        regime = {"reliability": 0.5, "regime": "unknown", "error": str(exc)}
    reliability = float(regime.get("reliability", 0.5))

    # 3. DE SELECTS the gene list (pluggable front-end)
    selected_idx, de_source, de_diag, de_stats = _resolve_de_membership(
        adata_res,
        all_results,
        de,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        padj_cutoff=padj_cutoff,
        logfc_cutoff=logfc_cutoff,
        logfc_direction=logfc_direction,
        de_logfc_col=de_logfc_col,
        de_padj_col=de_padj_col,
        de_kwargs=de_kwargs,
    )
    # For an external DE, the SELECTING DE is the source of truth for the DE
    # columns: write its stats onto all_results so reported values and the
    # mechanism direction match the gates that chose the genes (genes the external
    # DE did not test become NaN → mechanism "unknown", not the builtin DE's stats).
    # p_val is overwritten too (with the external raw p, or NaN) so no stale builtin
    # p_val sits next to the external logFC/p_adj.
    if de_stats is not None:
        all_results["logFC"] = de_stats["logFC"]
        all_results["p_adj"] = de_stats["p_adj"]
        all_results["p_val"] = de_stats["p_val"]

    # 3b. optional ADDITIVE active-transcription DETECTION columns (nascent Poisson-z
    #     + proxy-independent DE-reproducibility flag). This answers a DIFFERENT
    #     question than the mechanism partition and is DELIBERATELY DECOUPLED from it:
    #     the Poisson-z is an *absolute* nascent-increase (induction-coupled), so it
    #     recovers long/high-intron active genes for DETECTION but must NOT drive the
    #     transcription-vs-stabilization call — the mechanism signal stays the
    #     induction-normalized residual (using the z there collapses the ARE/stabilization
    #     signal). Columns are annotation-only; fail-soft.
    nascent_meta: dict[str, Any] = {"enabled": bool(add_nascent_score)}
    if add_nascent_score:
        try:
            nz = nascent_activity_score(
                adata_res,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                sample_col=sample_col,
            )
            for col in nz.columns:
                all_results[col] = nz[col].reindex(all_results.index)
            nascent_meta["status"] = "ok"
            nascent_meta["n_reproducible"] = int(all_results[REPRO_COL].sum())
        except Exception as exc:  # noqa: BLE001 — additive, keep the core result
            logger.warning("add_nascent_score failed (columns omitted): %s", exc)
            nascent_meta = {"enabled": True, "status": "error", "error": str(exc)}

    # 4. per-gene SOFT mechanism annotation on the FULL table (non-destructive;
    #    confidence scaled by the pre-flight reliability). Never gates membership.
    #    Mechanism ALWAYS uses the induction-normalized residual (never the nascent z).
    _, mech_diag = annotate_mechanism_class(
        all_results,
        class_threshold=class_threshold,
        preset=mechanism_preset,
        reliability=reliability,
        suppress_hard_labels_when_unreliable=suppress_hard_labels_when_unreliable,
        min_reliability_for_hard_labels=min_reliability_for_hard_labels,
        flag_induction_confound=flag_induction_confound,
        induction_confound_penalty=induction_confound_penalty,
        inplace=True,
    )
    selected = all_results.loc[selected_idx].copy()
    # stable, interpretable ordering: strongest significance first, then largest
    # effect (ascending p_adj, descending logFC — each independent of the other's
    # presence).
    sort_spec = [
        (c, asc) for c, asc in (("p_adj", True), ("logFC", False)) if c in selected.columns
    ]
    if sort_spec:
        selected = selected.sort_values(
            [c for c, _ in sort_spec], ascending=[asc for _, asc in sort_spec]
        )

    # 5. Program-level mechanism table (threshold-free pooling).
    programs = None
    programs_meta: dict[str, Any] = {"status": "not_requested"}
    programs_im = None
    programs_im_meta: dict[str, Any] = {"status": "not_requested"}
    if gene_sets:
        restrict = list(selected_idx) if program_restrict_to_selected else None
        try:
            programs = program_mechanism(
                all_results,
                gene_sets,
                restrict_index=restrict,
                min_genes=program_min_genes,
            )
            programs_meta = {
                "status": "ok" if len(programs) else "empty",
                "n_programs": int(len(programs)),
                "restrict_to_selected": bool(program_restrict_to_selected),
            }
            if len(programs) == 0:
                programs_meta["reason"] = (
                    "no program met min_genes, or the competitive background was empty "
                    "(all tested genes fall inside the provided gene sets)"
                )
        except Exception as exc:  # noqa: BLE001 — add-on, keep the core result
            logger.warning("program_mechanism skipped: %s", exc)
            programs_meta = {"status": "error", "error": str(exc)}

        if induction_matched:
            try:
                # Induction-matched tests target the induced (up) universe; curated
                # mechanism programs are up-regulated. logfc_min stays 0.0 (logFC > 0);
                # down-regulated programs are out of scope (see limitations).
                programs_im = program_mechanism_induction_matched(
                    all_results,
                    gene_sets,
                    min_genes=program_min_genes,
                    methods=induction_matched_methods,
                    padj_cutoff=padj_cutoff,
                    logfc_min=0.0,
                )
                programs_im_meta = {
                    "status": "ok" if programs_im is not None and len(programs_im) else "empty",
                    "n_programs": 0 if programs_im is None else int(len(programs_im)),
                    "methods": list(induction_matched_methods),
                }
            except Exception as exc:  # noqa: BLE001 — add-on
                logger.warning("program_mechanism_induction_matched skipped: %s", exc)
                programs_im_meta = {"status": "error", "error": str(exc)}

    # 6. optional GO enrichment on the SELECTED genes (DE list — never by class)
    enrichment = None
    if run_go_enrichment and len(selected) > 0:
        from ..enrich import run_enrichment

        enrichment = run_enrichment(
            selected.index.tolist(),
            gene_sets=go_gene_sets,
            organism=organism,
            adata=adata_res,
            pval_cutoff=0.05,
        )

    meta = {
        "scatrans_version": VERSION,
        "organism": organism,
        "de_source": de_source,
        "de": de_diag,
        "sample_col": sample_col,
        "pseudoreplication_warning": sample_col is None,
        "select": {
            "padj_cutoff": padj_cutoff,
            "logfc_cutoff": logfc_cutoff,
            "logfc_direction": logfc_direction,
            "n_selected": int(len(selected)),
        },
        "regime": regime,
        "nascent_score": nascent_meta,
        "mechanism": mech_diag,
        "programs": programs_meta,
        "programs_induction_matched": programs_im_meta,
    }
    logger.info(
        "partition_de_by_mechanism: de=%s selected=%d regime=%s reliability=%.2f "
        "class_counts(all)=%s programs=%s induction_matched=%s",
        de_source,
        len(selected),
        regime.get("regime"),
        reliability,
        mech_diag.get("class_counts"),
        programs_meta.get("status"),
        programs_im_meta.get("status"),
    )
    return PartitionResult(
        adata=adata_res,
        regime=regime,
        gene_table=all_results,
        selected=selected,
        programs=programs,
        programs_induction_matched=programs_im,
        enrichment=enrichment,
        meta=meta,
    )
