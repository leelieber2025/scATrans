"""scatrans.tl.compare — internal package module.

Compare differential expression across levels of a grouping variable (e.g.
cell type / cluster), running the same contrast (``groupby`` /
``target_group`` / ``reference_group``) independently within each level.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from .._utils import _categorical_or_unique
from ._common import COMPARE_LOGFC_CUTOFF
from .de import differential_expression
from .filter import filter_active_genes

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@dataclass
class CompareDEResult:
    """Result of :func:`compare_de_across_groups`.

    Attributes
    ----------
    results
        Full DE table (as returned by :func:`differential_expression`) for
        every level of ``split_by`` where DE ran successfully.
    up
        Upregulated genes per level, filtered with :func:`filter_active_genes`
        (``select_by="de"``, ``logfc_direction="up"``) using the thresholds
        recorded in ``meta``.
    down
        Downregulated genes per level, same filter with
        ``logfc_direction="down"``.
    summary
        One row per successful level of ``split_by`` with ``up`` / ``down`` /
        ``total`` gene counts, indexed by group name. Ready to hand to
        :func:`scatrans.pl.de_summary_barplot`.
    failed
        Mapping of ``split_by`` level -> error message, for levels skipped
        because DE could not be run (too few cells, no replicates, backend
        error, etc.).
    meta
        Thresholds, DE kwargs, and the contrast used, for provenance.
    """

    results: dict[str, pd.DataFrame]
    up: dict[str, pd.DataFrame]
    down: dict[str, pd.DataFrame]
    summary: pd.DataFrame
    failed: dict[str, str] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)


def compare_de_across_groups(
    adata: Any,
    split_by: str,
    groupby: str,
    target_group: str,
    reference_group: str,
    *,
    groups: Iterable[str] | None = None,
    min_cells: int = 10,
    padj_cutoff: float | None = 0.05,
    pval_cutoff: float | None = None,
    p_type: str = "auto",
    logfc_cutoff: float = COMPARE_LOGFC_CUTOFF,
    de_kwargs: Mapping[str, Any] | None = None,
    raise_on_error: bool = False,
) -> CompareDEResult:
    """
    Run the same DE contrast independently within every level of ``split_by``
    (e.g. cell type or cluster) and summarize up/down-regulated gene counts
    for comparison across levels.

    For each level ``g`` of ``adata.obs[split_by]``, this subsets
    ``adata[adata.obs[split_by] == g]`` and calls
    :func:`differential_expression` with ``groupby``/``target_group``/
    ``reference_group`` plus any ``de_kwargs`` (so any DE method the package
    supports — scanpy rank_genes_groups, pseudobulk PyDESeq2, mixed models,
    Memento — can be used; pass the relevant options through ``de_kwargs``,
    e.g. ``de_kwargs={"use_pseudobulk": True, "sample_col": "individual",
    "pseudobulk_de_backend": "pydeseq2"}``). Significant up/down genes are
    then selected with :func:`filter_active_genes` (``select_by="de"``) using
    ``padj_cutoff``/``pval_cutoff``/``p_type``/``logfc_cutoff``.

    Levels with fewer than ``min_cells`` cells, or where DE raises an
    exception (e.g. too few replicates for pseudobulk dispersion estimation),
    are skipped and recorded in ``result.failed`` rather than aborting the
    whole comparison — set ``raise_on_error=True`` to fail fast instead.

    Parameters
    ----------
    adata
        AnnData object containing all levels of ``split_by``.
    split_by
        ``adata.obs`` column to iterate over (e.g. ``"cell_type"``).
    groupby, target_group, reference_group
        The DE contrast, forwarded to :func:`differential_expression` for
        every ``split_by`` level (e.g. ``groupby="sample"``,
        ``target_group="RTX"``, ``reference_group="VEH"``).
    groups
        Explicit ordered subset of ``split_by`` levels to run (default: all
        categories/unique values, in their existing order).
    min_cells
        Minimum number of cells a ``split_by`` level must have to attempt DE.
    padj_cutoff, pval_cutoff, p_type, logfc_cutoff
        Significance thresholds forwarded to :func:`filter_active_genes`
        (``p_type="auto"`` prefers ``p_adj`` when present, else ``p_val``).
    de_kwargs
        Extra keyword arguments forwarded to :func:`differential_expression`
        (e.g. ``de_method``, ``use_pseudobulk``, ``sample_col``,
        ``pseudobulk_de_backend``, ``use_mixed_model``, ``use_memento_de``).
    raise_on_error
        If True, re-raise immediately on the first level that fails instead
        of skipping it and continuing.

    Returns
    -------
    CompareDEResult

    Examples
    --------
    >>> cmp = scat.compare_de_across_groups(
    ...     adata, split_by="cell_type", groupby="sample",
    ...     target_group="RTX", reference_group="VEH",
    ...     de_kwargs={"use_pseudobulk": True, "sample_col": "individual",
    ...                "pseudobulk_de_backend": "pydeseq2"},
    ... )
    >>> cmp.summary
    >>> scat.pl.de_summary_barplot(cmp.summary, mode="stacked")
    """
    if split_by not in adata.obs.columns:
        raise ValueError(f"split_by={split_by!r} not found in adata.obs")

    col = adata.obs[split_by]
    groups = _categorical_or_unique(col, drop_unused=False) if groups is None else list(groups)

    de_kwargs = dict(de_kwargs or {})
    filt_kwargs: dict[str, Any] = {
        "select_by": "de",
        "logfc_cutoff": logfc_cutoff,
        "p_type": p_type,
    }
    if padj_cutoff is not None:
        filt_kwargs["padj_cutoff"] = padj_cutoff
    if pval_cutoff is not None:
        filt_kwargs["pval_cutoff"] = pval_cutoff

    results: dict[str, pd.DataFrame] = {}
    up: dict[str, pd.DataFrame] = {}
    down: dict[str, pd.DataFrame] = {}
    failed: dict[str, str] = {}
    summary_rows: list[dict[str, Any]] = []

    for grp in groups:
        gname = str(grp)
        mask = (col == grp).to_numpy()
        n_cells = int(mask.sum())
        if n_cells < min_cells:
            msg = f"only {n_cells} cell(s) in this level (< min_cells={min_cells})"
            failed[gname] = msg
            if raise_on_error:
                raise RuntimeError(f"compare_de_across_groups: skipping {gname!r}: {msg}")
            logger.warning("compare_de_across_groups: skipping %r (%s)", gname, msg)
            continue

        logger.info("compare_de_across_groups: running DE for %r (%d cells)", gname, n_cells)
        try:
            adata_sub = adata[mask].copy()
            _, de_result = differential_expression(
                adata_sub,
                groupby=groupby,
                target_group=target_group,
                reference_group=reference_group,
                **de_kwargs,
            )
        except Exception as e:
            failed[gname] = str(e)
            if raise_on_error:
                raise
            logger.warning("compare_de_across_groups: DE failed for %r: %s", gname, e)
            continue

        results[gname] = de_result
        up_df = filter_active_genes(de_result, logfc_direction="up", **filt_kwargs)
        down_df = filter_active_genes(de_result, logfc_direction="down", **filt_kwargs)
        up[gname] = up_df
        down[gname] = down_df
        summary_rows.append(
            {
                "group": gname,
                "up": int(len(up_df)),
                "down": int(len(down_df)),
                "total": int(len(up_df) + len(down_df)),
            }
        )

    if summary_rows:
        summary = pd.DataFrame(summary_rows).set_index("group")
    else:
        summary = pd.DataFrame(columns=["up", "down", "total"]).rename_axis("group")

    if failed:
        logger.warning(
            "compare_de_across_groups: %d/%d level(s) of %r skipped: %s",
            len(failed),
            len(groups),
            split_by,
            failed,
        )

    return CompareDEResult(
        results=results,
        up=up,
        down=down,
        summary=summary,
        failed=failed,
        meta={
            "split_by": split_by,
            "groupby": groupby,
            "target_group": target_group,
            "reference_group": reference_group,
            "padj_cutoff": padj_cutoff,
            "pval_cutoff": pval_cutoff,
            "p_type": p_type,
            "logfc_cutoff": logfc_cutoff,
            "de_kwargs": de_kwargs,
        },
    )
