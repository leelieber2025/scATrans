"""scatrans.enrich.compare — internal package module."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from typing import Any, Callable

import numpy as np
import pandas as pd

from ._data import (  # noqa: F401 — explicit for type checkers
    BUNDLED_GENE_SET_PROVENANCE,
    GSEA_COLUMNS,
    ORA_COLUMNS,
    _apply_gene_case,
    _apply_p_adjust,
    _bh_p_adjust,
    _bundled_provenance_for,
    _clean_gene_list,
    _DeepcopyImmuneDict,
    _empty_gsea_result,
    _empty_ora_result,
    _expand_gene_list_input,
    _get_analysis_info,
    _load_gene_sets,
    _log_info,
    _open_package_data,
    _parse_gmt_content,
    _resolve_enrichment_padj_cutoff,
    _resolve_gene_set_name,
    _resolve_gseapy_weight,
    _try_load_bundled_gene_set,
    _warn_user,
    list_bundled_gene_sets,
)
from .ora import run_enrichment

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _clean_and_validate_gene_list_for_compare(
    genes: Iterable[Any], gene_case: str | None = None
) -> list[str]:
    """Lightweight cleaner used by compare helpers (re-uses the main cleaner)."""
    return _clean_gene_list(genes, gene_case=gene_case)


def _row_gene_ids_from_df(work: pd.DataFrame) -> list[str]:
    """Return per-row gene identifiers, preferring explicit DE columns over index."""
    return _expand_gene_list_input(work)


def extract_gene_lists(
    de_results: pd.DataFrame | Mapping[str, pd.DataFrame],
    *,
    logfc_cutoff: float = 0.5,
    pval_cutoff: float = 0.05,
    padj_cutoff: float | None = None,
    logfc_direction: str = "up",
    separate_directions: bool = False,
    name_prefix: str | None = None,
    gene_case: str | None = None,
) -> dict[str, list[str]]:
    """
    Extract named gene lists from one or more DE result DataFrames for downstream enrichment.

    This is the recommended way to prepare inputs for `compare_enrichment` when you
    have results from `differential_expression()` or the `all_results` from `active_score()`.

    Supports up / down / both, and can split directions into separate named sets
    (e.g. "GA_up", "GA_down") so you can enrich and visualize them distinctly
    (useful for upset plots and grouped dotplots).

    Parameters
    ----------
    de_results : DataFrame or dict[str, DataFrame]
        - Single DataFrame: treated as one unnamed contrast (you can use name_prefix).
        - dict {contrast_name: de_df}: each df is processed and keys become the cluster names.
    logfc_cutoff : float
        Minimum |logFC| (sign handled by direction).
    pval_cutoff : float
        Legacy name for the **adjusted** p-value cutoff (max p_adj, or p_val if no p_adj).
        Prefer ``padj_cutoff``.
    padj_cutoff : float or None
        Preferred name for the adjusted p-value cutoff. When set, overrides ``pval_cutoff``.
    logfc_direction : {"up", "down", "both"}
    separate_directions : bool
        If True and direction in {"both", "up", "down"} context, will emit separate entries
        "<name>_up" and "<name>_down". Great for "up vs down" enrichment comparison.
    name_prefix : str, optional
        Prepended to generated names when a single df is passed.
    gene_case : optional
        Passed through to gene cleaning.

    Returns
    -------
    dict[str, list[str]]
        Ready to pass to `compare_enrichment(gene_clusters=...)`.

    Example
    -------
    # Multiple contrasts
    de_dict = {"GA_vs_Ctrl": ga_res, "GB_vs_Ctrl": gb_res}
    gene_sets = scat.extract_gene_lists(
        de_dict, logfc_cutoff=0.5, padj_cutoff=0.05, logfc_direction="up"
    )
    comp = scat.compare_enrichment(gene_sets, gene_sets="GO_Biological_Process", organism="mouse")

    # Up and down as separate "clusters" for upset / grouped dotplot
    gene_sets = scat.extract_gene_lists(
        de_dict, logfc_direction="both", separate_directions=True
    )
    """
    # Normalize direction
    dir_raw = str(logfc_direction).lower()
    if dir_raw in ("up", "positive", "pos", "u"):
        direction = "up"
    elif dir_raw in ("down", "negative", "neg", "d"):
        direction = "down"
    elif dir_raw in ("both", "two_sided", "twosided", "abs", "either", "any", "b"):
        direction = "both"
    else:
        raise ValueError(f"logfc_direction={logfc_direction!r} not recognized.")

    # padj_cutoff preferred; pval_cutoff is the legacy alias for adjusted p.
    eff_padj_cutoff = float(padj_cutoff) if padj_cutoff is not None else float(pval_cutoff)

    def _padj_series(work: pd.DataFrame) -> pd.Series:
        """Resolve adjusted p from common DE/export columns; raw p only with warning."""
        adjusted_cols = (
            "p_adj",
            "padj",
            "pvals_adj",
            "p.adjust",
            "FDR",
            "fdr",
            "p_adj_BH",
            "FDR_qval",
        )
        raw_cols = ("p_val", "pval", "pvals", "pvalue", "p")
        for col in adjusted_cols:
            if col in work.columns:
                return pd.to_numeric(work[col], errors="coerce")
        for col in raw_cols:
            if col in work.columns:
                logger.warning(
                    "extract_gene_lists: no *adjusted* p-value column found "
                    "(looked for p_adj/padj/p.adjust/FDR); falling back to raw %r. "
                    "padj_cutoff/pval_cutoff=%.4g will be applied to unadjusted p-values, "
                    "which inflates false positives. Prefer tables with BH/FDR-adjusted p, "
                    "or recompute multipletests before extract_gene_lists.",
                    col,
                    eff_padj_cutoff,
                )
                return pd.to_numeric(work[col], errors="coerce")
        logger.warning(
            "extract_gene_lists: no p-value column found among "
            "p_adj/padj/p.adjust/FDR/p_val/pvalue; treating all p as 1.0 "
            "(likely empty gene lists unless padj_cutoff / pval_cutoff > 1)."
        )
        return pd.Series(1.0, index=work.index)

    def _lfc_series(work: pd.DataFrame) -> pd.Series:
        # Include Seurat FindMarkers avg_log2FC and common aliases.
        lfc_cols = (
            "logFC",
            "logfoldchanges",
            "log2FoldChange",
            "avg_log2FC",
            "avg_logFC",
            "avg_log2FoldChange",
            "lfc",
        )
        for col in lfc_cols:
            if col in work.columns:
                return pd.to_numeric(work[col], errors="coerce")
        logger.warning(
            "extract_gene_lists: no log-fold-change column found among "
            "logFC/logfoldchanges/log2FoldChange/avg_log2FC/lfc; treating logFC as 0.0 "
            "(gene lists will be empty unless logfc_cutoff <= 0). "
            "Rename your effect column or pass a table with a recognized logFC name."
        )
        return pd.Series(0.0, index=work.index)

    def _get_genes_from_df(df: pd.DataFrame) -> list[str]:
        if df is None or df.empty:
            return []
        work = df
        genes = _row_gene_ids_from_df(work)
        padj = _padj_series(work)
        lfc = _lfc_series(work)

        mask = padj < eff_padj_cutoff
        lc = float(logfc_cutoff)

        if direction == "up":
            mask = mask & (lfc > lc)
        elif direction == "down":
            mask = mask & (lfc < -lc)
        else:
            mask = mask & (lfc.abs() > lc)

        selected = [g for g, m in zip(genes, mask) if m]
        return _clean_and_validate_gene_list_for_compare(selected, gene_case=gene_case)

    result: dict[str, list[str]] = {}

    if isinstance(de_results, pd.DataFrame):
        genes = _get_genes_from_df(de_results)
        base_name = name_prefix or "contrast"
        if separate_directions:
            # Recompute up and down separately
            work = de_results
            padj = _padj_series(work)
            lfc = _lfc_series(work)
            sig = padj < eff_padj_cutoff
            lc = float(logfc_cutoff)
            genes_idx = _row_gene_ids_from_df(work)
            up = [g for g, keep, lf in zip(genes_idx, sig, lfc.fillna(0)) if keep and lf > lc]
            down = [g for g, keep, lf in zip(genes_idx, sig, lfc.fillna(0)) if keep and lf < -lc]
            if name_prefix:
                result[f"{name_prefix}_up"] = _clean_and_validate_gene_list_for_compare(
                    up, gene_case
                )
                result[f"{name_prefix}_down"] = _clean_and_validate_gene_list_for_compare(
                    down, gene_case
                )
            else:
                result["up"] = _clean_and_validate_gene_list_for_compare(up, gene_case)
                result["down"] = _clean_and_validate_gene_list_for_compare(down, gene_case)
        else:
            result[base_name] = genes
        return result

    # dict case
    for cname, df in de_results.items():
        if df is None or (hasattr(df, "empty") and df.empty):
            continue
        base = str(cname)
        if name_prefix:
            base = f"{name_prefix}_{base}"
        if separate_directions:
            # split up/down for this contrast
            work = df
            padj = _padj_series(work)
            lfc = _lfc_series(work)
            sig = padj < eff_padj_cutoff
            lc = float(logfc_cutoff)
            genes_idx = _row_gene_ids_from_df(work)
            up = [g for g, keep, lf in zip(genes_idx, sig, lfc.fillna(0)) if keep and lf > lc]
            down = [g for g, keep, lf in zip(genes_idx, sig, lfc.fillna(0)) if keep and lf < -lc]
            result[f"{base}_up"] = _clean_and_validate_gene_list_for_compare(up, gene_case)
            result[f"{base}_down"] = _clean_and_validate_gene_list_for_compare(down, gene_case)
        else:
            genes = _get_genes_from_df(df)
            result[base] = genes
    return result


def concat_compare_results(
    results: Mapping[str, pd.DataFrame] | list[tuple[str, pd.DataFrame]],
    cluster_col: str = "Cluster",
) -> pd.DataFrame:
    """
    Combine already-computed per-group enrichment results into a compare-style DataFrame.

    Each input must be a result from run_enrichment / run_kegg / run_go / run_gsea (or similar).
    The keys (or first element of tuple) become the "Cluster" value.

    This lets users who manually looped over contrasts easily get a single table
    with "Cluster" column for plotting.
    """
    items = list(results.items()) if isinstance(results, Mapping) else list(results)

    frames = []
    clusters_in_frames: list[str] = []
    per_cluster_attrs: dict[str, Any] = {}
    for name, df in items:
        if df is None or (hasattr(df, "empty") and df.empty):
            continue
        cname = str(name)
        out = df.copy()
        if cluster_col not in out.columns:
            out.insert(0, cluster_col, cname)
        frames.append(out)
        clusters_in_frames.append(cname)
        if hasattr(df, "attrs"):
            per_cluster_attrs[cname] = dict(df.attrs)

    if not frames:
        empty = _empty_ora_result()
        empty.attrs["method"] = "compare_concat"
        empty.attrs["cluster_col"] = cluster_col
        return empty

    combined = pd.concat(frames, ignore_index=True, sort=False)
    # Put Cluster first if present
    if cluster_col in combined.columns:
        cols = [cluster_col] + [c for c in combined.columns if c != cluster_col]
        combined = combined[cols]

    combined.attrs["method"] = "compare_concat"
    combined.attrs["cluster_col"] = cluster_col
    combined.attrs["per_cluster_attrs"] = per_cluster_attrs
    combined.attrs["n_clusters"] = len(frames)
    combined.attrs["clusters"] = clusters_in_frames
    return combined


def compare_enrichment(
    gene_clusters: Mapping[str, Iterable[Any]],
    *,
    fun: Callable | None = None,
    pval_cutoff: float | None = None,
    padj_cutoff: float | None = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    organism: str = "mouse",
    gene_sets: Any = "GO_Biological_Process",
    universe: Iterable[Any] | None = None,
    adata: Any | None = None,
    gene_set_source: str = "scatrans",
    gene_case: str | None = None,
    verbose: bool = True,
    return_all: bool = False,
    raise_on_error: bool = False,
    adjust_across_clusters: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Run enrichment analysis across multiple named groups/clusters/contrasts
    (clusterProfiler::compareCluster style).

    **Important default note on multiple-testing**:
    ``adjust_across_clusters=False`` (the default) adjusts p-values
    *separately inside each cluster* using ``p_adjust_method`` (default BH FDR)
    when the per-cluster callable supports it (e.g. :func:`run_enrichment`).
    This is conservative per group but when you have many clusters the *overall*
    false discovery rate across the whole table can be higher than a single
    global correction. If you intend to compare significance across clusters,
    use ``adjust_across_clusters=True``; per-cluster calls then use
    ``return_all=True`` internally (all size-eligible terms, including
    zero-overlap) before a single global re-adjustment with the **same**
    ``p_adjust_method`` (not a silent switch to BH), and the final table is
    filtered by ``padj_cutoff``/``pval_cutoff`` unless you also pass
    ``return_all=True``.

    Biological comparability notes (inspired by clusterProfiler best practices):
    - When `adata=` or `universe=` is supplied, the **same background** is used for every cluster.
      This is critical for fair comparison across groups.
    - Term filtering (min_size etc.) and p-adjustment are done per-cluster (standard and
      conservative). You can post-filter the returned table.
    - A 'Cluster' column is always added as the first column.
    - Rich per-cluster metadata (including failures and empty) is stored under
      `.attrs["per_cluster"]` and also under `.attrs["scatrans"]["per_cluster"]`.

    Parameters
    ----------
    ...
    raise_on_error : bool, default False
        If True, any exception from a per-cluster enrichment call will be re-raised
        immediately (good for debugging). If False (default), the cluster is recorded
        as failed and execution continues (good for large batch runs).

    Returns
    -------
    pd.DataFrame
        Concatenated results with a leading "Cluster" column.
        Compatible with the enhanced `scat.pl.enrich_dotplot(..., cluster_col="Cluster")`
        and `scat.pl.enrich_upsetplot(...)` / `enrich_vennplot`.

    Example
    -------
    ...
    """
    if not isinstance(gene_clusters, Mapping) or len(gene_clusters) == 0:
        raise ValueError("gene_clusters must be a non-empty dict-like {name: gene_list}")

    if fun is None:
        fun = run_enrichment

    # Resolve a shared universe once if possible (key for biological comparability)
    shared_universe = universe
    shared_adata = adata

    # If user gave adata but no explicit universe, we let the subcalls use the smart default,
    # which is the desired "same background" behavior when the same adata object is passed.
    if shared_universe is None and shared_adata is None and verbose:
        # Nothing shared provided - each cluster will use its own default (usually gene set genes).
        _log_info(
            "compare_enrichment: no shared `universe` or `adata` provided. "
            "Each cluster will determine its own background. "
            "For best cross-cluster comparability, pass the same adata (after store_raw_counts) or explicit universe."
        )

    eff_cut = _resolve_enrichment_padj_cutoff(pval_cutoff, padj_cutoff)
    # Same method as per-cluster enrichment when fun=run_enrichment (via kwargs).
    p_adjust_method = str(kwargs.get("p_adjust_method", "fdr_bh"))
    eff_return_all = return_all
    if adjust_across_clusters:
        if not return_all:
            _warn_user(
                "compare_enrichment: adjust_across_clusters=True requires every tested term "
                "from each cluster (including zero-overlap and non-significant) before global "
                f"p-value adjustment (p_adjust_method={p_adjust_method!r}). "
                "Forcing return_all=True on per-cluster calls; significance "
                "filtering is applied after global p.adjust."
            )
        eff_return_all = True

    frames = []
    clusters_with_data: list[str] = []
    per_cluster: dict[str, Any] = {}
    cluster_names = list(gene_clusters.keys())

    for cname in cluster_names:
        cname_str = str(cname)
        # Materialize once: gene_clusters[cname] may be a generator/iterator.
        # Consuming it in _clean... and then again for n_raw would give n_raw=0.
        raw_items = list(gene_clusters[cname])
        n_raw = len(raw_items)
        genes = _clean_and_validate_gene_list_for_compare(raw_items, gene_case=gene_case)
        n_clean = len(genes)

        if not genes:
            per_cluster[cname_str] = {
                "n_input_genes_raw": n_raw,
                "n_input_genes_clean": 0,
                "skipped": "no genes after cleaning",
            }
            if verbose:
                _log_info(
                    f"compare_enrichment: cluster '{cname}' has no genes after cleaning; skipping"
                )
            continue

        try:
            res = fun(
                gene_list=genes,
                pval_cutoff=pval_cutoff,
                padj_cutoff=padj_cutoff,
                min_size=min_size,
                max_size=max_size,
                restrict_background_to_gene_sets=restrict_background_to_gene_sets,
                force_universe=force_universe,
                organism=organism,
                gene_sets=gene_sets,
                universe=shared_universe,
                adata=shared_adata,
                gene_set_source=gene_set_source,
                gene_case=gene_case,
                return_all=eff_return_all,
                verbose=verbose,
                **kwargs,
            )
        except Exception as e:
            per_cluster[cname_str] = {
                "n_input_genes_raw": n_raw,
                "n_input_genes_clean": n_clean,
                "failed": True,
                "error": str(e),
            }
            if raise_on_error:
                raise
            _warn_user(f"compare_enrichment: failed on cluster '{cname}': {e}")
            continue

        if res is None or (hasattr(res, "empty") and res.empty):
            per_cluster[cname_str] = {
                "n_input_genes_raw": n_raw,
                "n_input_genes_clean": n_clean,
                "empty": True,
                "attrs": getattr(res, "attrs", {}) if res is not None else {},
            }
            if verbose:
                _log_info(f"compare_enrichment: cluster '{cname}' returned empty result")
            continue

        out = res.copy()
        if "Cluster" not in out.columns:
            out.insert(0, "Cluster", cname_str)

        frames.append(out)
        clusters_with_data.append(cname_str)
        per_cluster[cname_str] = {
            "n_input_genes_raw": n_raw,
            "n_input_genes_clean": n_clean,
            "n_terms": len(out),
            "attrs": dict(getattr(res, "attrs", {})),
        }

    if not frames:
        # Return a properly attributed empty (consistent columns when possible)
        empty = _empty_ora_result()
        empty.insert(0, "Cluster", pd.Series(dtype=str))  # keep Cluster column for API stability
        empty.attrs.update(
            {
                "method": "compare_enrichment",
                "cluster_col": "Cluster",
                "n_clusters_attempted": len(cluster_names),
                "per_cluster": per_cluster,
                "shared_universe_used": bool(shared_universe or shared_adata),
                "clusterprofiler_aligned": True,
            }
        )
        empty.attrs.setdefault("scatrans", {})
        empty.attrs["scatrans"].update(
            {
                "method": "compare_enrichment",
                "per_cluster": per_cluster,
                "n_clusters_attempted": len(cluster_names),
                "adjust_across_clusters": bool(adjust_across_clusters),
                "multiple_testing": {
                    "scope": "per_cluster",
                    "method": p_adjust_method,
                    "n_tests": 0,
                },
            }
        )
        return empty

    combined = pd.concat(frames, ignore_index=True, sort=False)

    # Ensure Cluster is the very first column for convenience
    if "Cluster" in combined.columns:
        other_cols = [c for c in combined.columns if c != "Cluster"]
        combined = combined[["Cluster"] + other_cols]

    # Rich metadata
    combined.attrs["method"] = "compare_enrichment"
    combined.attrs["cluster_col"] = "Cluster"
    combined.attrs["n_clusters"] = len(frames)
    combined.attrs["clusters"] = list(clusters_with_data)
    combined.attrs["per_cluster"] = per_cluster
    combined.attrs["shared_universe_used"] = bool(shared_universe or shared_adata)
    combined.attrs["clusterprofiler_aligned"] = True
    combined.attrs["p_adjust_method"] = p_adjust_method

    # Multi-cluster multiple testing (analogous to adjust_across_all in run_go).
    # Use the same p_adjust_method as per-cluster runs — never hardcode BH.
    multiple_testing_scope = "per_cluster"
    if adjust_across_clusters and not combined.empty and "pvalue" in combined.columns:
        combined = combined.copy()
        if "p.adjust" in combined.columns:
            combined["p.adjust.within_cluster"] = combined["p.adjust"]
        combined["p.adjust"] = _apply_p_adjust(combined["pvalue"].values, method=p_adjust_method)
        combined["neg_log10_padj"] = -np.log10(
            combined["p.adjust"].astype(float).clip(lower=1e-300)
        )
        combined = combined.sort_values("p.adjust").reset_index(drop=True)
        multiple_testing_scope = "all_clusters"
        if verbose:
            _log_info(
                "compare_enrichment: p.adjust re-computed across all clusters "
                f"(adjust_across_clusters=True, p_adjust_method={p_adjust_method!r}; "
                f"n_tests={len(combined)})"
            )

    # Structured under scatrans (fixed from previous "scatans" note)
    combined.attrs.setdefault("scatrans", {})
    combined.attrs["scatrans"].update(
        {
            "method": "compare_enrichment",
            "per_cluster": per_cluster,
            "n_clusters_attempted": len(cluster_names),
            "raise_on_error": raise_on_error,
            "adjust_across_clusters": bool(adjust_across_clusters),
            "multiple_testing": {
                "scope": multiple_testing_scope,
                "method": p_adjust_method,
                "n_tests": len(combined) if "pvalue" in combined.columns else 0,
            },
        }
    )
    # Also keep flat for compatibility
    combined.attrs["multiple_testing_scope"] = multiple_testing_scope

    if verbose:
        _log_info(
            f"compare_enrichment: {len(frames)} clusters, total {len(combined)} rows. "
            f"Shared background: {bool(shared_universe or shared_adata)}"
        )

    if not return_all and not combined.empty and "p.adjust" in combined.columns:
        sig = combined[combined["p.adjust"] < eff_cut].copy().reset_index(drop=True)
        sig.attrs.update(combined.attrs)
        if "Cluster" in sig.columns:
            clusters_returned = sig["Cluster"].astype(str).unique().tolist()
            sig.attrs["clusters"] = clusters_returned
            sig.attrs["n_clusters"] = len(clusters_returned)
        else:
            sig.attrs["clusters"] = []
            sig.attrs["n_clusters"] = 0
        if verbose:
            _log_info(
                f"compare_enrichment: returning {len(sig)} significant rows at padj < {eff_cut}"
            )
        return sig

    return combined
