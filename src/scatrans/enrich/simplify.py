"""scatrans.enrich.simplify — internal package module."""

from __future__ import annotations

import logging
import math
import re
from collections.abc import Iterable, Mapping
from fractions import Fraction
from typing import Any

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

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _comb_fraction(n: int, k: int) -> Fraction:
    """Binomial coefficient C(n, k) as an exact Fraction (small-n / tests).

    Prefer :func:`_comb_comb_comb` (scipy hypergeom) for production PathwayDenester
    calls — full factorials are O(n) multiplications and hang on large GO terms.
    """
    if k > n or k < 0:
        return Fraction(0)
    # math.comb is iterative product (avoids computing n! fully when k is small)
    return Fraction(math.comb(n, k))


def _comb_comb_comb(
    degs_in_test: int,
    degs_in_intersection: int,
    intersection_size: int,
    size_test: int,
) -> float:
    """PathwayDenester independence test (upper tail of a hypergeometric).

    Equivalent to the original exact combinatorial sum over
    ``C(K,k)*C(N-K,n-k)/C(N,n)`` for ``k >= degs_in_intersection``, but uses
    :func:`scipy.stats.hypergeom.sf` so large pathway sizes stay tractable.
    """
    if intersection_size <= 0 or size_test <= 0 or degs_in_test <= 0:
        return 1.0
    if degs_in_intersection <= 0:
        return 1.0
    # Guard invalid hypergeometric parameterizations
    if (
        degs_in_test > size_test
        or intersection_size > size_test
        or degs_in_intersection > min(intersection_size, degs_in_test)
    ):
        return 1.0
    try:
        from scipy.stats import hypergeom
    except ImportError:
        # Fallback: exact Fraction sum (only for tiny n if scipy missing)
        denominator = _comb_fraction(size_test, intersection_size)
        if denominator == 0:
            return 1.0
        p_sum = Fraction(0)
        upper = min(intersection_size, degs_in_test)
        for desired_number_of_degs in range(degs_in_intersection, upper + 1):
            numerator = _comb_fraction(
                size_test - degs_in_test, intersection_size - desired_number_of_degs
            ) * _comb_fraction(degs_in_test, desired_number_of_degs)
            p_sum += numerator / denominator
        return float(p_sum)

    # P(X >= k) = sf(k - 1) for X ~ Hypergeometric(N, K, n)
    return float(
        hypergeom.sf(
            degs_in_intersection - 1,
            size_test,
            degs_in_test,
            intersection_size,
        )
    )


def _parse_gene_overlap_field(genes_str: str) -> list[str]:
    return [g.strip() for g in re.split(r"[;,]+", str(genes_str)) if g.strip()]


def _resolve_gene_sets_for_simplify(
    enrich_df: pd.DataFrame,
    gene_sets: Mapping[str, Iterable[Any]] | str | None,
    organism: str = "mouse",
    verbose: bool = True,
) -> dict[str, set]:
    """Resolve full pathway gene memberships for PathwayDenester."""
    if gene_sets is not None:
        term_to_genes, _, _ = _load_gene_sets(gene_sets, organism=organism, verbose=verbose)
        return term_to_genes

    attrs = getattr(enrich_df, "attrs", {}) or {}
    gene_set_info = attrs.get("gene_set_info") or {}
    for key in ("resolved", "requested", "library_name"):
        candidate = gene_set_info.get(key)
        if candidate and candidate != "<dict>":
            try:
                term_to_genes, _, _ = _load_gene_sets(candidate, organism=organism, verbose=verbose)
                if verbose:
                    _log_info(
                        f"PathwayDenester: loaded gene sets from enrichment attrs "
                        f"({key}='{candidate}')"
                    )
                return term_to_genes
            except Exception:
                continue

    requested = gene_set_info.get("requested")
    extra = ""
    if requested == "<dict>" or isinstance(gene_sets, dict):
        extra = (
            " Custom dict-based enrichment does not store full pathway gene memberships in "
            "the result table — pass the same `gene_sets=` dict (or GMT/library name) to "
            "simplify_enrichment(..., method='pathway_denester', gene_sets=...)."
        )
    raise ValueError(
        "method='pathway_denester' requires `gene_sets` (GMT path, bundled name, or dict) "
        "or an enrichment result whose .attrs['gene_set_info'] can be resolved." + extra
    )


def _simplify_by_pathway_denester(
    df: pd.DataFrame,
    term_to_genes: dict[str, set],
    *,
    gene_col: str,
    p_col: str,
    to_test_threshold: float = 0.0,
    pval_threshold: float = 0.05,
    term_size_limit: int = 0,
    show_excluded: bool = False,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    PathwayDenester-style nested pathway filtering.

    Adapted from PathwayDenester v3.7 (Helmy-Lab):
    https://github.com/Helmy-Lab/PathwayDenester
    """
    if not (0 <= to_test_threshold <= 1):
        raise ValueError("to_test_threshold must be between 0 and 1")
    if not (0 <= pval_threshold <= 1):
        raise ValueError("pval_threshold must be between 0 and 1")

    work = df.copy()
    if "Term" not in work.columns:
        raise ValueError("PathwayDenester requires a 'Term' column with pathway IDs.")

    name_col = "Description" if "Description" in work.columns else "Term"
    size_col = "TermSize" if "TermSize" in work.columns else None

    work["_pd_p_value"] = pd.to_numeric(work[p_col], errors="coerce")
    work = work.dropna(subset=["_pd_p_value"]).copy()
    if work.empty:
        return work.drop(columns=["_pd_p_value"], errors="ignore")

    if "Count" in work.columns:
        work["_pd_intersection_size"] = pd.to_numeric(work["Count"], errors="coerce")
    else:
        work["_pd_intersection_size"] = work[gene_col].map(
            lambda x: len(_parse_gene_overlap_field(x))
        )

    if size_col is not None:
        work["_pd_term_size"] = pd.to_numeric(work[size_col], errors="coerce")
        work["_pd_ratio"] = work["_pd_intersection_size"] / work["_pd_term_size"].replace(0, np.nan)
        work = work.sort_values(
            by=["_pd_intersection_size", "_pd_ratio", "_pd_p_value"],
            ascending=[False, False, True],
            kind="stable",
        )
    else:
        work = work.sort_values(
            by=["_pd_intersection_size", "_pd_p_value"],
            ascending=[False, True],
            kind="stable",
        )
    work = work.reset_index(drop=True)

    missing_terms = work.loc[~work["Term"].isin(term_to_genes), "Term"].tolist()
    if missing_terms and verbose:
        preview = ", ".join(str(t) for t in missing_terms[:5])
        suffix = "..." if len(missing_terms) > 5 else ""
        _warn_user(
            f"PathwayDenester: {len(missing_terms)} terms not found in gene_sets and will be skipped. "
            f"Examples: {preview}{suffix}"
        )
    work = work[work["Term"].isin(term_to_genes)].copy()
    if work.empty:
        return work.drop(
            columns=["_pd_p_value", "_pd_intersection_size", "_pd_term_size", "_pd_ratio"],
            errors="ignore",
        )

    if term_size_limit > 1:
        work = work[
            work["Term"].map(lambda t: len(term_to_genes.get(str(t), set()))) <= term_size_limit
        ].copy()
        if work.empty:
            return work.drop(
                columns=["_pd_p_value", "_pd_intersection_size", "_pd_term_size", "_pd_ratio"],
                errors="ignore",
            )

    pathways: list[dict[str, Any]] = []
    for row_idx, row in work.iterrows():
        term_id = str(row["Term"])
        all_genes = set(term_to_genes.get(term_id, set()))
        entry: dict[str, Any] = {
            "row_index": int(row_idx),
            "id": term_id,
            "name": str(row.get(name_col, term_id)),
            "p_value": float(row["_pd_p_value"]),
            "all_genes": all_genes,
            "deg_list": [],
            "degs": set(),
            "density": 0.0,
            "result": 1.0,
            "filter": "exclude",
            "vs": "itself",
            "vs_name": "",
            "reciprocal": 0.0,
            "intersection_size": 0,
            "degs_in_intersection": 0,
            "testable": False,
        }
        if not all_genes:
            pathways.append(entry)
            continue
        raw_degs = _parse_gene_overlap_field(row.get(gene_col, ""))
        deg_list = [g for g in raw_degs if g in all_genes]
        if not deg_list:
            pathways.append(entry)
            continue
        unexpected = set(raw_degs) - all_genes
        if unexpected:
            ratio_unexpected = len(unexpected) / max(len(raw_degs), 1)
            if ratio_unexpected > 0.1:
                if verbose:
                    _warn_user(
                        f"PathwayDenester: term '{term_id}' has "
                        f"{ratio_unexpected:.1%} DEGs absent from gene_sets; skipped."
                    )
                pathways.append(entry)
                continue
        entry.update(
            {
                "deg_list": deg_list,
                "degs": set(deg_list),
                "density": len(deg_list) / len(all_genes),
                "filter": "keep",
                "testable": True,
            }
        )
        pathways.append(entry)

    testable_indices = [i for i, p in enumerate(pathways) if p["testable"]]

    for pos, current_idx in enumerate(testable_indices):
        if pos == 0:
            continue
        current = pathways[current_idx]
        degs_in_current = len(current["degs"])
        size_current = len(current["all_genes"])
        for test_pos in range(pos):
            test = pathways[testable_indices[test_pos]]
            if test["filter"] != "keep":
                continue
            degs_in_test = len(test["degs"])
            shared_genes = test["all_genes"] & current["all_genes"]
            degs_in_intersection_current = len(current["degs"] & shared_genes)
            if degs_in_intersection_current <= to_test_threshold * min(
                degs_in_current, degs_in_test
            ):
                continue
            intersection_size = len(shared_genes)
            size_test = len(test["all_genes"])
            current_result = _comb_comb_comb(
                degs_in_current,
                degs_in_intersection_current,
                intersection_size,
                size_current,
            )
            reverse_result = _comb_comb_comb(
                degs_in_test,
                degs_in_intersection_current,
                intersection_size,
                size_test,
            )
            if current_result < pval_threshold and reverse_result > pval_threshold:
                current["result"] = current_result
                current["reciprocal"] = reverse_result
                current["filter"] = "exclude"
                current["vs"] = test["id"]
                current["vs_name"] = test["name"]
                current["intersection_size"] = intersection_size
                current["degs_in_intersection"] = degs_in_intersection_current
                break
            if current_result < current["result"] and current["filter"] == "keep":
                current["vs"] = test["id"]
                current["vs_name"] = test["name"]
                current["result"] = current_result
                current["reciprocal"] = reverse_result
                current["intersection_size"] = intersection_size
                current["degs_in_intersection"] = degs_in_intersection_current

    out = work.copy()
    out["Denester_filter"] = [p["filter"] for p in pathways]
    out["Denester_result"] = [p["result"] for p in pathways]
    out["Denester_reciprocal"] = [p["reciprocal"] for p in pathways]
    out["Denester_vs"] = [p["vs"] for p in pathways]
    out["Denester_vs_name"] = [p["vs_name"] for p in pathways]
    out["Denester_intersection_size"] = [p["intersection_size"] for p in pathways]
    out["Denester_degs_in_intersection"] = [p["degs_in_intersection"] for p in pathways]

    if show_excluded:
        result = out.drop(
            columns=["_pd_p_value", "_pd_intersection_size", "_pd_term_size", "_pd_ratio"],
            errors="ignore",
        )
    else:
        kept_mask = out["Denester_filter"] == "keep"
        result = out.loc[kept_mask].drop(
            columns=["_pd_p_value", "_pd_intersection_size", "_pd_term_size", "_pd_ratio"],
            errors="ignore",
        )

    if verbose:
        n_kept = int((result["Denester_filter"] == "keep").sum()) if show_excluded else len(result)
        _log_info(
            f"Simplified from {len(df)} to {n_kept} terms "
            f"(PathwayDenester, pval_threshold={pval_threshold}, "
            f"to_test_threshold={to_test_threshold})"
        )

    result = result.reset_index(drop=True)
    if hasattr(df, "attrs"):
        result.attrs.update(dict(df.attrs))
    result.attrs["simplify_method"] = "pathway_denester"
    result.attrs["pathway_denester_params"] = {
        "to_test_threshold": float(to_test_threshold),
        "pval_threshold": float(pval_threshold),
        "term_size_limit": int(term_size_limit),
        "show_excluded": bool(show_excluded),
        "p_col": p_col,
    }
    return result


def simplify_enrichment(
    enrich_df: pd.DataFrame,
    similarity_cutoff: float = 0.5,
    by: str | None = None,
    ascending: bool = True,
    min_count: int = 3,
    gene_col: str | None = None,
    method: str = "jaccard",
    obo_file: str | None = None,
    verbose: bool = True,
    gene_sets: Mapping[str, Iterable[Any]] | str | None = None,
    organism: str = "mouse",
    to_test_threshold: float = 0.0,
    pval_threshold: float = 0.05,
    term_size_limit: int = 0,
    show_excluded: bool = False,
) -> pd.DataFrame:
    """
    Redundancy reduction for enrichment results.

    Parameters
    ----------
    method : {"jaccard", "pathway_denester", "goatools"}, default "jaccard"
        - ``jaccard``: greedy filtering by Jaccard overlap of enriched gene sets.
        - ``pathway_denester``: combinatorial nested-pathway test from
          `PathwayDenester <https://github.com/Helmy-Lab/PathwayDenester>`_.
          Requires full pathway gene memberships via ``gene_sets`` (or resolvable
          from ``enrich_df.attrs['gene_set_info']``).
        - ``goatools``: not implemented.

    gene_sets : dict, GMT path, or bundled library name, optional
        Full pathway definitions used only for ``method="pathway_denester"``.
        If omitted, the function tries to reload the library recorded in
        ``enrich_df.attrs['gene_set_info']``.

    to_test_threshold : float, default 0.0
        PathwayDenester only. Minimum fraction of shared DEGs (relative to the
        smaller pathway) before testing nested enrichment.

    pval_threshold : float, default 0.05
        PathwayDenester only. Independence p-value cutoff for excluding a term.

    term_size_limit : int, default 0
        PathwayDenester only. Drop pathways larger than this size before testing.
        ``0`` or negative values keep all pathways.

    show_excluded : bool, default False
        PathwayDenester only. If True, return all terms with diagnostic
        ``Denester_*`` columns; otherwise return only kept terms.
    """
    if not (0 <= similarity_cutoff <= 1):
        raise ValueError("similarity_cutoff must be between 0 and 1")
    if min_count < 1:
        raise ValueError("min_count must be >= 1")
    if enrich_df is None or enrich_df.empty:
        return enrich_df

    method_norm = str(method).lower().replace("-", "_")
    if method_norm in {"denester", "pathwaydenester"}:
        method_norm = "pathway_denester"

    df = enrich_df.copy()
    if gene_col is None:
        for c in ["Genes", "Lead_genes", "leadingEdge", "leading_edge"]:
            if c in df.columns:
                gene_col = c
                break
    if gene_col is None or gene_col not in df.columns:
        if verbose:
            _log_info("No suitable gene column found. Returning original DataFrame.")
        return df
    if by is None:
        for c in ["p.adjust", "FDR_qval", "pvalue"]:
            if c in df.columns:
                by = c
                break
    if by and by in df.columns and method_norm == "jaccard":
        df = df.sort_values(by, ascending=ascending).reset_index(drop=True)
    size_col = "Count" if "Count" in df.columns else "Size"
    if size_col in df.columns:
        df = df[df[size_col] >= min_count].copy()
    if df.empty:
        return df

    cluster_col = None
    for cand in ("Cluster", "cluster", "group", "Group"):
        if cand in df.columns:
            cluster_col = cand
            break

    if method_norm == "jaccard":
        if cluster_col and df[cluster_col].nunique() > 1:
            # Simplify within each cluster to preserve multi-group structure (clusterProfiler-like)
            kept_parts = []
            for _cl, sub in df.groupby(cluster_col, sort=False):
                kept = []
                kept_sets: list[set[str]] = []
                for idx, row in sub.iterrows():
                    genes_str = str(row.get(gene_col, ""))
                    current = {g.strip() for g in re.split(r"[;,]+", genes_str) if g.strip()}
                    if not current:
                        continue
                    redundant = any(
                        (len(current & s) / len(current | s) if len(current | s) > 0 else 0)
                        >= similarity_cutoff
                        for s in kept_sets
                    )
                    if not redundant:
                        kept.append(idx)
                        kept_sets.append(current)
                if kept:
                    kept_parts.append(df.loc[kept])
            result = pd.concat(kept_parts, ignore_index=True) if kept_parts else df.iloc[0:0].copy()
            if verbose:
                _log_info(
                    f"Simplified (per-cluster Jaccard >= {similarity_cutoff}) from {len(df)} to {len(result)} terms"
                )
        else:
            kept = []
            kept_sets = []
            for idx, row in df.iterrows():
                genes_str = str(row.get(gene_col, ""))
                current = {g.strip() for g in re.split(r"[;,]+", genes_str) if g.strip()}
                if not current:
                    continue
                redundant = any(
                    (len(current & s) / len(current | s) if len(current | s) > 0 else 0)
                    >= similarity_cutoff
                    for s in kept_sets
                )
                if not redundant:
                    kept.append(idx)
                    kept_sets.append(current)
            if verbose:
                _log_info(
                    f"Simplified from {len(df)} to {len(kept)} terms (Jaccard >= {similarity_cutoff})"
                )
            result = df.loc[kept].reset_index(drop=True)
        if hasattr(enrich_df, "attrs"):
            result.attrs.update(dict(enrich_df.attrs))
        result.attrs["simplify_method"] = "jaccard"
        if cluster_col:
            result.attrs["simplify_per_cluster"] = bool(
                cluster_col and df[cluster_col].nunique() > 1
            )
        return result
    if method_norm == "pathway_denester":
        if by is None or by not in df.columns:
            raise ValueError(
                "method='pathway_denester' requires a p-value column (e.g. 'p.adjust' or 'pvalue')."
            )
        term_to_genes = _resolve_gene_sets_for_simplify(
            df, gene_sets=gene_sets, organism=organism, verbose=verbose
        )
        if cluster_col and df[cluster_col].nunique() > 1:
            parts = []
            for _, sub in df.groupby(cluster_col, sort=False):
                if sub.empty:
                    continue
                s = _simplify_by_pathway_denester(
                    sub,
                    term_to_genes,
                    gene_col=gene_col,
                    p_col=by,
                    to_test_threshold=to_test_threshold,
                    pval_threshold=pval_threshold,
                    term_size_limit=term_size_limit,
                    show_excluded=show_excluded,
                    verbose=verbose,
                )
                if not s.empty:
                    parts.append(s)
            result = pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0].copy()
            if hasattr(enrich_df, "attrs"):
                result.attrs.update(dict(enrich_df.attrs))
            result.attrs["simplify_method"] = "pathway_denester"
            result.attrs["simplify_per_cluster"] = True
            return result
        return _simplify_by_pathway_denester(
            df,
            term_to_genes,
            gene_col=gene_col,
            p_col=by,
            to_test_threshold=to_test_threshold,
            pval_threshold=pval_threshold,
            term_size_limit=term_size_limit,
            show_excluded=show_excluded,
            verbose=verbose,
        )
    if method_norm == "goatools":
        raise NotImplementedError(
            "method='goatools' is not implemented yet. "
            "Use method='jaccard' or method='pathway_denester'."
        )
    raise ValueError("method must be 'jaccard', 'pathway_denester', or 'goatools'")
