"""scatrans.enrich.gsea — internal package module."""

from __future__ import annotations

import logging
import warnings
from collections.abc import Iterable, Mapping
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
    _check_gene_set_mapping_rate,
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

# Prefer *signed* ranking metrics for GSEA. Non-negative heuristics such as
# active_score must not be auto-selected: preranked GSEA needs bidirectional
# ranks so NES can be negative (depletion / down-enrichment).
_GSEA_SCORE_COLUMN_PRIORITY: tuple[str, ...] = (
    "logFC",
    "logfoldchanges",
    "log2FoldChange",
    "lfc",
    "t_stat",
    "t",
    "stat",
    "wald",
    "score",
    "NES",
)

# Known non-negative / one-sided scATrans columns — never auto-pick for GSEA.
_GSEA_UNSIGNED_SCORE_COLUMNS: frozenset[str] = frozenset(
    {
        "active_score",
        "active_score_pval",
        "active_score_fdr",
        "unspliced_excess_residual",
        "unspliced_excess_delta",
        "unspliced_excess_pval",
        "unspliced_excess_fdr",
        "velocity_residual",
        "velocity_delta_raw",
        "p_adj",
        "p_val",
        "padj",
        "pvalue",
        "total_us_counts",
        "gene_length",
        "intron_number",
        "effective_gamma",
    }
)


def _labels_look_like_gene_symbols(labels: pd.Index | pd.Series) -> bool:
    """True when index/labels look usable as gene IDs (symbols *or* Entrez-style numerics).

    Rejects only a default :class:`~pandas.RangeIndex` (row numbers, not genes).
    Purely numeric string IDs (Entrez) are accepted so DE tables indexed by
    Entrez IDs still coerce correctly when a score column is present.
    """
    if len(labels) == 0:
        return False
    # Non-default indexes (including all-numeric Entrez-like labels) are gene IDs.
    # RangeIndex is the default empty/positional index, not gene IDs.
    return not isinstance(labels, pd.RangeIndex)


def _is_known_unsigned_gsea_column(name: str) -> bool:
    low = str(name).lower()
    if low in {c.lower() for c in _GSEA_UNSIGNED_SCORE_COLUMNS}:
        return True
    # residual / fdr / p-value style names that should not rank GSEA by default
    return any(
        token in low
        for token in (
            "active_score",
            "residual",
            "fdr",
            "pval",
            "p_val",
            "p_adj",
            "padj",
            "pvalue",
        )
    )


def _warn_if_one_sided_gsea_ranking(scores: pd.Series, score_col: str) -> None:
    """Warn when the ranking metric cannot support bidirectional NES.

    Classic preranked GSEA needs signed scores (e.g. logFC). Non-negative
    metrics such as ``active_score`` make negative NES / depletion impossible.
    """
    vals = pd.to_numeric(scores, errors="coerce")
    vals = vals[np.isfinite(vals.to_numpy(dtype=float, na_value=np.nan))]
    if len(vals) < 5:
        return
    known_unsigned = _is_known_unsigned_gsea_column(score_col)
    arr = vals.to_numpy(dtype=float)
    all_nonneg = bool(np.all(arr >= 0.0))
    all_nonpos = bool(np.all(arr <= 0.0))
    if not (known_unsigned or all_nonneg or all_nonpos):
        return
    # All zeros is degenerate but not specifically a "unsigned active_score" path.
    if np.allclose(arr, 0.0):
        return
    _warn_user(
        f"run_gsea: ranking metric {score_col!r} appears one-sided / non-negative "
        f"(n={len(arr)}, min={float(np.min(arr)):.4g}, max={float(np.max(arr)):.4g}). "
        "Preranked GSEA requires a *signed* metric (e.g. logFC) so NES can be "
        "negative for depleted / down-regulated sets. Prefer score_column='logFC' "
        "or pass a signed pd.Series. Continuing with the selected metric; results "
        "will only capture one enrichment direction."
    )


def _pick_gsea_score_column(df: pd.DataFrame, *, prefer: str | None = None) -> str:
    if prefer is not None:
        if prefer not in df.columns:
            raise ValueError(
                f"score_column={prefer!r} not found in ranked_genes DataFrame columns: "
                f"{list(df.columns)}"
            )
        return prefer
    for col in _GSEA_SCORE_COLUMN_PRIORITY:
        if col in df.columns:
            return col
    numeric_cols = [
        c
        for c in df.columns
        if pd.api.types.is_numeric_dtype(df[c]) and not _is_known_unsigned_gsea_column(c)
    ]
    if not numeric_cols:
        # Last resort: any numeric column (may be unsigned — caller will warn).
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if len(numeric_cols) == 1:
        return numeric_cols[0]
    if len(numeric_cols) > 1:
        for col in numeric_cols:
            low = str(col).lower()
            if "logfc" in low or "log_fc" in low or low in {"lfc", "t_stat", "stat", "wald"}:
                return col
        # Prefer a column with both positive and negative values (signed).
        for col in numeric_cols:
            vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if len(vals) >= 5 and (vals > 0).any() and (vals < 0).any():
                return col
        return numeric_cols[0]
    raise ValueError(
        "Could not infer a numeric ranking column from ranked_genes DataFrame. "
        f"Columns: {list(df.columns)}. Pass score_column= explicitly (prefer a "
        "signed metric such as 'logFC'), or provide a pd.Series indexed by gene names."
    )


def _coerce_ranked_genes_dataframe(
    df: pd.DataFrame,
    *,
    score_column: str | None = None,
) -> pd.Series:
    """Normalize DE / active_score result tables to gene-indexed pd.Series."""
    if df.empty:
        raise ValueError("ranked_genes DataFrame is empty.")

    gene_name_cols = [c for c in ("gene", "genes", "names", "symbol", "Gene") if c in df.columns]
    score_col: str | None = None

    # Standard scATrans output: gene IDs in the index (not a default RangeIndex).
    if _labels_look_like_gene_symbols(df.index):
        score_col = _pick_gsea_score_column(df, prefer=score_column)
        if (
            score_column is None
            and "logFC" in df.columns
            and score_col == "logFC"
            and "active_score" in df.columns
        ):
            logger.info(
                "run_gsea: ranking by signed 'logFC' (default). "
                "Pass score_column= explicitly for another metric; "
                "active_score is non-negative and is not auto-selected."
            )
        scores = pd.to_numeric(df[score_col], errors="coerce")
        out = pd.Series(scores.values, index=df.index.astype(str))
        _warn_if_one_sided_gsea_ranking(out, score_col)
        return out

    # Explicit gene column + numeric score columns (e.g. CSV export with RangeIndex).
    if gene_name_cols:
        gene_col = gene_name_cols[0]
        genes = df[gene_col].astype(str)
        if score_column is not None:
            score_col = score_column
        elif len(gene_name_cols) == 1 and df.shape[1] == 2:
            other = [c for c in df.columns if c != gene_col][0]
            score_col = other
        else:
            score_col = _pick_gsea_score_column(df, prefer=None)
        scores = pd.to_numeric(df[score_col], errors="coerce")
        out = pd.Series(scores.values, index=genes)
        _warn_if_one_sided_gsea_ranking(out, score_col)
        return out

    # Legacy: column 0 = gene names, column 1 = scores (no meaningful index).
    if df.shape[1] >= 2:
        col0 = df.iloc[:, 0]
        col1 = df.iloc[:, 1]
        if (
            _labels_look_like_gene_symbols(col0)
            and pd.to_numeric(col1, errors="coerce").notna().any()
        ):
            out = pd.Series(
                pd.to_numeric(col1, errors="coerce").values,
                index=col0.astype(str).values,
            )
            _warn_if_one_sided_gsea_ranking(out, str(df.columns[1]))
            return out

    if df.shape[1] == 1 and _labels_look_like_gene_symbols(df.index):
        scores = pd.to_numeric(df.iloc[:, 0], errors="coerce")
        out = pd.Series(scores.values, index=df.index.astype(str))
        _warn_if_one_sided_gsea_ranking(out, str(df.columns[0]))
        return out

    raise ValueError(
        "ranked_genes DataFrame format not recognized. Expected one of:\n"
        "  - gene symbols as index + numeric score column (prefer signed logFC)\n"
        "  - columns ['gene', <score>] or legacy [gene_names, scores]\n"
        "  - pd.Series indexed by gene names (recommended)\n"
        f"Got index type {type(df.index).__name__}, columns={list(df.columns)}. "
        "Pass score_column= to select the ranking metric explicitly."
    )


def run_gsea(
    ranked_genes: pd.Series | Mapping[str, float] | Iterable[str] | pd.DataFrame,
    gene_sets: Mapping[str, Iterable[Any]] | str,
    min_size: int = 15,
    max_size: int = 500,
    nperm: int = 1000,
    organism: str = "mouse",
    gene_case: str | None = None,
    gene_set_source: str = "scatrans",
    verbose: bool = True,
    seed: int = 42,
    threads: int = 4,
    ascending: bool = False,
    weight: float | None = None,
    weighted_score_type: str | float | None = None,
    score_column: str | None = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Pre-ranked Gene Set Enrichment Analysis (GSEA) using gseapy.prerank.

    This implements the classic GSEA algorithm on a user-provided ranked gene list
    (e.g. logFC, t-statistic or custom score from active_score / DE results).
    It is the Python equivalent of clusterProfiler::GSEA / Broad GSEA Preranked.

    Parameters
    ----------
    ranked_genes : pd.Series, dict, DataFrame, or list-like

        - Preferred: pd.Series with gene names as index and numeric scores as values.
          Higher score = more "up" in target group (e.g. logFC).
          The function will sort internally if needed.

        - pd.DataFrame from ``active_score`` / ``differential_expression`` ``all_results``:
          gene symbols in the **index**, numeric score in a column (auto-prefers signed
          ``logFC`` / t-stat style columns). Non-negative metrics such as
          ``active_score`` are **not** auto-selected (GSEA needs signed ranks for
          bidirectional NES); pass ``score_column=`` to force them (emits a warning).

        - Legacy DataFrame: two columns ``[gene_names, scores]`` with a default RangeIndex.
        - dict: gene -> score
        - list of genes: treated as pre-sorted from high to low (scores assigned decreasing).

        Gene names will be cleaned according to gene_case.
    score_column : str, optional
        When ``ranked_genes`` is a DataFrame, which column holds the ranking metric.
        Defaults to signed columns (``logFC``, then t-stat-like names). Prefer a
        signed metric; non-negative columns trigger a warning.
    gene_sets : str, dict or list
        Same as run_enrichment: bundled name (e.g. "GO_Biological_Process"), GMT path,
        dict of term->genes, or Enrichr library name.
    min_size, max_size : int
        Minimum / maximum number of genes in a gene set to consider.
    nperm : int
        Number of permutations for p-value estimation.
    organism : str
        "mouse" or "human" (used for Enrichr/gseapy library lookup).
    gene_case : {"upper", "lower", None}
        Case normalization for gene symbols (same as other enrichment functions).
    gene_set_source : {"scatrans", "enrichr"}
        Control source preference (same semantics as run_enrichment).
    verbose : bool
        Print progress.
    seed : int
        Random seed forwarded to ``gseapy.prerank`` (reproducible permutations).
    threads : int
        CPU threads for gseapy prerank.
    ascending : bool
        If True, lower ranked metric = more enriched (gseapy convention).
    weight : float, optional
        GSEA enrichment weight passed to gseapy (Broad ``p`` exponent; default ``1.0`` = weighted).
    weighted_score_type : str or float, optional
        Deprecated alias for ``weight``. Broad/GSEA naming: ``"classic"`` → ``0.0`` (unweighted KS);
        ``"weighted"`` → ``1.0``. When omitted, defaults to weighted (``1.0``).
    **kwargs
        Additional arguments forwarded to ``gseapy.prerank`` (e.g. ``graph_num``).

    Returns
    -------
    pd.DataFrame
        GSEA results with columns including Term, Description, ES, NES, pvalue, p.adjust,
        neg_log10_padj, leading_edge, etc. Sorted by ``|NES|`` (absolute value) descending
        so that the strongest magnitude effects (positive or negative) appear first.
        Rich metadata in .attrs (method="gsea_prerank", gene_set_info, nperm, gsea_info, analysis_info).

    Notes
    -----

    - Unlike ORA, GSEA does not use an explicit "universe" in the same way; the ranked
      list itself defines the background. min_size/max_size still apply.

    - **Mapping check (same as ORA):** after loading gene sets, ranked genes are
      intersected with gene-set members. Mapping rate &lt; 20% emits a UserWarning
      with symbol examples; zero overlap returns an empty frame with
      ``reason="no_ranked_genes_mapped"`` (avoids opaque gseapy filter failures from
      case/ID mismatches). Prefer ``gene_case="upper"`` for Enrichr libraries.

    - Prefer **signed** ranking metrics (logFC). Non-negative columns such as
      ``active_score`` are not auto-selected and trigger a warning if forced.

    - Requires gseapy. Install via `pip install gseapy` or `pip install "scatrans[gsea]"`.
    - ``all_results`` from active_score can be passed directly::

        res = scat.run_gsea(all_results, gene_sets="GO_Biological_Process", score_column="logFC")

      or as a Series: ``all_results["logFC"]`` (index = gene names).
    """
    try:
        import gseapy as gp
    except ImportError as e:
        raise ImportError(
            "run_gsea requires the 'gseapy' package. "
            "Please install it with: pip install gseapy or pip install 'scatrans[gsea]'"
        ) from e

    analysis_info = _get_analysis_info()
    organism_norm = str(organism).lower()

    # Normalize ranked_genes input to pd.Series (gene -> score).
    # Flags avoid noisy warnings for synthetic ranks (gene list → n..1) and
    # double-warnings after DataFrame coercion (already warned there).
    _skip_unsigned_rank_warning = False
    if isinstance(ranked_genes, pd.DataFrame):
        ranked_genes = _coerce_ranked_genes_dataframe(ranked_genes, score_column=score_column)
        _skip_unsigned_rank_warning = True
    if isinstance(ranked_genes, (list, tuple)):
        # treat as pre-ordered high->low, assign descending ranks
        genes = _apply_gene_case([str(g).strip() for g in ranked_genes], gene_case)
        scores = list(range(len(genes), 0, -1))
        ranked = pd.Series(scores, index=genes)
        _skip_unsigned_rank_warning = True
    elif isinstance(ranked_genes, Mapping):
        ranked = pd.Series(ranked_genes)
    elif isinstance(ranked_genes, pd.Series):
        ranked = ranked_genes.copy()
    else:
        raise ValueError(
            "ranked_genes must be a pd.Series (gene->score), dict, or list of genes (sorted high->low)"
        )

    if not isinstance(ranked, pd.Series):
        raise ValueError(
            "ranked_genes must be a pd.Series indexed by gene names with numeric scores."
        )
    if ranked.index is None or len(ranked.index) == 0:
        raise ValueError("ranked_genes must have a non-empty gene index.")

    ranked.index = pd.Index(_apply_gene_case(ranked.index.astype(str).tolist(), gene_case))
    n_before_numeric = len(ranked)
    ranked = pd.to_numeric(ranked, errors="coerce")
    if n_before_numeric > 0 and ranked.notna().sum() == 0:
        raise ValueError(
            "ranked_genes scores must be numeric (gene → score pd.Series). "
            f"All {n_before_numeric} values became NaN after coercion."
        )
    ranked = ranked.dropna()
    if not _skip_unsigned_rank_warning:
        _warn_if_one_sided_gsea_ranking(ranked, score_column or "ranked_genes")
    # Duplicate IDs (case-folding / multi-mapped symbols): keep the score with
    # largest absolute value so ranking is not arbitrarily first-row dependent.
    if ranked.index.duplicated().any():
        n_dup = int(ranked.index.duplicated().sum())
        collapsed: dict[str, float] = {}
        for gene, val in ranked.items():
            v = float(val)
            key = str(gene)
            if key not in collapsed or abs(v) > abs(collapsed[key]):
                collapsed[key] = v
        ranked = pd.Series(collapsed, dtype=float)
        logger.info(
            "run_gsea: collapsed %d duplicate gene IDs keeping max |score| per gene.",
            n_dup,
        )
    if len(ranked) == 0:
        if verbose:
            _log_info("ranked_genes is empty after cleaning")
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="ranked_genes_empty",
            gene_set_info=None,
            analysis_info=analysis_info,
        )

    if min_size < 1 or max_size < min_size:
        raise ValueError("min_size and max_size must be positive with max_size >= min_size")

    requested_gene_sets = gene_sets if isinstance(gene_sets, str) else "<dict>"
    resolved_gene_sets = requested_gene_sets
    if isinstance(gene_sets, str):
        resolved_gene_sets = _resolve_gene_set_name(gene_sets, gene_set_source, organism_norm)
        gene_sets = resolved_gene_sets

    term_to_genes, term_to_desc, load_info = _load_gene_sets(
        gene_sets, organism=organism_norm, verbose=verbose, gene_case=gene_case
    )
    gene_set_info = {
        "requested": requested_gene_sets,
        "resolved": resolved_gene_sets,
        "requested_source": gene_set_source,
        "actual_source": load_info.get("actual_source"),
        "library_name": load_info.get("library_name"),
        "n_terms": int(len(term_to_genes)),
        "n_unique_genes": int(len(set().union(*term_to_genes.values()))) if term_to_genes else 0,
    }

    if not term_to_genes:
        if verbose:
            _log_info("No gene sets loaded")
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="no_gene_sets",
            gene_set_info=gene_set_info,
            analysis_info=analysis_info,
        )

    # Input sanity: same mapping-rate gate as ORA (species / case / ID mismatches).
    all_gs_genes: set[str] = set().union(*(set(gs) for gs in term_to_genes.values()))
    mapping_info = _check_gene_set_mapping_rate(
        ranked.index.astype(str).tolist(),
        all_gs_genes,
        context="run_gsea",
        threshold=0.2,
        gene_case=gene_case,
    )
    gene_set_info["mapping"] = mapping_info
    if verbose:
        _log_info(
            f"Ranked genes mapped to gene sets: {mapping_info['n_mapped']}/{mapping_info['n_input']} "
            f"({mapping_info['mapping_rate']:.1%})"
        )
    if mapping_info["n_mapped"] == 0:
        msg = (
            "run_gsea: zero ranked genes overlap any gene-set member "
            f"(n_ranked={mapping_info['n_input']}). "
            f"Input examples: {mapping_info['example_input']}; "
            f"gene-set examples: {mapping_info['example_reference']}. "
            "Check organism, gene ID type, and gene_case "
            "(Enrichr libraries are typically UPPERCASE — try gene_case='upper')."
        )
        _warn_user(msg)
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="no_ranked_genes_mapped",
            n_genes_ranked=int(mapping_info["n_input"]),
            n_genes_overlap=0,
            gene_set_info=gene_set_info,
            analysis_info=analysis_info,
        )

    # Prepare for gseapy: dict of str -> list
    gene_sets_for_gp = {term: list(genes) for term, genes in term_to_genes.items()}

    if verbose:
        _log_info(
            f"Running gseapy.prerank on {len(ranked)} genes with {len(gene_sets_for_gp)} gene sets "
            f"(min_size={min_size}, max_size={max_size}, nperm={nperm})"
        )

    gsea_weight = _resolve_gseapy_weight(
        weight=weight,
        weighted_score_type=weighted_score_type,
    )

    try:
        prerank_kwargs = dict(kwargs)
        if "weighted_score_type" in prerank_kwargs:
            warnings.warn(
                "weighted_score_type is deprecated for run_gsea; use weight instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            if weight is None:
                gsea_weight = _resolve_gseapy_weight(
                    weight=None,
                    weighted_score_type=prerank_kwargs.pop("weighted_score_type"),
                )
            else:
                prerank_kwargs.pop("weighted_score_type")
        if "weight" in prerank_kwargs:
            gsea_weight = float(prerank_kwargs.pop("weight"))
        pre_res = gp.prerank(
            rnk=ranked,
            gene_sets=gene_sets_for_gp,
            min_size=min_size,
            max_size=max_size,
            permutation_num=nperm,
            outdir=None,
            no_plot=True,
            verbose=False,  # we control logging
            seed=prerank_kwargs.pop("seed", seed),
            threads=prerank_kwargs.pop("threads", threads),
            ascending=prerank_kwargs.pop("ascending", ascending),
            weight=gsea_weight,
            **prerank_kwargs,
        )
        res_df = pre_res.res2d.copy()
    except Exception as e:
        overlap = int(mapping_info["n_mapped"])
        msg = (
            f"gseapy.prerank failed: {e} "
            f"(ranked={len(ranked)}, overlap_with_gene_sets={overlap}, "
            f"mapping_rate={mapping_info['mapping_rate']:.1%}, gene_case={gene_case!r}). "
            "If overlap is low, check gene_case/organism (Enrichr uses UPPERCASE symbols)."
        )
        if verbose:
            _log_info(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="gseapy_error",
            error=str(e),
            n_genes_ranked=int(len(ranked)),
            n_genes_overlap=overlap,
            gene_set_info=gene_set_info,
            analysis_info=analysis_info,
        )

    if res_df is None or res_df.empty:
        overlap = int(mapping_info["n_mapped"])
        msg = (
            "gseapy returned no results (all gene sets filtered out?). "
            f"Ranked genes={len(ranked)}, overlap with gene sets={overlap} "
            f"({mapping_info['mapping_rate']:.1%}). "
            "Check gene symbols/IDs match the library; try lowering min_size "
            "or gene_case='upper' for Enrichr libraries."
        )
        if verbose:
            _log_info(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="no_results_after_filters",
            n_genes_ranked=int(len(ranked)),
            n_genes_overlap=int(overlap),
            gene_set_info=gene_set_info,
            analysis_info=analysis_info,
        )

    # Normalize column names for scATrans consistency
    rename = {
        "NOM p-val": "pvalue",
        "FDR q-val": "p.adjust",
        "FWER p-val": "FWER_pval",
        "Lead_genes": "leading_edge",
        "Tag %": "Tag_percent",
        "Gene %": "Gene_percent",
    }
    res_df = res_df.rename(columns={k: v for k, v in rename.items() if k in res_df.columns})

    # Ensure standard names
    if "pvalue" not in res_df.columns and "pval" in res_df.columns:
        res_df = res_df.rename(columns={"pval": "pvalue"})
    if "p.adjust" not in res_df.columns and "padj" in res_df.columns:
        res_df = res_df.rename(columns={"padj": "p.adjust"})

    # gseapy's res2d frequently returns numeric columns as dtype "object" (a column
    # of plain Python floats, not strings, but not cast to a numpy numeric dtype
    # either). Downstream consumers (matplotlib scatter coloring in enrich_dotplot,
    # sort_values, filter_active_genes-style cutoffs) assume real numeric dtypes;
    # an object-dtype color array in particular makes matplotlib fall back to its
    # (very slow, effectively hanging) per-point color-spec parsing path instead of
    # numeric colormap normalization. Coerce explicitly so the dtype is never object.
    _numeric_gsea_cols = [
        "ES",
        "NES",
        "pvalue",
        "p.adjust",
        "FWER_pval",
        "Tag_percent",
        "Gene_percent",
    ]
    for _col in _numeric_gsea_cols:
        if _col in res_df.columns:
            res_df[_col] = pd.to_numeric(res_df[_col], errors="coerce")

    # Add neg_log10_padj for compatibility with filters/plots
    if "p.adjust" in res_df.columns:
        res_df["neg_log10_padj"] = -np.log10(
            pd.to_numeric(res_df["p.adjust"], errors="coerce").clip(lower=1e-300)
        )
    if "pvalue" in res_df.columns and "neg_log10_padj" not in res_df.columns:
        res_df["neg_log10_pval"] = -np.log10(
            pd.to_numeric(res_df["pvalue"], errors="coerce").clip(lower=1e-300)
        )

    # Add TermSize if possible (from original gene sets)
    if "TermSize" not in res_df.columns:
        term_sizes = {t: len(gs) for t, gs in term_to_genes.items()}
        res_df["TermSize"] = res_df["Term"].map(term_sizes).fillna(0).astype(int)

    # Sort by |NES| descending (strongest effects first, regardless of sign) if present.
    if "NES" in res_df.columns:
        res_df = res_df.sort_values("NES", ascending=False, key=abs, na_position="last")
    elif "p.adjust" in res_df.columns:
        res_df = res_df.sort_values("p.adjust", ascending=True)
    res_df = res_df.reset_index(drop=True)

    # Reorder columns preferring GSEA_COLUMNS
    col_order = [c for c in GSEA_COLUMNS if c in res_df.columns]
    other_cols = [c for c in res_df.columns if c not in col_order]
    res_df = res_df[col_order + other_cols]

    # Attach rich diagnostics (consistent with ORA)
    res_df.attrs.update(
        {
            "method": "gsea_prerank",
            "organism": organism_norm,
            "gene_case": gene_case,
            "gene_set_info": gene_set_info,
            "nperm": int(nperm),
            "analysis_info": analysis_info,
            "clusterprofiler_aligned": True,
        }
    )
    res_df.attrs["gsea_info"] = {
        "n_genes_ranked": int(len(ranked)),
        "n_genes_overlap": int(mapping_info["n_mapped"]),
        "mapping_rate": float(mapping_info["mapping_rate"]),
        "mapping_info": mapping_info,
        "score_min": float(ranked.min()),
        "score_max": float(ranked.max()),
        "score_median": float(ranked.median()),
    }
    # Store gseapy internals for accurate gseaplot (RES curve + hits per term).
    # Wrapped in _DeepcopyImmuneDict: see its docstring for why (large nested
    # payload + pandas' per-operation .attrs deepcopy is a serious perf trap).
    if hasattr(pre_res, "results"):
        res_df.attrs["gsea_details"] = _DeepcopyImmuneDict(pre_res.results)
    if hasattr(pre_res, "ranking"):
        res_df.attrs["ranking"] = _DeepcopyImmuneDict(pre_res.ranking.to_dict())

    if verbose:
        n_sig = (
            int((res_df.get("p.adjust", pd.Series(1)) < 0.05).sum())
            if "p.adjust" in res_df
            else len(res_df)
        )
        _log_info(
            f"GSEA completed: {len(res_df)} terms tested, {n_sig} with p.adjust < 0.05 (sorted by |NES|)"
        )

    return res_df
