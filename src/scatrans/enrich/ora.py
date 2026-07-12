"""scatrans.enrich.ora — internal package module."""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

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


def run_enrichment(
    gene_list: Iterable[Any],
    gene_sets: Mapping[str, Iterable[Any]] | str,
    universe: Iterable[Any] | None = None,
    background: Iterable[Any] | None = None,
    adata: Any
    | None = None,  # NEW: if provided and no explicit universe, we try to use the preserved raw_gene_list
    pval_cutoff: float | None = None,
    padj_cutoff: float | None = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    organism: str = "mouse",
    gene_case: str | None = None,
    gene_set_source: str = "scatrans",
    include_gene_list: bool = False,
    p_adjust_method: str = "fdr_bh",
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Hypergeometric over-representation analysis (clusterProfiler-style ORA).

    Background / universe handling is designed to be close to clusterProfiler's
    `enricher` / `enrichGO` default conservative behavior:
    - If you provide a `universe` (or `background` for compat), by default it is
      intersected with the genes that appear in the gene_sets (i.e. have annotation).
      This matches clusterProfiler's default (see issues #283/#636).
    - Set `restrict_background_to_gene_sets=False` or `force_universe=True` to
      use the user-provided list untouched (analogous to clusterProfiler's
      `options(enrichment_force_universe = TRUE)`).
    - If neither provided, universe = union of all genes present in the gene_sets
      (safe default, similar to clusterProfiler when no universe given).

    New smart default (recommended):
    - If you do not pass `universe` or `background`, but you pass an `adata` on which
      `scat.store_raw_counts(adata)` was previously called, `run_enrichment` will
      automatically use the preserved full measured gene list (`adata.uns["scatrans"]["raw_gene_list"]`)
      as the background. This is the safest and most convenient behavior for
      single-cell data.
    - Explicit `universe=...` or `background=...` always takes precedence.

    Historical note on `universe`:
    Passing `universe=adata.var_names.tolist()` after HVG subsetting is usually wrong.
    The background should be the genes that were actually measured in the experiment.

    Returned DataFrame is rich: clusterProfiler-compatible columns + RichFactor,
    string helpers, TermSize, neg_log10_padj, plus detailed `.attrs["universe_info"]`
    and other diagnostics (including `gene_set_info` and `reason` on empty results).

    gene_list : list-like, pd.Series, or pd.DataFrame
        Genes to test for over-representation. Besides plain lists, accepts DE / filter
        output tables with gene symbols in the **index** (``all_results``, ``significant``,
        ``filter_active_genes(...)``), or explicit ``gene`` / ``names`` columns.

    pval_cutoff / padj_cutoff : float
        Cutoff applied to **adjusted p-values** (`p.adjust` column), **NOT** raw p-values.
        IMPORTANT: Despite the name, pval_cutoff filters on the BH-adjusted p-value.
        - Preferred: use `padj_cutoff` explicitly.
        - `pval_cutoff` is deprecated for new code (warning emitted when used alone).
        Default 0.05. If both passed, padj_cutoff wins.

    gene_set_source : {"scatrans", "enrichr"}, default "scatrans"
        Explicit override for which family to use.
        - "scatrans" (default): Prefer the bundled scATrans / clusterProfiler-derived sets.
        - "enrichr": Force the original Enrichr/gseapy libraries.

        In most cases you do **not** need this parameter:
        - Default behavior uses the package's bundled sets (only organism needed for run_kegg).
        - To pick a specific Enrichr historical version, just write the exact name
          (e.g. gene_sets="GO_Biological_Process_2021" or kegg_library="KEGG_2021").
          Names containing year suffixes are automatically treated as Enrichr requests.
        The parameter is mainly for forcing one side when the auto-detection would
        choose differently.

    include_gene_list : bool, default False
        If True, include an additional "Genes_list" column containing Python lists
        of the overlapping genes (in addition to the semicolon-joined "Genes" string).
        Useful for in-memory Python workflows; "Genes" remains for CSV/export compat.

    p_adjust_method : {"fdr_bh", "bonferroni", "none"}, default "fdr_bh"
        Multiple-testing correction applied **across all tested terms** in this call
        (including when a custom dict mixes GO/KEGG/custom pathways). For multi-cluster
        comparisons use ``compare_enrichment(..., adjust_across_clusters=True)``.

    background : optional
        Deprecated alias of `universe`. Use `universe` instead.
        If both are provided, raises ValueError.
        In docs and examples we strongly prefer the term `universe`.
    """
    # Normalize organism early for consistency (attrs, resolver, etc.)
    organism_norm = str(organism).lower()

    # Reproducibility info for manuscript supplementary materials
    analysis_info = _get_analysis_info()

    # Resolve cutoff: prefer explicit padj_cutoff, fall back to pval_cutoff (legacy name).
    # Both are applied to the adjusted p-value column ("p.adjust").
    cutoff = _resolve_enrichment_padj_cutoff(pval_cutoff, padj_cutoff)

    genes = _clean_gene_list(gene_list, gene_case=gene_case)
    if not genes:
        if verbose:
            _log_info("gene_list is empty")
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            method="ora",
            organism=organism_norm,
            gene_case=gene_case,
            reason="gene_list_empty",
            gene_set_info=None,
            universe_info=None,
            pval_cutoff=cutoff,
            clusterprofiler_aligned=True,
            analysis_info=analysis_info,
        )

    # Early validation (clear errors, consistent with review feedback)
    if min_size < 1:
        raise ValueError("min_size must be >= 1")
    if max_size < min_size:
        raise ValueError("max_size must be >= min_size")
    if not (0 <= cutoff <= 1):
        raise ValueError("pval_cutoff / padj_cutoff must be between 0 and 1")
    if gene_case is not None and gene_case not in (None, "upper", "lower"):
        raise ValueError("gene_case must be None, 'upper' or 'lower'")

    # Record the originally requested gene_sets name (for attrs / reproducibility)
    requested_gene_sets = gene_sets if isinstance(gene_sets, str) else "<dict>"

    # Resolve gene set name based on explicit source choice (new clean API)
    resolved_gene_sets = requested_gene_sets
    if isinstance(gene_sets, str):
        resolved_gene_sets = _resolve_gene_set_name(gene_sets, gene_set_source, organism_norm)
        gene_sets = resolved_gene_sets

    term_to_genes, term_to_desc, load_info = _load_gene_sets(
        gene_sets, organism=organism_norm, verbose=verbose, gene_case=gene_case
    )
    all_gs_genes = set().union(*term_to_genes.values()) if term_to_genes else set()

    # --- gene_set_info for reproducibility and diagnostics (always populated) ---
    gene_set_info = {
        "requested": requested_gene_sets,
        "resolved": resolved_gene_sets,
        "requested_source": gene_set_source,
        "actual_source": load_info.get("actual_source"),
        "library_name": load_info.get("library_name"),
        "n_terms": int(len(term_to_genes)),
        "n_unique_genes": int(len(all_gs_genes)),
    }
    if load_info.get("provenance"):
        gene_set_info["provenance"] = load_info["provenance"]
    if load_info.get("resolved_name"):
        gene_set_info["resolved_name"] = load_info["resolved_name"]

    if not term_to_genes:
        _warn_user(
            "No valid gene sets loaded after resolving `gene_sets`. "
            "Check the gene_sets name/path/dict or gene_set_source."
        )
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            method="ora",
            organism=organism_norm,
            gene_case=gene_case,
            reason="no_gene_sets_loaded",
            gene_set_info=gene_set_info,
            universe_info=None,
            pval_cutoff=cutoff,
            clusterprofiler_aligned=True,
            analysis_info=analysis_info,
        )

    # --- Universe / background resolution (clusterProfiler-aligned conservative default) ---
    # `universe` is the preferred name; `background` kept only as deprecated alias.
    # If both are given, raise immediately.
    if universe is not None and background is not None:
        raise ValueError(
            "Please provide only one of `universe` or `background`, not both. Use `universe`."
        )
    if background is not None and universe is None:
        _warn_user(
            "`background` parameter is deprecated, please use `universe` instead. "
            "`background` will be removed in a future version."
        )
        universe = background

    # Smart default: if the user did not explicitly pass universe/background,
    # and they pass an `adata` on which `store_raw_counts` was previously called,
    # we automatically use the preserved full measured gene list.
    # This is much more robust than relying on adata.var_names after HVG subsetting.
    provided = universe if universe is not None else background

    if provided is None and adata is not None:
        try:
            if "scatrans" in adata.uns:
                sc_meta = adata.uns["scatrans"]
                # Prefer full pre-HVG universe when store_raw_counts was re-run after
                # gene subsetting (sticky raw_gene_list_full).
                preserved = sc_meta.get("raw_gene_list_full")
                source_key = "raw_gene_list_full"
                if preserved is None or (hasattr(preserved, "__len__") and len(preserved) == 0):
                    preserved = sc_meta.get("raw_gene_list")
                    source_key = "raw_gene_list"
                # `preserved` may be a plain list (fresh run) or a numpy array /
                # pandas Index (after a round-trip through .h5ad, where h5py/anndata
                # deserializes stored lists as arrays). `if preserved:` raises
                # ValueError ("truth value of an array... is ambiguous") for those,
                # which used to be swallowed below and silently fall back to the
                # full GO/KEGG gene universe instead of the real measured genes.
                # Use len() instead, which works for all of these container types.
                if preserved is not None and len(preserved) > 0:
                    provided = list(preserved)
                    if verbose:
                        _log_info(
                            f"Using preserved {source_key} from adata.uns['scatrans'] "
                            f"({len(provided)} genes) as universe (from previous store_raw_counts)."
                        )
        except Exception as _e:
            # This should now only trigger for genuinely malformed data; surface it
            # as a warning (not just INFO) since it silently changes the background
            # gene universe used for the hypergeometric test.
            logger.warning(
                "Could not read raw_gene_list from adata.uns['scatrans']; "
                f"falling back to the full gene-set universe as background. Error: {_e}"
            )
            # leave provided as-is

    provided_is_str_all = False
    bg_set: set = set()
    if provided is not None:
        if isinstance(provided, str) and provided.lower() == "all":
            provided_is_str_all = True
            effective_universe = all_gs_genes
        else:
            bg_set = set(_clean_gene_list(provided, gene_case=gene_case))
            if force_universe:
                effective_universe = bg_set
            else:
                effective_universe = (
                    bg_set & all_gs_genes if restrict_background_to_gene_sets else bg_set
                )
    else:
        effective_universe = all_gs_genes

    universe_set = effective_universe  # rename for clarity inside func
    N = len(universe_set)

    # Rich diagnostics so users understand effective N (why it may be < provided background)
    provided_size = len(bg_set) if bg_set else (len(all_gs_genes) if provided_is_str_all else 0)
    dropped_by_restrict = provided_size - len(universe_set) if bg_set and not force_universe else 0
    restricted = bool(
        bg_set
        and not force_universe
        and restrict_background_to_gene_sets
        and dropped_by_restrict > 0
    )

    if N == 0:
        if verbose:
            _log_info(
                "Universe is empty after intersecting user background with annotated genes. "
                "Try restrict_background_to_gene_sets=False (or force_universe=True) if this is unexpected. "
                "Also check that your background contains the genes you are testing."
            )
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            method="ora",
            organism=organism_norm,
            gene_case=gene_case,
            reason="universe_empty",
            universe_info={
                "provided_size": int(provided_size),
                "gene_sets_genes": int(len(all_gs_genes)),
                "effective_universe_size": 0,
            },
            gene_set_info=gene_set_info,
            pval_cutoff=cutoff,
            clusterprofiler_aligned=True,
            analysis_info=analysis_info,
        )

    genes_in_universe = [g for g in genes if g in universe_set]
    n = len(genes_in_universe)
    if verbose:
        _log_info(f"Input genes: {len(genes)}, mapped to universe: {n}, effective universe: {N}")
        if bg_set or provided_is_str_all:
            _log_info(
                f"  Background provided size: {provided_size}, gene_sets total genes: {len(all_gs_genes)}, "
                f"effective (after intersect): {N} (restricted={restricted}, dropped_by_no_annotation={dropped_by_restrict})"
            )
        if force_universe:
            _log_info(
                "  force_universe=True → using raw user background (no forced intersect with gene sets)"
            )
    mapping_rate = n / max(len(genes), 1)
    if mapping_rate < 0.2:
        example_input = genes[:5]
        example_gs = list(all_gs_genes)[:5]
        _warn_user(
            f"Low mapping rate ({mapping_rate:.1%}). "
            f"Input examples: {example_input}; gene set examples: {example_gs}. "
            "Check gene ID type, organism and gene_case."
        )
    if n == 0:
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            method="ora",
            organism=organism_norm,
            gene_case=gene_case,
            reason="no_input_genes_mapped",
            universe_info={
                "provided_size": int(provided_size),
                "gene_sets_genes": int(len(all_gs_genes)),
                "effective_universe_size": int(N),
                "n_input_mapped": 0,
                "n_input_raw": int(len(genes)),
            },
            gene_set_info=gene_set_info,
            pval_cutoff=cutoff,
            clusterprofiler_aligned=True,
            analysis_info=analysis_info,
        )

    results = []
    n_terms_tested = 0  # size-eligible terms (includes zero-overlap)
    n_terms_size_excluded = 0  # failed min_size/max_size
    for term, term_genes in term_to_genes.items():
        term_genes_in_universe = term_genes & universe_set
        K = len(term_genes_in_universe)
        if min_size > K or max_size < K:
            n_terms_size_excluded += 1
            continue
        n_terms_tested += 1
        overlap = set(genes_in_universe) & term_genes_in_universe
        k = len(overlap)
        # Include zero-overlap terms (p=1) so BH/FDR denominator m equals all size-eligible
        # gene sets tested — clusterProfiler / standard ORA convention.
        pval = float(hypergeom.sf(k - 1, N, K, n)) if k > 0 else 1.0
        GeneRatio = k / n if n > 0 else 0.0
        BgRatio = K / N if N > 0 else 0.0
        row = {
            "Term": term,
            "Description": term_to_desc.get(term, ""),
            "Count": k,
            "GeneRatio": GeneRatio,
            "GeneRatio_str": f"{k}/{n}",
            "BgRatio": BgRatio,
            "BgRatio_str": f"{K}/{N}",
            "FoldEnrichment": GeneRatio / BgRatio if BgRatio > 0 else 0,
            "RichFactor": k / K if K > 0 else 0,
            "Overlap": f"{k}/{K}",
            "pvalue": pval,
            "Genes": ";".join(sorted(overlap)),
            "TermSize": K,
        }
        if include_gene_list:
            row["Genes_list"] = sorted(overlap)
        results.append(row)

    # Build result or empty with rich attrs
    base_attrs = {
        "method": "ora",
        "organism": organism_norm,
        "gene_case": gene_case,
        "gene_set_info": gene_set_info,
        "clusterprofiler_aligned": True,
        "analysis_info": analysis_info,
    }

    if not results:
        if verbose:
            _log_info(
                f"No terms passed min_size={min_size}/max_size={max_size} "
                f"(of {len(term_to_genes)} gene sets loaded)"
            )
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            reason="no_terms_after_size_filters",
            n_tested_terms=0,
            n_terms_tested=0,
            n_terms_size_excluded=int(n_terms_size_excluded),
            # legacy alias (was misnamed; historically counted *tested* terms)
            n_terms_size_filtered=0,
            pval_cutoff=cutoff,
            universe_info=base_attrs.get("universe_info"),  # may be None here but ok
            **base_attrs,
        )

    res_df = pd.DataFrame(results)
    res_df["p.adjust"] = _apply_p_adjust(res_df["pvalue"].values, method=p_adjust_method)
    res_df["neg_log10_padj"] = -np.log10(
        res_df["p.adjust"].astype(float).clip(lower=1e-300)
    )  # p.adjust values below 1e-300 are clipped to avoid -inf in neg_log10_padj
    res_df = res_df.sort_values("p.adjust").reset_index(drop=True)
    # Reorder to match declared ORA_COLUMNS for consistency (new columns included)
    col_order = [c for c in ORA_COLUMNS if c in res_df.columns]
    res_df = res_df[col_order + [c for c in res_df.columns if c not in col_order]]

    # Enriched diagnostics in attrs (universe_info + clusterProfiler alignment note)
    universe_info = {
        "provided_size": int(provided_size),
        "gene_sets_genes": int(len(all_gs_genes)),
        "effective_universe_size": int(N),
        "restricted_to_gene_sets": bool(restricted),
        "dropped_by_annotation_filter": int(dropped_by_restrict),
        "force_universe": bool(force_universe),
        "n_input_mapped": int(n),
        "n_input_raw": int(len(genes)),
    }
    res_df.attrs.update(base_attrs)
    res_df.attrs["universe_info"] = universe_info
    res_df.attrs["pval_cutoff"] = cutoff  # record the effective cutoff used
    res_df.attrs["p_adjust_method"] = p_adjust_method
    res_df.attrs["multiple_testing"] = {
        "scope": "all_tested_terms",
        "method": p_adjust_method,
        "n_tests": int(len(res_df)),
        "n_terms_tested": int(n_terms_tested),
        "n_terms_size_excluded": int(n_terms_size_excluded),
        # legacy key: previously misnamed; value is n size-eligible (tested) terms
        "n_terms_size_filtered": int(n_terms_tested),
        "includes_zero_overlap_terms": True,
    }

    if return_all:
        if verbose:
            _log_info(
                f"Tested {len(res_df)} terms (incl. zero-overlap); returning all (return_all=True)"
            )
        return res_df

    sig = res_df[res_df["p.adjust"] < cutoff].copy().reset_index(drop=True)
    sig.attrs.update(res_df.attrs)
    if verbose:
        _log_info(
            f"Tested {len(res_df)} terms (incl. zero-overlap); "
            f"found {len(sig)} significant at padj < {cutoff}"
        )
    return sig


def run_kegg(
    gene_list: Iterable[Any],
    organism: str = "mouse",
    universe: Iterable[Any] | None = None,
    background: Iterable[Any] | None = None,
    adata: Any | None = None,  # forwarded to run_enrichment for smart universe default
    pval_cutoff: float | None = None,
    padj_cutoff: float | None = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: str | None = None,
    kegg_library: str | None = None,
    gene_set_source: str = "scatrans",
    include_gene_list: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    KEGG pathway enrichment (wrapper around run_enrichment).

    Defaults to the bundled scATrans (clusterProfiler-derived) gene set.
    You only need to specify the organism.

    To use a specific original Enrichr version, just pass the full name:
        kegg_library="KEGG_2021"   # or KEGG_2019, KEGG_2016, etc.

    `gene_set_source` can be used as an explicit override ("scatrans" or "enrichr")
    if needed.

    Supports the same universe/background controls as run_enrichment for
    clusterProfiler-like conservative behavior by default.

    Note: internal organism is normalized to lowercase (e.g. "mouse") for attrs.
    """
    org_lower = str(organism).lower()
    org_map = {
        "mouse": "Mouse",
        "mm": "Mouse",
        "mmu": "Mouse",
        "human": "Human",
        "hs": "Human",
        "hsa": "Human",
    }
    gseapy_org = org_map.get(org_lower)
    if gseapy_org is None:
        raise ValueError(f"Unsupported organism '{organism}' for run_kegg")

    # Default is now the organism-specific built-in library (Hs_KEGG_2026 / Mm_KEGG_2026)
    # added to data/. User only needs organism.
    # Specific historical Enrichr names (e.g. "KEGG_2021") will be resolved accordingly.
    if kegg_library is None:
        kegg_library = (
            "KEGG"  # resolver will turn this into the correct Hs/Mm_2026 based on organism
        )

    return run_enrichment(
        gene_list=gene_list,
        gene_sets=kegg_library,
        universe=universe,
        background=background,
        adata=adata,
        pval_cutoff=pval_cutoff,
        padj_cutoff=padj_cutoff,
        min_size=min_size,
        max_size=max_size,
        restrict_background_to_gene_sets=restrict_background_to_gene_sets,
        force_universe=force_universe,
        return_all=return_all,
        verbose=verbose,
        organism=gseapy_org,  # run_enrichment will normalize again to lower for attrs
        gene_case=gene_case,
        gene_set_source=gene_set_source,
        include_gene_list=include_gene_list,
        **kwargs,
    )


def run_go(
    gene_list: Iterable[Any],
    ontology: str = "BP",
    organism: str = "mouse",
    universe: Iterable[Any] | None = None,
    background: Iterable[Any] | None = None,
    adata: Any | None = None,
    pval_cutoff: float | None = None,
    padj_cutoff: float | None = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: str | None = None,
    gene_set_source: str = "scatrans",
    include_gene_list: bool = False,
    adjust_across_all: bool = False,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    GO enrichment (BP / CC / MF) wrapper around run_enrichment.

    Similar to clusterProfiler::enrichGO(ont = "BP").

    Parameters
    ----------
    ontology : {"BP", "CC", "MF", "ALL"}, default "BP"
        Which GO subtree to test.
        - "BP": Biological Process (default; also the only one with bundled 2026 sets)
        - "CC": Cellular Component
        - "MF": Molecular Function
        - "ALL": run BP + CC + MF and concatenate results (adds an "Ontology" column)

    adjust_across_all : bool, default False
        Only relevant for ontology="ALL".
        - False (default): p.adjust is the one computed *within each ontology* separately
          (each sub-run uses its own ``p_adjust_method``). Documented behavior for most users.
        - True: After concatenating BP+CC+MF, recompute a single correction across
          *all* tested GO terms using their raw pvalues and the same
          ``p_adjust_method`` passed through to :func:`run_enrichment` (default
          ``fdr_bh``; also honors ``bonferroni`` / ``none``). Recommended when you
          want one multiple-testing control across the entire GO.

    Only `organism` + `ontology` are usually required when using the built-in
    libraries. Other parameters (universe/background/adata, cutoffs, etc.)
    are forwarded to run_enrichment.

    When ontology="ALL", results are concatenated and an "Ontology" column
    ("BP"/"CC"/"MF") is prepended for disambiguation.

    The returned DataFrame.attrs contains rich diagnostics:
      - method, ontology="ALL", organism, gene_case, pval_cutoff, adjust_across_all
      - per_ontology_attrs: dict mapping "BP"/"CC"/"MF" -> the full .attrs from each sub-run
        (includes each ontology's own gene_set_info, universe_info, analysis_info etc.)
      - analysis_info (package/timestamp/version)
      - If adjust_across_all=True, the main "p.adjust" is the cross-GO re-adjustment
        using ``p_adjust_method``; the column "p.adjust.within_ontology" preserves
        the original per-ontology adjusted values.

    Note on p.adjust for "ALL":
        For ontology="ALL", p.adjust is computed within each ontology separately
        before concatenation (unless adjust_across_all=True). Use per_ontology_attrs
        or the within_ontology column for full transparency. The cross-ontology step
        never silently falls back to BH when another ``p_adjust_method`` was requested.

    Examples
    --------
    >>> # Basic BP with bundled sets
    >>> res = run_go(markers, ontology="BP", organism="mouse", return_all=True)
    >>> # ALL ontologies + cross-ontology correction + full provenance in attrs
    >>> res = run_go(markers, ontology="ALL", organism="mouse",
    ...              universe=background, return_all=True, adjust_across_all=True)
    >>> print(res.attrs["per_ontology_attrs"]["BP"]["actual_source"])  # or gene_set_info etc.
    """
    analysis_info = _get_analysis_info()
    ont = str(ontology).upper().strip()
    if ont not in ("BP", "GO_BP"):
        logger.warning(
            "run_go(ontology=%s): only BP has fully bundled offline 2026 data; "
            "CC/MF will use gseapy/Enrichr (possible network access).",
            ontology,
        )
    ont_map = {
        "BP": "GO_Biological_Process",
        "GO_BP": "GO_Biological_Process",
        "CC": "GO_Cellular_Component",
        "GO_CC": "GO_Cellular_Component",
        "MF": "GO_Molecular_Function",
        "GO_MF": "GO_Molecular_Function",
    }

    if ont == "ALL":
        ont_list = ["BP", "CC", "MF"]
    else:
        if ont not in ont_map and ont not in {
            "GO_BIOLOGICAL_PROCESS",
            "GO_CELLULAR_COMPONENT",
            "GO_MOLECULAR_FUNCTION",
        }:
            # allow passing the full Enrichr-style name directly
            gs_name = ont
        else:
            gs_name = ont_map.get(ont, ont)

        if ont not in ("BP", "GO_BP"):
            # Only BP has bundled offline data. CC/MF will fall back to gseapy/Enrichr (needs optional dep + network).
            logger.info(
                "run_go(ontology=%s): CC/MF are not bundled; will use gseapy/Enrichr (requires internet if not cached).",
                ont,
            )

        return run_enrichment(
            gene_list=gene_list,
            gene_sets=gs_name,
            universe=universe,
            background=background,
            adata=adata,
            pval_cutoff=pval_cutoff,
            padj_cutoff=padj_cutoff,
            min_size=min_size,
            max_size=max_size,
            restrict_background_to_gene_sets=restrict_background_to_gene_sets,
            force_universe=force_universe,
            return_all=return_all,
            verbose=verbose,
            organism=organism,
            gene_case=gene_case,
            gene_set_source=gene_set_source,
            include_gene_list=include_gene_list,
            **kwargs,
        )  # adjust_across_all only affects the ALL branch below

    # "ALL" case: run three and concat
    frames = []
    per_ontology_attrs: dict[str, dict] = {}
    eff_cut = _resolve_enrichment_padj_cutoff(pval_cutoff, padj_cutoff)
    # Honor the same p_adjust_method as per-ontology run_enrichment calls (via kwargs).
    p_adjust_method = str(kwargs.get("p_adjust_method", "fdr_bh"))

    for o in ont_list:
        gs_name = ont_map[o]
        df = None
        try:
            df = run_enrichment(
                gene_list=gene_list,
                gene_sets=gs_name,
                universe=universe,
                background=background,
                adata=adata,
                pval_cutoff=pval_cutoff,
                padj_cutoff=padj_cutoff,
                min_size=min_size,
                max_size=max_size,
                restrict_background_to_gene_sets=restrict_background_to_gene_sets,
                force_universe=force_universe,
                return_all=True,  # we will filter at the end if needed
                verbose=verbose,
                organism=organism,
                gene_case=gene_case,
                gene_set_source=gene_set_source,
                include_gene_list=include_gene_list,
                **kwargs,
            )
        except Exception as e:
            if verbose:
                _log_info(f"GO {o} skipped/failed: {e}")
            per_ontology_attrs[o] = {"error": str(e)}
            continue

        # Always record attrs for this ontology (even if empty result)
        if df is not None:
            per_ontology_attrs[o] = dict(df.attrs) if hasattr(df, "attrs") else {}

        if df is not None and not df.empty:
            df = df.copy()
            df.insert(0, "Ontology", o)
            frames.append(df)

    if not frames:
        # return a properly attributed empty
        empty = _empty_ora_result(
            include_gene_list=include_gene_list,
            method="ora_go_all",
            ontology="ALL",
            organism=str(organism).lower(),
            gene_case=gene_case,
            reason="no_go_terms_for_all",
            gene_set_info=None,
            universe_info=None,
            pval_cutoff=eff_cut,
            clusterprofiler_aligned=True,
            analysis_info=analysis_info,
            adjust_across_all=bool(adjust_across_all),
            per_ontology_attrs=per_ontology_attrs,
        )
        return empty

    combined = pd.concat(frames, ignore_index=True, sort=False)
    # Try to keep a stable column order
    preferred = ["Ontology"] + [c for c in ORA_COLUMNS if c in combined.columns]
    other = [c for c in combined.columns if c not in preferred]
    combined = combined[preferred + other]

    # Build rich combined attrs for ALL (preserve per-ontology diagnostics)
    combined.attrs["method"] = "ora_go_all"
    combined.attrs["ontology"] = "ALL"
    combined.attrs["organism"] = str(organism).lower()
    combined.attrs["gene_case"] = gene_case
    combined.attrs["pval_cutoff"] = eff_cut
    combined.attrs["adjust_across_all"] = bool(adjust_across_all)
    combined.attrs["p_adjust_method"] = p_adjust_method
    combined.attrs["per_ontology_attrs"] = per_ontology_attrs
    combined.attrs["analysis_info"] = analysis_info
    # Also carry forward gene_set_info etc from first if useful, but per_ontology is authoritative
    if frames:
        first_attrs = frames[0].attrs if hasattr(frames[0], "attrs") else {}
        for k in ("gene_set_info", "universe_info", "clusterprofiler_aligned"):
            if k in first_attrs and k not in combined.attrs:
                combined.attrs[k] = first_attrs[k]

    if adjust_across_all and not combined.empty and "pvalue" in combined.columns:
        # Preserve original within-ontology corrected values for transparency.
        # Use the same method as per-ontology runs (do not hardcode BH).
        combined = combined.copy()
        combined["p.adjust.within_ontology"] = combined["p.adjust"]
        combined["p.adjust"] = _apply_p_adjust(combined["pvalue"].values, method=p_adjust_method)
        combined["neg_log10_padj"] = -np.log10(
            combined["p.adjust"].astype(float).clip(lower=1e-300)
        )
        combined = combined.sort_values("p.adjust").reset_index(drop=True)
        if verbose:
            _log_info(
                "GO ALL: re-adjusted p.adjust across BP+CC+MF "
                f"(adjust_across_all=True, p_adjust_method={p_adjust_method!r}); "
                "original per-ontology p.adjust saved in 'p.adjust.within_ontology'"
            )

    if return_all:
        if verbose:
            _log_info(f"GO ALL: combined {len(combined)} terms (BP+CC+MF)")
        return combined

    # Apply cutoff to the combined table (using the effective cutoff)
    sig = combined[combined["p.adjust"] < eff_cut].copy().reset_index(drop=True)
    sig.attrs.update(combined.attrs)
    if verbose:
        _log_info(
            f"GO ALL: tested/combined {len(combined)} terms; found {len(sig)} significant at padj < {eff_cut}"
        )
    return sig
