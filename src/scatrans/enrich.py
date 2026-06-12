import logging
import os
import re
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --- importlib.resources compatibility (py>=3.9 stdlib, else backport) ---
# Follows the same pattern used in pp_bias.py for robust package data access.
if sys.version_info >= (3, 9):
    from importlib.resources import as_file, files
else:
    from importlib_resources import as_file, files


@contextmanager
def _open_package_data(filename: str) -> Iterator[Path]:
    """Yield a real filesystem Path for a file inside scatrans/data/.

    Safe for wheels, sdists, and editable installs.
    """
    ref = files("scatrans.data") / filename
    with as_file(ref) as concrete:
        yield Path(concrete)


def _parse_gmt_content(
    content: str, gene_case: Optional[str] = None
) -> Tuple[Dict[str, set], Dict[str, str]]:
    """Parse GMT text content (term<TAB>desc<TAB>gene1<TAB>gene2...)."""
    term_to_genes: Dict[str, set] = {}
    term_to_desc: Dict[str, str] = {}
    for line in content.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        term = parts[0].strip()
        if not term:
            continue
        desc = parts[1].strip() if len(parts) > 1 else ""
        cleaned = _clean_gene_list(parts[2:], gene_case=gene_case)
        if not cleaned:
            continue
        if term in term_to_genes:
            term_to_genes[term].update(cleaned)
            if not term_to_desc.get(term):
                term_to_desc[term] = desc
        else:
            term_to_genes[term] = set(cleaned)
            term_to_desc[term] = desc
    return term_to_genes, term_to_desc


def list_bundled_gene_sets(verbose: bool = False) -> List[str]:
    """
    List gene set files (.gmt and similar) that are bundled inside the package
    under scatrans/data/.

    These are typically ClusterProfiler-derived GO/KEGG sets that you can use
    for better consistency with clusterProfiler results (same gene-term mappings).

    Users can pass the filename (with or without extension) to
    `run_enrichment(gene_sets=...)` or `run_kegg(...)`.

    Example:
        sets = scat.list_bundled_gene_sets()
        # Recommended: just use base names + organism (auto-resolves to the Hs/Mm 2026 built-ins)
        res = scat.run_kegg(genes, organism="mouse")
        res = scat.run_enrichment(genes, gene_sets="GO_Biological_Process", organism="mouse")
        # Legacy names are also accepted and mapped to the 2026 sets.
    """
    discovered: List[str] = []
    try:
        data_traversable = files("scatrans.data")
        for item in data_traversable.iterdir():
            name = getattr(item, "name", str(item))
            if name.endswith((".gmt", ".tsv", ".txt")) or "gene_set" in name.lower():
                discovered.append(name)
    except Exception:
        discovered = []

    # Fallback known names. New organism-specific defaults (2026) are the built-in
    # libraries when organism is given and no specific version is requested.
    known = [
        "Hs_GO_Biological_Process_2026.txt",
        "Hs_KEGG_2026.txt",
        "Mm_GO_Biological_Process_2026.txt",
        "Mm_KEGG_2026.txt",
        "GO_Biological_Process_scATrans.gmt",
        "KEGG_scATrans.gmt",
    ]
    files_list = sorted(set(discovered)) if discovered else known

    if verbose:
        logger.info("Available bundled gene sets in package data:")
        for f in files_list:
            logger.info(f"   • {f}")
        if not files_list:
            logger.info(
                "   (no extra .gmt files found yet — add your ClusterProfiler-derived sets to src/scatrans/data/)"
            )

    return files_list


def _try_load_bundled_gene_set(
    name: str, gene_case: Optional[str] = None
) -> Optional[Tuple[Dict[str, set], Dict[str, str]]]:
    """
    Try to load a gene set from the files bundled in scatrans/data/.

    Accepts:
      - exact filename present in data/   e.g. "KEGG_2026_cp.gmt"
      - basename without extension        e.g. "KEGG_2026_cp"
      - name that matches a file when .gmt is appended
    """
    candidates = [name]
    if not name.endswith((".gmt", ".tsv", ".txt")):
        candidates.append(name + ".txt")
        candidates.append(name + ".gmt")
        candidates.append(name + ".tsv")

    for cand in candidates:
        try:
            with _open_package_data(cand) as p:
                content = p.read_text(encoding="utf-8")
            return _parse_gmt_content(content, gene_case=gene_case)
        except Exception:
            continue

    return None


def _resolve_gene_set_name(requested: str, source: str, organism: str = "mouse") -> str:
    """
    Resolve gene set name.

    New default (per latest requirement):
    - If no specific version/year is requested, default to the organism-specific
      built-in libraries added in data/ (Hs_*/Mm_*_2026.txt for human/mouse).
    - User only needs to specify organism (for run_kegg especially).
    - If user writes a specific Enrichr-style name with year (e.g. GO_Biological_Process_2023
      or KEGG_2021), it is passed through (to gseapy/Enrichr) unless it matches a bundled name.
    - Explicit gene_set_source="enrichr" forces original Enrichr.

    The 4 new built-in files are the default when nothing more specific is given.
    """
    if source == "enrichr":
        return requested

    prefix = "Mm"
    o = str(organism).lower()
    if o in ("human", "hs", "hsa"):
        prefix = "Hs"
    elif o in ("mouse", "mm", "mmu"):
        prefix = "Mm"

    # New organism-specific 2026 defaults (the 4 files added to data/)
    # These become the default built-in library if user does not specify a year/version.
    # Only map names for which we actually ship bundled files (BP + KEGG for each organism).
    default_map = {
        # GO base names -> organism 2026 built-in (only BP is bundled)
        "GO_Biological_Process": f"{prefix}_GO_Biological_Process_2026",
        "GO_Biological_Process_2023": f"{prefix}_GO_Biological_Process_2026",
        "GO_Biological_Process_2026": f"{prefix}_GO_Biological_Process_2026",
        "GO_BP": f"{prefix}_GO_Biological_Process_2026",
        # KEGG
        "KEGG": f"{prefix}_KEGG_2026",
        "KEGG_2026": f"{prefix}_KEGG_2026",
        "KEGG_2021": f"{prefix}_KEGG_2026",  # map old popular names to new default
    }

    mapped = default_map.get(requested)
    if mapped is not None:
        return mapped

    # Legacy *_scATrans (and .gmt) names from older examples/docs.
    # Map them to the current organism's 2026 built-in for seamless backward compat
    # (old .gmt files are no longer shipped; the 2026 txt files are the supported bundled sets).
    if (
        "scATrans" in requested
        or "_scATrans" in requested.lower()
        or "_scatrans" in requested.lower()
    ):
        # strip extension and legacy suffix to get base
        base = requested
        for ext in (".gmt", ".txt", ".tsv"):
            if base.lower().endswith(ext):
                base = base[: -len(ext)]
        base = base.replace("_scATrans", "").replace("_scatrans", "").strip("._ ").lower()
        if base in (
            "go_biological_process",
            "go_bp",
            "go_biological_process_2023",
            "go_biological_process_2026",
        ):
            return f"{prefix}_GO_Biological_Process_2026"
        if "kegg" in base:
            return f"{prefix}_KEGG_2026"
        # unknown legacy: fall through (will likely fail bundled and gseapy with helpful list)

    # If the user wrote a very specific year that we have a bundled for, use it
    # (e.g. they can still ask for 2023 if the file exists, but we now default to 2026)
    # Otherwise fall through to the name as-is (for custom or Enrichr)

    # Final fallback: use as written (allows Enrichr historical or user .gmt)
    return requested


ORA_COLUMNS = [
    "Term",
    "Description",
    "Count",
    "GeneRatio",
    "GeneRatio_str",
    "BgRatio",
    "BgRatio_str",
    "FoldEnrichment",
    "RichFactor",
    "Overlap",
    "pvalue",
    "p.adjust",
    "neg_log10_padj",
    "Genes",
    "TermSize",
]


def _log_info(msg: str, verbose: bool = True) -> None:
    if verbose:
        logger.info(msg)


def _warn_user(msg: str) -> None:
    warnings.warn(msg, UserWarning, stacklevel=2)
    logger.warning(msg)


def _apply_gene_case(genes: Iterable[Any], gene_case: Optional[str]) -> List[str]:
    if gene_case is not None:
        gene_case = str(gene_case).lower()
    if gene_case is None:
        return [str(g).strip() for g in genes]
    if gene_case == "upper":
        return [str(g).strip().upper() for g in genes]
    if gene_case == "lower":
        return [str(g).strip().lower() for g in genes]
    raise ValueError("gene_case must be None, 'upper', or 'lower'")


def _clean_gene_list(
    gene_list: Optional[Iterable[Any]], gene_case: Optional[str] = None
) -> List[str]:
    if gene_list is None:
        return []
    cleaned = _apply_gene_case(gene_list, gene_case)
    s = pd.Series(cleaned)
    s = s.dropna().astype(str).str.strip()
    s = s[(s != "") & (s.str.lower() != "nan")]
    return s.drop_duplicates().tolist()


def _load_gene_sets(
    gene_sets_input: Union[Mapping[str, Iterable[Any]], str, None],
    organism: str = "mouse",
    verbose: bool = True,
    gene_case: Optional[str] = None,
) -> Tuple[Dict[str, set], Dict[str, str]]:
    if gene_sets_input is None:
        raise ValueError("gene_sets cannot be None")
    if isinstance(gene_sets_input, Mapping):
        term_to_genes: Dict[str, set] = {}
        term_to_desc: Dict[str, str] = {}
        for k, v in gene_sets_input.items():
            cleaned = _clean_gene_list(v, gene_case=gene_case)
            if cleaned:
                term = str(k).strip()
                if term in term_to_genes:
                    term_to_genes[term].update(cleaned)
                else:
                    term_to_genes[term] = set(cleaned)
                    term_to_desc[term] = ""
        return term_to_genes, term_to_desc
    if isinstance(gene_sets_input, str):
        looks_like_path = (
            os.path.exists(gene_sets_input)
            or os.path.isabs(gene_sets_input)
            or "/" in gene_sets_input
            or "\\" in gene_sets_input
        )
        if looks_like_path:
            if not os.path.exists(gene_sets_input):
                raise FileNotFoundError(f"GMT file not found: {gene_sets_input}")
            with open(gene_sets_input, encoding="utf-8") as f:
                content = f.read()
            term_to_genes, term_to_desc = _parse_gmt_content(content, gene_case=gene_case)
            return term_to_genes, term_to_desc

        # 1. Try bundled package data first (ClusterProfiler-derived GO/KEGG etc.)
        bundled = _try_load_bundled_gene_set(gene_sets_input, gene_case=gene_case)
        if bundled is not None:
            term_to_genes, term_to_desc = bundled
            if verbose:
                _log_info(f"Loaded bundled gene set '{gene_sets_input}' from package data")
            return term_to_genes, term_to_desc

        # 2. Fall back to gseapy / Enrichr (original behavior)
        try:
            import gseapy as gp

            try:
                gene_dict = gp.get_library(name=gene_sets_input, organism=organism)
            except Exception:
                gene_dict = gp.get_library(name=gene_sets_input)
            term_to_genes = {}
            term_to_desc = {}
            for term, genes in gene_dict.items():
                cleaned = _clean_gene_list(genes, gene_case=gene_case)
                if cleaned:
                    t = str(term).strip()
                    term_to_genes[t] = set(cleaned)
                    term_to_desc[t] = ""
            if verbose:
                _log_info(f"Loaded gene set library '{gene_sets_input}' via gseapy")
            return term_to_genes, term_to_desc
        except Exception as e:
            # Final attempt: maybe the user meant a bundled set but gseapy was tried first
            bundled2 = _try_load_bundled_gene_set(gene_sets_input, gene_case=gene_case)
            if bundled2 is not None:
                if verbose:
                    _log_info(
                        f"Loaded bundled gene set '{gene_sets_input}' from package data (after gseapy fallback)"
                    )
                return bundled2
            raise ValueError(
                f"Failed to load '{gene_sets_input}' via gseapy or as bundled package data: {e}\n"
                f"Available bundled sets: {list_bundled_gene_sets()}"
            ) from e

    raise ValueError(
        "gene_sets must be dict, GMT path or gseapy library name (or a bundled set name)"
    )


def _bh_p_adjust(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    if not np.isfinite(pvalues).all():
        raise ValueError("pvalues must contain only finite values.")
    n = len(pvalues)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(pvalues)
    sorted_p = pvalues[sorted_idx]
    adjusted = sorted_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    out = np.empty(n)
    out[sorted_idx] = adjusted
    return out


def run_enrichment(
    gene_list: Iterable[Any],
    gene_sets: Union[Mapping[str, Iterable[Any]], str],
    universe: Optional[Iterable[Any]] = None,
    background: Optional[Iterable[Any]] = None,
    adata: Optional[
        Any
    ] = None,  # NEW: if provided and no explicit universe, we try to use the preserved raw_gene_list
    pval_cutoff: float = 0.05,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    organism: str = "mouse",
    gene_case: Optional[str] = None,
    gene_set_source: str = "scatrans",
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
    and other diagnostics.

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
    """
    genes = _clean_gene_list(gene_list, gene_case=gene_case)
    if not genes:
        if verbose:
            _log_info("gene_list is empty")
        return pd.DataFrame(columns=ORA_COLUMNS)
    if min_size < 1 or max_size < min_size or not (0 <= pval_cutoff <= 1):
        raise ValueError("Invalid min_size, max_size or pval_cutoff")

    # Resolve gene set name based on explicit source choice (new clean API)
    if isinstance(gene_sets, str):
        gene_sets = _resolve_gene_set_name(gene_sets, gene_set_source, organism)

    term_to_genes, term_to_desc = _load_gene_sets(
        gene_sets, organism=organism, verbose=verbose, gene_case=gene_case
    )
    all_gs_genes = set().union(*term_to_genes.values()) if term_to_genes else set()

    # --- Universe / background resolution (clusterProfiler-aligned conservative default) ---
    # `universe` is now the preferred name; `background` kept for full backward compat.
    #
    # Smart default: if the user did not explicitly pass universe/background,
    # and they pass an `adata` on which `store_raw_counts` was previously called,
    # we automatically use the preserved full measured gene list.
    # This is much more robust than relying on adata.var_names after HVG subsetting.
    provided = universe if universe is not None else background

    if provided is None and adata is not None:
        try:
            if "scatrans" in adata.uns and "raw_gene_list" in adata.uns["scatrans"]:
                preserved = adata.uns["scatrans"]["raw_gene_list"]
                if preserved:
                    provided = preserved
                    if verbose:
                        _log_info(
                            "Using preserved raw_gene_list from adata.uns['scatrans'] "
                            f"({len(preserved)} genes) as universe (from previous store_raw_counts)."
                        )
        except Exception:
            pass  # be defensive

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

    universe = effective_universe  # used by the rest of the function and per-term calcs
    N = len(universe)

    # Rich diagnostics so users understand effective N (why it may be < provided background)
    provided_size = len(bg_set) if bg_set else (len(all_gs_genes) if provided_is_str_all else 0)
    dropped_by_restrict = provided_size - len(universe) if bg_set and not force_universe else 0
    restricted = bool(
        bg_set
        and not force_universe
        and restrict_background_to_gene_sets
        and dropped_by_restrict > 0
    )

    if N == 0:
        if verbose:
            _log_info("Universe is empty")
        return pd.DataFrame(columns=ORA_COLUMNS)
    genes_in_universe = [g for g in genes if g in universe]
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
        _warn_user(
            f"Low mapping rate ({mapping_rate:.1%}). Check gene ID type, organism and gene_case."
        )
    if n == 0:
        return pd.DataFrame(columns=ORA_COLUMNS)
    results = []
    for term, term_genes in term_to_genes.items():
        term_genes_in_universe = term_genes & universe
        K = len(term_genes_in_universe)
        if min_size > K or max_size < K:
            continue
        overlap = set(genes_in_universe) & term_genes_in_universe
        k = len(overlap)
        if k == 0:
            continue
        pval = hypergeom.sf(k - 1, N, K, n)
        GeneRatio = k / n
        BgRatio = K / N
        results.append(
            {
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
        )
    if not results:
        return pd.DataFrame(columns=ORA_COLUMNS)
    res_df = pd.DataFrame(results)
    res_df["p.adjust"] = _bh_p_adjust(res_df["pvalue"].values)
    res_df["neg_log10_padj"] = -np.log10(res_df["p.adjust"].astype(float).clip(lower=1e-300))
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
    res_df.attrs.update(
        {
            "method": "ora",
            "organism": organism,
            "gene_case": gene_case,
            "universe_info": universe_info,
            "clusterprofiler_aligned": True,
        }
    )
    if return_all:
        return res_df
    sig = res_df[res_df["p.adjust"] < pval_cutoff].copy().reset_index(drop=True)
    sig.attrs.update(res_df.attrs)
    if verbose:
        _log_info(f"Found {len(sig)} significant terms")
    return sig


def run_kegg(
    gene_list: Iterable[Any],
    organism: str = "mouse",
    universe: Optional[Iterable[Any]] = None,
    background: Optional[Iterable[Any]] = None,
    adata: Optional[Any] = None,  # forwarded to run_enrichment for smart universe default
    pval_cutoff: float = 0.05,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: Optional[str] = None,
    kegg_library: Optional[str] = None,
    gene_set_source: str = "scatrans",
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
    """
    org_map = {"mouse": "Mouse", "mmu": "Mouse", "human": "Human", "hsa": "Human"}
    gseapy_org = org_map.get(str(organism).lower())
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
        min_size=min_size,
        max_size=max_size,
        restrict_background_to_gene_sets=restrict_background_to_gene_sets,
        force_universe=force_universe,
        return_all=return_all,
        verbose=verbose,
        organism=gseapy_org,
        gene_case=gene_case,
        gene_set_source=gene_set_source,
        **kwargs,
    )


def simplify_enrichment(
    enrich_df: pd.DataFrame,
    similarity_cutoff: float = 0.5,
    by: Optional[str] = None,
    ascending: bool = True,
    min_count: int = 3,
    gene_col: Optional[str] = None,
    method: str = "jaccard",
    obo_file: Optional[str] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Greedy redundancy reduction using Jaccard gene overlap.
    """
    if not (0 <= similarity_cutoff <= 1):
        raise ValueError("similarity_cutoff must be between 0 and 1")
    if min_count < 1:
        raise ValueError("min_count must be >= 1")
    if enrich_df is None or enrich_df.empty:
        return enrich_df
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
    if by and by in df.columns:
        df = df.sort_values(by, ascending=ascending).reset_index(drop=True)
    size_col = "Count" if "Count" in df.columns else "Size"
    if size_col in df.columns:
        df = df[df[size_col] >= min_count].copy()
    if df.empty:
        return df
    if method == "jaccard":
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
        return df.loc[kept].reset_index(drop=True)
    elif method == "goatools":
        raise NotImplementedError(
            "goatools semantic simplification is not implemented in this version."
        )
    else:
        raise ValueError("method must be 'jaccard' or 'goatools'")
