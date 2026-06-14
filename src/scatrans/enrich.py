import json
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

# --- importlib.resources compatibility (py>=3.10 stdlib, else backport)
# Using the backport on 3.9 avoids spec.origin=None issues in some install scenarios (editable/wheel in CI).
# Follows the same pattern used in pp_bias.py for robust package data access.
if sys.version_info >= (3, 10):
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
    """Parse GMT text content (term<TAB>desc<TAB>gene1<TAB>gene2...).

    The second column (description) is retained when present and stored
    under the "Description" output column. Many bundled sets ship with an
    empty description column (two tabs); this is handled gracefully.
    """
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


def _get_analysis_info() -> Dict[str, Any]:
    """Return reproducibility metadata for attrs / save_enrichment_report."""
    from datetime import datetime as _dt

    try:
        import scatrans as _scat

        _pkg_version = getattr(_scat, "__version__", None)
    except Exception:
        _pkg_version = None
    return {
        "package": "scatrans",
        "package_version": _pkg_version,
        "timestamp": _dt.now().isoformat(),
        "module": "scatrans.enrich",
    }


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
    s = pd.Series(list(gene_list))
    s = s.dropna()
    s = s.astype(str).str.strip()
    s = s[(s != "") & (s.str.lower() != "nan")]
    cleaned = _apply_gene_case(s.tolist(), gene_case)
    return pd.Series(cleaned).drop_duplicates().tolist()


def _load_gene_sets(
    gene_sets_input: Union[Mapping[str, Iterable[Any]], str, None],
    organism: str = "mouse",
    verbose: bool = True,
    gene_case: Optional[str] = None,
) -> Tuple[Dict[str, set], Dict[str, str], Dict[str, Any]]:
    """Load gene sets and return provenance info for reproducibility.

    Returns
    -------
    term_to_genes, term_to_desc, load_info
        load_info contains:
            actual_source: "dict" | "gmt" | "bundled" | "gseapy"
            library_name: basename (for files) or "<dict>"
            resolved_name: original/resolved input
            path: abspath for direct GMT files, else None
    """
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
        load_info = {
            "actual_source": "dict",
            "library_name": "<dict>",
            "resolved_name": "<dict>",
            "path": None,
        }
        return term_to_genes, term_to_desc, load_info
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
            load_info = {
                "actual_source": "gmt",
                "library_name": os.path.basename(gene_sets_input),
                "resolved_name": gene_sets_input,
                "path": os.path.abspath(gene_sets_input),
            }
            return term_to_genes, term_to_desc, load_info

        # 1. Try bundled package data first (ClusterProfiler-derived GO/KEGG etc.)
        bundled = _try_load_bundled_gene_set(gene_sets_input, gene_case=gene_case)
        if bundled is not None:
            term_to_genes, term_to_desc = bundled
            if verbose:
                _log_info(f"Loaded bundled gene set '{gene_sets_input}' from package data")
            load_info = {
                "actual_source": "bundled",
                "library_name": gene_sets_input,
                "resolved_name": gene_sets_input,
                "path": None,
            }
            return term_to_genes, term_to_desc, load_info

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
            load_info = {
                "actual_source": "gseapy",
                "library_name": gene_sets_input,
                "resolved_name": gene_sets_input,
                "path": None,
            }
            return term_to_genes, term_to_desc, load_info
        except Exception as e:
            # Final attempt: maybe the user meant a bundled set but gseapy was tried first
            bundled2 = _try_load_bundled_gene_set(gene_sets_input, gene_case=gene_case)
            if bundled2 is not None:
                if verbose:
                    _log_info(
                        f"Loaded bundled gene set '{gene_sets_input}' from package data (after gseapy fallback)"
                    )
                load_info = {
                    "actual_source": "bundled",
                    "library_name": gene_sets_input,
                    "resolved_name": gene_sets_input,
                    "path": None,
                }
                b_terms, b_desc = bundled2 if isinstance(bundled2, tuple) else (bundled2, {})
                return b_terms, b_desc, load_info
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


def _empty_ora_result(include_gene_list: bool = False, **attrs) -> pd.DataFrame:
    """Return an empty ORA result DataFrame with ORA_COLUMNS (and optional Genes_list) and .attrs populated.

    Ensures diagnostics are available even on empty results.
    """
    cols = list(ORA_COLUMNS)
    if include_gene_list and "Genes_list" not in cols:
        cols = cols + ["Genes_list"]
    df = pd.DataFrame(columns=cols)
    df.attrs.update(attrs)
    return df


def run_enrichment(
    gene_list: Iterable[Any],
    gene_sets: Union[Mapping[str, Iterable[Any]], str],
    universe: Optional[Iterable[Any]] = None,
    background: Optional[Iterable[Any]] = None,
    adata: Optional[
        Any
    ] = None,  # NEW: if provided and no explicit universe, we try to use the preserved raw_gene_list
    pval_cutoff: float = 0.05,
    padj_cutoff: Optional[float] = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    organism: str = "mouse",
    gene_case: Optional[str] = None,
    gene_set_source: str = "scatrans",
    include_gene_list: bool = False,
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

    pval_cutoff / padj_cutoff : float
        Cutoff applied to **adjusted p-values** (`p.adjust` column), not raw p-values.
        - `padj_cutoff` is the preferred modern name.
        - `pval_cutoff` is kept for backward compatibility and behaves identically.
        Default 0.05. If both are passed, `padj_cutoff` takes precedence.

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

    background : optional
        Deprecated alias of `universe`. Use `universe` instead.
        If both `universe` and `background` are provided, a ValueError is raised.
    """
    # Normalize organism early for consistency (attrs, resolver, etc.)
    organism_norm = str(organism).lower()

    # Reproducibility info for manuscript supplementary materials
    analysis_info = _get_analysis_info()

    # Resolve cutoff: prefer explicit padj_cutoff, fall back to pval_cutoff (legacy name)
    # Both are applied to the adjusted p-value column ("p.adjust").
    if padj_cutoff is not None:
        cutoff = float(padj_cutoff)
        if pval_cutoff != 0.05:
            _warn_user(
                "Both pval_cutoff and padj_cutoff were provided. "
                "Using padj_cutoff; pval_cutoff is ignored."
            )
    else:
        cutoff = float(pval_cutoff)

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

    if min_size < 1 or max_size < min_size or not (0 <= cutoff <= 1):
        raise ValueError("Invalid min_size, max_size or pval_cutoff/padj_cutoff")

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
    # If both are given, raise immediately to avoid silent "only universe wins" behavior.
    if universe is not None and background is not None:
        raise ValueError("Please provide only one of `universe` or `background`, not both. Use `universe`.")

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
    for term, term_genes in term_to_genes.items():
        term_genes_in_universe = term_genes & universe_set
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
            _log_info(f"Tested {len(term_to_genes)} terms; found 0 terms with overlap (after size filters)")
        return _empty_ora_result(
            include_gene_list=include_gene_list,
            reason="no_term_overlap_after_filters",
            n_tested_terms=len(term_to_genes),
            pval_cutoff=cutoff,
            universe_info=base_attrs.get("universe_info"),  # may be None here but ok
            **base_attrs,
        )

    res_df = pd.DataFrame(results)
    res_df["p.adjust"] = _bh_p_adjust(res_df["pvalue"].values)
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

    if return_all:
        if verbose:
            _log_info(f"Tested {len(res_df)} terms; returning all (return_all=True)")
        return res_df

    sig = res_df[res_df["p.adjust"] < cutoff].copy().reset_index(drop=True)
    sig.attrs.update(res_df.attrs)
    if verbose:
        _log_info(f"Tested {len(res_df)} terms; found {len(sig)} significant terms at padj < {cutoff}")
    return sig


def run_kegg(
    gene_list: Iterable[Any],
    organism: str = "mouse",
    universe: Optional[Iterable[Any]] = None,
    background: Optional[Iterable[Any]] = None,
    adata: Optional[Any] = None,  # forwarded to run_enrichment for smart universe default
    pval_cutoff: float = 0.05,
    padj_cutoff: Optional[float] = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: Optional[str] = None,
    kegg_library: Optional[str] = None,
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
    org_map = {"mouse": "Mouse", "mmu": "Mouse", "human": "Human", "hsa": "Human"}
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
    universe: Optional[Iterable[Any]] = None,
    background: Optional[Iterable[Any]] = None,
    adata: Optional[Any] = None,
    pval_cutoff: float = 0.05,
    padj_cutoff: Optional[float] = None,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    force_universe: bool = False,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: Optional[str] = None,
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
          (each sub-run does its own BH correction). Documented behavior for most users.
        - True: After concatenating BP+CC+MF, recompute a single BH correction across
          *all* tested GO terms using their raw pvalues. This is stricter/more conservative
          for cross-ontology multiple testing. Recommended when you want a single
          family-wise error control across the entire GO.

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
      - If adjust_across_all=True, the main "p.adjust" is the cross-GO BH; the column
        "p.adjust.within_ontology" preserves the original per-ontology adjusted values.

    Note on p.adjust for "ALL":
        For ontology="ALL", p.adjust is computed within each ontology separately
        before concatenation (unless adjust_across_all=True). Use per_ontology_attrs
        or the within_ontology column for full transparency.

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
        if ont not in ont_map and ont not in {"GO_BIOLOGICAL_PROCESS", "GO_CELLULAR_COMPONENT", "GO_MOLECULAR_FUNCTION"}:
            # allow passing the full Enrichr-style name directly
            gs_name = ont
        else:
            gs_name = ont_map.get(ont, ont)

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
    per_ontology_attrs: Dict[str, dict] = {}
    eff_cut = padj_cutoff if padj_cutoff is not None else pval_cutoff

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
    combined.attrs["per_ontology_attrs"] = per_ontology_attrs
    combined.attrs["analysis_info"] = analysis_info
    # Also carry forward gene_set_info etc from first if useful, but per_ontology is authoritative
    if frames:
        first_attrs = frames[0].attrs if hasattr(frames[0], "attrs") else {}
        for k in ("gene_set_info", "universe_info", "clusterprofiler_aligned"):
            if k in first_attrs and k not in combined.attrs:
                combined.attrs[k] = first_attrs[k]

    if adjust_across_all and not combined.empty and "pvalue" in combined.columns:
        # Preserve original within-ontology corrected values for transparency
        combined = combined.copy()
        combined["p.adjust.within_ontology"] = combined["p.adjust"]
        combined["p.adjust"] = _bh_p_adjust(combined["pvalue"].values)
        combined["neg_log10_padj"] = -np.log10(
            combined["p.adjust"].astype(float).clip(lower=1e-300)
        )
        combined = combined.sort_values("p.adjust").reset_index(drop=True)
        if verbose:
            _log_info("GO ALL: re-adjusted p.adjust across BP+CC+MF (adjust_across_all=True); "
                      "original per-ontology p.adjust saved in 'p.adjust.within_ontology'")

    if return_all:
        if verbose:
            _log_info(f"GO ALL: combined {len(combined)} terms (BP+CC+MF)")
        return combined

    # Apply cutoff to the combined table (using the effective cutoff)
    sig = combined[combined["p.adjust"] < eff_cut].copy().reset_index(drop=True)
    sig.attrs.update(combined.attrs)
    if verbose:
        _log_info(f"GO ALL: tested/combined {len(combined)} terms; found {len(sig)} significant at padj < {eff_cut}")
    return sig


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


# --- New export helpers for manuscript-ready supplementary tables (per review) ---


def _flatten_metadata(d: Any, prefix: str = "") -> list:
    """Recursively flatten nested dicts/lists into rows for human-readable Excel metadata sheet.

    Example:
        {"gene_set_info": {"actual_source": "bundled"}, "universe_info": {"effective_universe_size": 123}}
    becomes rows like:
        {"key": "gene_set_info.actual_source", "value": "bundled"}
        {"key": "universe_info.effective_universe_size", "value": "123"}
    """
    rows = []
    if isinstance(d, dict):
        for k, v in d.items():
            key = f"{prefix}.{k}" if prefix else k
            if isinstance(v, dict):
                rows.extend(_flatten_metadata(v, key))
            elif isinstance(v, (list, tuple)):
                try:
                    rows.append({"key": key, "value": json.dumps(v, ensure_ascii=False, default=str)})
                except Exception:
                    rows.append({"key": key, "value": str(v)})
            else:
                rows.append({"key": key, "value": str(v)})
    else:
        # top-level non-dict
        rows.append({"key": prefix or "value", "value": str(d)})
    return rows


def expand_enrichment_genes(res: pd.DataFrame) -> pd.DataFrame:
    """
    Expand semicolon-joined Genes column into a long term-gene table.

    One row per Term-Gene pair. Extremely convenient for supplementary
    materials or downstream gene-level analysis.

    If the input comes from run_go(ontology="ALL"), the "Ontology" column
    (BP/CC/MF) is preserved and placed first for clarity.

    Returns
    -------
    DataFrame with columns including [Ontology], Term, Description, Gene, plus
    the main stats columns from the original result.
    """
    has_ontology = "Ontology" in (res.columns if res is not None else [])

    if res is None or res.empty or "Genes" not in res.columns:
        base_cols = [
            "Term",
            "Description",
            "Gene",
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
            "TermSize",
        ]
        if has_ontology:
            base_cols = ["Ontology"] + base_cols
        return pd.DataFrame(columns=base_cols)

    rows = []
    for _, row in res.iterrows():
        genes_str = str(row.get("Genes", "") or "")
        genes = [g.strip() for g in genes_str.split(";") if g and g.strip() and g.strip().lower() != "nan"]

        for gene in genes:
            rec = {}
            if has_ontology:
                rec["Ontology"] = row.get("Ontology", "")
            rec.update({
                "Term": row.get("Term", ""),
                "Description": row.get("Description", ""),
                "Gene": gene,
                "Count": row.get("Count", None),
                "GeneRatio": row.get("GeneRatio", None),
                "GeneRatio_str": row.get("GeneRatio_str", ""),
                "BgRatio": row.get("BgRatio", None),
                "BgRatio_str": row.get("BgRatio_str", ""),
                "FoldEnrichment": row.get("FoldEnrichment", None),
                "RichFactor": row.get("RichFactor", None),
                "Overlap": row.get("Overlap", ""),
                "pvalue": row.get("pvalue", None),
                "p.adjust": row.get("p.adjust", None),
                "TermSize": row.get("TermSize", None),
            })
            rows.append(rec)

    df = pd.DataFrame(rows)
    # Ensure Ontology first if present
    if "Ontology" in df.columns:
        cols = ["Ontology"] + [c for c in df.columns if c != "Ontology"]
        df = df[cols]
    return df


def save_enrichment_report(
    res: pd.DataFrame,
    prefix: str = "enrichment",
    save_excel: bool = True,
    save_csv: bool = True,
    save_tsv: bool = False,
    save_metadata: bool = True,
    save_term_gene_table: bool = True,
    index: bool = False,
) -> Dict[str, str]:
    """
    Save enrichment results in formats friendly for manuscripts and supplementary materials.

    Produces a combination of:
      - {prefix}_results.csv / .tsv / .xlsx   (main table; Genes column is semicolon-joined)
      - {prefix}_term_gene_table.csv / .tsv / (in xlsx)   (long format: one row per term-gene pair)
      - {prefix}_metadata.json   (and a "metadata" sheet in xlsx)  (res.attrs + analysis provenance)

    List-typed columns (e.g. Genes_list when include_gene_list=True) are automatically
    converted to semicolon-joined strings for clean CSV/Excel/TSV export.

    The parent directory of `prefix` is created if it does not exist.

    Returns a dict with the written file paths, e.g.:
        {
            "results_csv": "..._results.csv",
            "term_gene_table_csv": "..._term_gene_table.csv",
            "metadata_json": "..._metadata.json",
            "results_xlsx": "..._results.xlsx",
            # plus _tsv variants if save_tsv=True
        }

    Examples
    --------
    res = run_kegg(genes, organism="mouse", return_all=True, include_gene_list=True)
    saved = save_enrichment_report(res, prefix="cluster1_kegg")
    # saved keys typically: 'results_csv', 'term_gene_table_csv', 'metadata_json', 'results_xlsx'

    # With TSV (great for Excel locale issues) and auto-created subdir:
    saved = save_enrichment_report(res, prefix="results/suppl/cluster1_kegg", save_tsv=True)

    # Long table for network analysis or gene-level follow-up:
    long_table = expand_enrichment_genes(res)
    # If from run_go(ontology="ALL"), long_table will have 'Ontology' as first column.
    """
    import json
    from pathlib import Path

    if res is None:
        raise ValueError("res cannot be None")

    prefix_path = Path(str(prefix))
    prefix_path.parent.mkdir(parents=True, exist_ok=True)
    outputs: Dict[str, str] = {}

    # Prepare a copy safe for export (convert list columns like Genes_list to joined strings)
    res_export = res.copy()
    for col in list(res_export.columns):
        try:
            if res_export[col].apply(lambda x: isinstance(x, (list, tuple))).any():
                res_export[col] = res_export[col].apply(
                    lambda x: ";".join(map(str, x)) if isinstance(x, (list, tuple)) else x
                )
        except Exception:
            pass  # be defensive for weird columns

    term_gene_df = None
    if save_term_gene_table:
        term_gene_df = expand_enrichment_genes(res)  # use original (with ; Genes) for expansion

    if save_csv:
        result_csv = f"{prefix_path}_results.csv"
        res_export.to_csv(result_csv, index=index)
        outputs["results_csv"] = result_csv

        if save_term_gene_table and term_gene_df is not None:
            term_gene_csv = f"{prefix_path}_term_gene_table.csv"
            term_gene_df.to_csv(term_gene_csv, index=False)
            outputs["term_gene_table_csv"] = term_gene_csv

    if save_tsv:
        result_tsv = f"{prefix_path}_results.tsv"
        res_export.to_csv(result_tsv, index=index, sep="\t")
        outputs["results_tsv"] = result_tsv

        if save_term_gene_table and term_gene_df is not None:
            term_gene_tsv = f"{prefix_path}_term_gene_table.tsv"
            term_gene_df.to_csv(term_gene_tsv, index=False, sep="\t")
            outputs["term_gene_table_tsv"] = term_gene_tsv

    if save_metadata:
        metadata_json = f"{prefix_path}_metadata.json"
        try:
            with open(metadata_json, "w", encoding="utf-8") as f:
                json.dump(res.attrs if hasattr(res, "attrs") else {}, f, indent=2, ensure_ascii=False, default=str)
            outputs["metadata_json"] = metadata_json
        except Exception as e:
            _warn_user(f"Could not write metadata JSON: {e}")

    if save_excel:
        excel_file = f"{prefix_path}_results.xlsx"
        try:
            with pd.ExcelWriter(excel_file) as writer:
                res_export.to_excel(writer, sheet_name="enrichment_results", index=index)

                if save_term_gene_table and term_gene_df is not None:
                    term_gene_df.to_excel(writer, sheet_name="term_gene_table", index=False)

                if save_metadata:
                    attrs = res.attrs if hasattr(res, "attrs") else {}
                    # Flattened for easy human reading in Excel (nested keys)
                    meta_rows = _flatten_metadata(attrs)
                    meta_df = pd.DataFrame(meta_rows)
                    meta_df.to_excel(writer, sheet_name="metadata", index=False)

            outputs["results_xlsx"] = excel_file
        except Exception as e:
            _warn_user(f"Could not write Excel report: {e}")

    return outputs
