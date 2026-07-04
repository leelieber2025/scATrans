from __future__ import annotations

import functools
import json
import logging
import math
import os
import re
import sys
import warnings
from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Reproducibility metadata for bundled 2026 gene-set libraries (see data/README.md).
BUNDLED_GENE_SET_PROVENANCE: dict[str, dict[str, Any]] = {
    "Hs_GO_Biological_Process_2026": {
        "bundled_file": "Hs_GO_Biological_Process_2026.txt",
        "species": "Homo sapiens",
        "collection": "GO Biological Process",
        "source_pipeline": "clusterProfiler::buildGOmap + org.Hs.eg.db + GO.db (Bioconductor)",
        "data_license": "GO / Bioconductor-derived (not Apache-2.0); see scatrans/data/DATA_LICENSES.md",
        "license_url": "http://geneontology.org/docs/go-licenses/",
        "package_version_tag": "scATrans 0.9.x data bundle",
        "extracted_date": "2026-06",
        "format": "GMT-like tab-separated (.txt)",
        "gene_id_type": "SYMBOL",
        "n_terms": 14208,
    },
    "Mm_GO_Biological_Process_2026": {
        "bundled_file": "Mm_GO_Biological_Process_2026.txt",
        "species": "Mus musculus",
        "collection": "GO Biological Process",
        "source_pipeline": "clusterProfiler::buildGOmap + org.Mm.eg.db + GO.db (Bioconductor)",
        "data_license": "GO / Bioconductor-derived (not Apache-2.0); see scatrans/data/DATA_LICENSES.md",
        "license_url": "http://geneontology.org/docs/go-licenses/",
        "package_version_tag": "scATrans 0.9.x data bundle",
        "extracted_date": "2026-06",
        "format": "GMT-like tab-separated (.txt)",
        "gene_id_type": "SYMBOL",
        "n_terms": 14956,
    },
    "Hs_KEGG_2026": {
        "bundled_file": "Hs_KEGG_2026.txt",
        "species": "Homo sapiens",
        "collection": "KEGG pathways",
        "source_pipeline": "clusterProfiler KEGG cache (organism hsa) via enrichKEGG mappings",
        "data_license": "KEGG (not Apache-2.0); academic use with attribution; commercial requires KEGG license",
        "license_url": "https://www.kegg.jp/kegg/legal.html",
        "package_version_tag": "scATrans 0.9.x data bundle",
        "extracted_date": "2026-06",
        "format": "GMT-like tab-separated (.txt)",
        "gene_id_type": "SYMBOL",
        "n_terms": 222,
    },
    "Mm_KEGG_2026": {
        "bundled_file": "Mm_KEGG_2026.txt",
        "species": "Mus musculus",
        "collection": "KEGG pathways",
        "source_pipeline": "clusterProfiler KEGG cache (organism mmu) via enrichKEGG mappings",
        "data_license": "KEGG (not Apache-2.0); academic use with attribution; commercial requires KEGG license",
        "license_url": "https://www.kegg.jp/kegg/legal.html",
        "package_version_tag": "scATrans 0.9.x data bundle",
        "extracted_date": "2026-06",
        "format": "GMT-like tab-separated (.txt)",
        "gene_id_type": "SYMBOL",
        "n_terms": 218,
    },
}

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
    content: str, gene_case: str | None = None
) -> tuple[dict[str, set], dict[str, str]]:
    """Parse GMT text content (term<TAB>desc<TAB>gene1<TAB>gene2...).

    The second column (description) is retained when present and stored
    under the "Description" output column. Many bundled sets ship with an
    empty description column (two tabs); this is handled gracefully.
    """
    term_to_genes: dict[str, set] = {}
    term_to_desc: dict[str, str] = {}
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


def list_bundled_gene_sets(verbose: bool = False) -> list[str]:
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
    discovered: list[str] = []
    try:
        data_traversable = files("scatrans.data")
        for item in data_traversable.iterdir():
            name = getattr(item, "name", str(item))
            if name.endswith((".gmt", ".tsv", ".txt")) or "gene_set" in name.lower():
                discovered.append(name)
    except Exception as exc:
        logger.debug("list_bundled_gene_sets discovery failed (best-effort): %s", exc)
        discovered = []

    # Fallback known names. These must match files actually present under src/scatrans/data/
    # (packaged via pyproject.toml package-data). Keep in sync with shipped files.
    known = [
        "Hs_GO_Biological_Process_2026.txt",
        "Hs_KEGG_2026.txt",
        "Mm_GO_Biological_Process_2026.txt",
        "Mm_KEGG_2026.txt",
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


def _bundled_provenance_for(resolved_name: str) -> dict[str, Any]:
    """Return provenance dict for a bundled library basename (without extension)."""
    base = resolved_name
    for ext in (".txt", ".gmt", ".tsv"):
        if base.lower().endswith(ext):
            base = base[: -len(ext)]
    return dict(BUNDLED_GENE_SET_PROVENANCE.get(base, {}))


@functools.lru_cache(maxsize=16)
def _try_load_bundled_gene_set(
    name: str, gene_case: str | None = None
) -> tuple[dict[str, set], dict[str, str], str] | None:
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
            terms, desc = _parse_gmt_content(content, gene_case=gene_case)
            loaded_as = cand
            if not loaded_as.endswith((".txt", ".gmt", ".tsv")):
                loaded_as = cand
            base = loaded_as
            for ext in (".txt", ".gmt", ".tsv"):
                if base.lower().endswith(ext):
                    base = base[: -len(ext)]
            return terms, desc, base
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
        # Base names (no year) -> organism 2026 built-in (only BP + KEGG are bundled)
        "GO_Biological_Process": f"{prefix}_GO_Biological_Process_2026",
        "GO_Biological_Process_2026": f"{prefix}_GO_Biological_Process_2026",
        "GO_BP": f"{prefix}_GO_Biological_Process_2026",
        "KEGG": f"{prefix}_KEGG_2026",
        "KEGG_2026": f"{prefix}_KEGG_2026",
        # Historical Enrichr names with explicit years (e.g. GO_Biological_Process_2023,
        # KEGG_2021) are intentionally NOT mapped here — they pass through to gseapy/Enrichr.
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

GSEA_COLUMNS = [
    "Term",
    "Description",
    "ES",
    "NES",
    "pvalue",
    "p.adjust",
    "neg_log10_padj",
    "leading_edge",
    "Tag_percent",
    "Gene_percent",
    "TermSize",
]


def _get_analysis_info() -> dict[str, Any]:
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


def _resolve_enrichment_padj_cutoff(
    pval_cutoff: float | None,
    padj_cutoff: float | None,
) -> float:
    """Resolve the effective adjusted-p cutoff; warn on legacy ``pval_cutoff`` usage."""
    if padj_cutoff is not None:
        cutoff = float(padj_cutoff)
        if pval_cutoff is not None:
            _warn_user(
                "Both pval_cutoff and padj_cutoff were provided. "
                "Using padj_cutoff; pval_cutoff is ignored."
            )
        return cutoff
    if pval_cutoff is not None:
        _warn_user(
            "`pval_cutoff` is deprecated (it applies to *adjusted* p-values, not raw p-values). "
            "Please use `padj_cutoff` instead for clarity."
        )
        return float(pval_cutoff)
    return 0.05


def _apply_gene_case(genes: Iterable[Any], gene_case: str | None) -> list[str]:
    if gene_case is not None:
        gene_case = str(gene_case).lower()
    if gene_case is None:
        return [str(g).strip() for g in genes]
    if gene_case == "upper":
        return [str(g).strip().upper() for g in genes]
    if gene_case == "lower":
        return [str(g).strip().lower() for g in genes]
    raise ValueError("gene_case must be None, 'upper', or 'lower'")


def _resolve_gseapy_weight(
    *,
    weight: float | None = None,
    weighted_score_type: str | float | None = None,
) -> float:
    """Map Broad/GSEA score-type names to gseapy's numeric ``weight`` (p exponent) parameter."""
    if weight is not None:
        return float(weight)
    if weighted_score_type is None:
        return 1.0
    if isinstance(weighted_score_type, (int, float)):
        return float(weighted_score_type)
    wst = str(weighted_score_type).strip().lower()
    if wst in {"classic", "0", "0.0", "unweighted", "none"}:
        return 0.0
    if wst in {"weighted", "1", "1.0"}:
        return 1.0
    try:
        return float(wst)
    except ValueError as e:
        raise ValueError(
            "weighted_score_type must be numeric or one of "
            f"'classic'/'weighted'/'unweighted', got {weighted_score_type!r}"
        ) from e


def _clean_gene_list(gene_list: Iterable[Any] | None, gene_case: str | None = None) -> list[str]:
    if gene_list is None:
        return []
    s = pd.Series(list(gene_list))
    s = s.dropna()
    s = s.astype(str).str.strip()
    s = s[(s != "") & (s.str.lower() != "nan")]
    cleaned = _apply_gene_case(s.tolist(), gene_case)
    return pd.Series(cleaned).drop_duplicates().tolist()


def _load_gene_sets(
    gene_sets_input: Mapping[str, Iterable[Any]] | str | None,
    organism: str = "mouse",
    verbose: bool = True,
    gene_case: str | None = None,
) -> tuple[dict[str, set], dict[str, str], dict[str, Any]]:
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
        term_to_genes: dict[str, set] = {}
        term_to_desc: dict[str, str] = {}
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
        p = Path(gene_sets_input).expanduser()
        looks_like_path = (
            p.exists()
            or os.path.isabs(gene_sets_input)
            or gene_sets_input.startswith("~")
            or "/" in gene_sets_input
            or "\\" in gene_sets_input
        )
        if looks_like_path:
            if not p.exists():
                raise FileNotFoundError(f"GMT file not found: {gene_sets_input}")
            with open(p, encoding="utf-8") as f:
                content = f.read()
            term_to_genes, term_to_desc = _parse_gmt_content(content, gene_case=gene_case)
            load_info = {
                "actual_source": "gmt",
                "library_name": os.path.basename(gene_sets_input),
                "resolved_name": gene_sets_input,
                "path": str(p.resolve()) if p.exists() else os.path.abspath(gene_sets_input),
            }
            return term_to_genes, term_to_desc, load_info

        # 1. Try bundled package data first (ClusterProfiler-derived GO/KEGG etc.)
        bundled = _try_load_bundled_gene_set(gene_sets_input, gene_case=gene_case)
        if bundled is not None:
            term_to_genes, term_to_desc, bundled_key = bundled
            if verbose:
                _log_info(f"Loaded bundled gene set '{gene_sets_input}' from package data")
            prov = _bundled_provenance_for(bundled_key)
            load_info = {
                "actual_source": "bundled",
                "library_name": gene_sets_input,
                "resolved_name": bundled_key,
                "path": None,
                "provenance": prov,
            }
            return term_to_genes, term_to_desc, load_info

        # 2. Fall back to gseapy / Enrichr (original behavior)
        try:
            import gseapy as gp

            try:
                gene_dict = gp.get_library(name=gene_sets_input, organism=organism)
            except Exception as exc:
                logger.debug("gseapy get_library with organism failed, retrying without: %s", exc)
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
                b_terms, b_desc, bundled_key = bundled2
                prov = _bundled_provenance_for(bundled_key)
                load_info = {
                    "actual_source": "bundled",
                    "library_name": gene_sets_input,
                    "resolved_name": bundled_key,
                    "path": None,
                    "provenance": prov,
                }
                return b_terms, b_desc, load_info
            raise ValueError(
                f"Failed to load '{gene_sets_input}' via gseapy or as bundled package data.\n"
                f"Error ({type(e).__name__}): {e}\n"
                f"Available bundled sets: {list_bundled_gene_sets()}"
            ) from e

    raise ValueError(
        "gene_sets must be dict, GMT path or gseapy library name (or a bundled set name)"
    )


def _bh_p_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR control.

    NaN-tolerant: non-finite p-values are left as NaN in the output (consistent
    with statsmodels.stats.multitest.multipletests behavior used elsewhere).
    """
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return np.array([])
    finite_mask = np.isfinite(pvalues)
    out = np.full(n, np.nan, dtype=float)
    if not finite_mask.any():
        return out
    finite_vals = pvalues[finite_mask]
    m = len(finite_vals)
    sorted_idx_local = np.argsort(finite_vals)
    sorted_p = finite_vals[sorted_idx_local]
    adjusted = sorted_p * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    # map adjusted values back to original positions among the finite ones
    adjusted_full = np.empty(m)
    adjusted_full[sorted_idx_local] = adjusted
    out[finite_mask] = adjusted_full
    return out


def _apply_p_adjust(pvalues: np.ndarray, method: str = "fdr_bh") -> np.ndarray:
    """Apply multiple-testing correction across all tested terms in one ORA call.

    For standard fdr_bh / bonferroni we delegate to statsmodels.multipletests
    when available (consistent with _de.py and _permutation.py). Falls back to
    our implementations otherwise.
    """
    method_norm = str(method).lower().replace("-", "_")
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return np.array([])
    if method_norm in ("none", "raw", "off"):
        return pvalues.copy()

    if method_norm in ("fdr_bh", "bh", "fdr", "padj"):
        try:
            from statsmodels.stats.multitest import multipletests

            # multipletests returns (reject, pvals_corrected, alphacSidak, alphacBonf)
            _, padj, _, _ = multipletests(pvalues, alpha=0.05, method="fdr_bh", is_sorted=False)
            return np.asarray(padj, dtype=float)
        except Exception:
            return _bh_p_adjust(pvalues)

    if method_norm in ("bonferroni", "bonf"):
        try:
            from statsmodels.stats.multitest import multipletests

            _, padj, _, _ = multipletests(pvalues, alpha=0.05, method="bonferroni")
            return np.asarray(padj, dtype=float)
        except Exception:
            return np.minimum(pvalues * n, 1.0)
    raise ValueError(
        f"p_adjust_method={method!r} not recognized. Use one of: 'fdr_bh', 'bonferroni', 'none'."
    )


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


def _empty_gsea_result(**attrs) -> pd.DataFrame:
    """Return an empty GSEA result DataFrame with GSEA_COLUMNS and .attrs populated."""
    df = pd.DataFrame(columns=list(GSEA_COLUMNS))
    df.attrs.update(attrs)
    return df


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
            if "scatrans" in adata.uns and "raw_gene_list" in adata.uns["scatrans"]:
                preserved = adata.uns["scatrans"]["raw_gene_list"]
                if preserved:
                    provided = preserved
                    if verbose:
                        _log_info(
                            "Using preserved raw_gene_list from adata.uns['scatrans'] "
                            f"({len(preserved)} genes) as universe (from previous store_raw_counts)."
                        )
        except Exception as _e:
            if verbose:
                _log_info(f"Could not read raw_gene_list from adata: {_e}")
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
    n_terms_size_filtered = 0
    for term, term_genes in term_to_genes.items():
        term_genes_in_universe = term_genes & universe_set
        K = len(term_genes_in_universe)
        if min_size > K or max_size < K:
            continue
        n_terms_size_filtered += 1
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
        "n_terms_size_filtered": int(n_terms_size_filtered),
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
            _log_info(
                "GO ALL: re-adjusted p.adjust across BP+CC+MF (adjust_across_all=True); "
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


def _comb_fraction(n: int, k: int) -> Fraction:
    """Binomial coefficient C(n, k) as an exact Fraction (PathwayDenester)."""
    if k > n or k < 0:
        return Fraction(0)
    return Fraction(math.factorial(n), math.factorial(k) * math.factorial(n - k))


def _comb_comb_comb(
    degs_in_test: int,
    degs_in_intersection: int,
    intersection_size: int,
    size_test: int,
) -> float:
    """PathwayDenester independence test (combinatorial hypergeometric sum)."""
    if intersection_size <= 0 or size_test <= 0 or degs_in_test <= 0:
        return 1.0
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
                kept_sets = []
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
                    rows.append(
                        {"key": key, "value": json.dumps(v, ensure_ascii=False, default=str)}
                    )
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

    rows = []
    for _, row in res.iterrows():
        genes_str = str(row.get("Genes", "") or "")
        genes = [
            g.strip()
            for g in genes_str.split(";")
            if g and g.strip() and g.strip().lower() != "nan"
        ]

        for gene in genes:
            rec = {}
            if has_ontology:
                rec["Ontology"] = row.get("Ontology", "")
            rec.update(
                {
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
                }
            )
            rows.append(rec)

    df = pd.DataFrame(rows, columns=base_cols) if not rows else pd.DataFrame(rows)
    # Ensure Ontology first if present
    if "Ontology" in df.columns:
        cols = ["Ontology"] + [c for c in df.columns if c != "Ontology"]
        df = df[cols]
    return df


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
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Pre-ranked Gene Set Enrichment Analysis (GSEA) using gseapy.prerank.

    This implements the classic GSEA algorithm on a user-provided ranked gene list
    (e.g. logFC, t-statistic or custom score from active_score / DE results).
    It is the Python equivalent of clusterProfiler::GSEA / Broad GSEA Preranked.

    Parameters
    ----------
    ranked_genes : pd.Series, dict, or list-like
        - Preferred: pd.Series with gene names as index and numeric scores as values.
          Higher score = more "up" in target group (e.g. logFC).
          The function will sort internally if needed.
        - dict: gene -> score
        - list of genes: treated as pre-sorted from high to low (scores assigned decreasing).
        Gene names will be cleaned according to gene_case.
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
        neg_log10_padj, leading_edge, etc. Sorted by |NES| (absolute value) descending
        so that the strongest magnitude effects (positive or negative) appear first.
        Rich metadata in .attrs (method="gsea_prerank", gene_set_info, nperm, gsea_info, analysis_info).

    Notes
    -----
    - Unlike ORA, GSEA does not use an explicit "universe" in the same way; the ranked
      list itself defines the background. min_size/max_size still apply.
    - Requires gseapy. Install via `pip install gseapy` or `pip install "scatrans[gsea]"`.
    - For best results with scATrans, pass a Series derived from active_score results, e.g.:
        ranked = all_results.set_index("gene")["logFC"]   # or suitable score
        res = scat.run_gsea(ranked, gene_sets="GO_Biological_Process")
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

    # Normalize ranked_genes input to pd.Series (gene -> score)
    if isinstance(ranked_genes, pd.DataFrame):
        if ranked_genes.shape[1] >= 2:
            # Column 0 = gene names, Column 1 = scores (original intent)
            gene_names = ranked_genes.iloc[:, 0].astype(str).values
            scores = ranked_genes.iloc[:, 1].values
            ranked_genes = pd.Series(scores, index=gene_names)
        else:
            ranked_genes = ranked_genes.iloc[:, 0]
    if isinstance(ranked_genes, (list, tuple)):
        # treat as pre-ordered high->low, assign descending ranks
        genes = _apply_gene_case([str(g).strip() for g in ranked_genes], gene_case)
        scores = list(range(len(genes), 0, -1))
        ranked = pd.Series(scores, index=genes)
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
    ranked = ranked[~ranked.index.duplicated(keep="first")]  # keep first occurrence
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
        msg = f"gseapy.prerank failed: {e}"
        if verbose:
            _log_info(msg)
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return _empty_gsea_result(
            method="gsea_prerank",
            organism=organism_norm,
            gene_case=gene_case,
            reason="gseapy_error",
            error=str(e),
            gene_set_info=gene_set_info,
            analysis_info=analysis_info,
        )

    if res_df is None or res_df.empty:
        overlap = len(set(ranked.index) & set().union(*term_to_genes.values()))
        msg = (
            "gseapy returned no results (all gene sets filtered out?). "
            f"Ranked genes={len(ranked)}, overlap with gene sets={overlap}. "
            "Check gene symbols/IDs match the library; try lowering min_size."
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
        "score_min": float(ranked.min()),
        "score_max": float(ranked.max()),
        "score_median": float(ranked.median()),
    }
    # Store gseapy internals for accurate gseaplot (RES curve + hits per term)
    if hasattr(pre_res, "results"):
        res_df.attrs["gsea_details"] = pre_res.results
    if hasattr(pre_res, "ranking"):
        res_df.attrs["ranking"] = pre_res.ranking.to_dict()

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


def save_enrichment_report(
    res: pd.DataFrame,
    prefix: str = "enrichment",
    save_excel: bool = True,
    save_csv: bool = True,
    save_tsv: bool = False,
    save_metadata: bool = True,
    save_term_gene_table: bool = True,
    index: bool = False,
) -> dict[str, str]:
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
    outputs: dict[str, str] = {}

    _MAX_EXPORT_CELL_CHARS = 30_000  # stay below Excel's ~32k cell limit

    def _join_list_for_export(x: Any) -> Any:
        if isinstance(x, (list, tuple)):
            joined = ";".join(map(str, x))
            if len(joined) > _MAX_EXPORT_CELL_CHARS:
                n = len(x)
                joined = joined[:_MAX_EXPORT_CELL_CHARS] + f"...(truncated, {n} genes total)"
            return joined
        return x

    # Prepare a copy safe for export (convert list columns like Genes_list to joined strings)
    res_export = res.copy()
    for col in list(res_export.columns):
        try:
            if res_export[col].apply(lambda x: isinstance(x, (list, tuple))).any():
                res_export[col] = res_export[col].apply(_join_list_for_export)
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
                json.dump(
                    res.attrs if hasattr(res, "attrs") else {},
                    f,
                    indent=2,
                    ensure_ascii=False,
                    default=str,
                )
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


# =============================================================================
# Multi-group / compareCluster-style support
# =============================================================================


def _clean_and_validate_gene_list_for_compare(
    genes: Iterable[Any], gene_case: str | None = None
) -> list[str]:
    """Lightweight cleaner used by compare helpers (re-uses the main cleaner)."""
    return _clean_gene_list(genes, gene_case=gene_case)


def _row_gene_ids_from_df(work: pd.DataFrame) -> list[str]:
    """Return per-row gene identifiers, preferring explicit DE columns over index."""
    if work is None or work.empty:
        return []
    if "gene" in work.columns:
        return work["gene"].astype(str).tolist()
    if "names" in work.columns:
        return work["names"].astype(str).tolist()
    if (
        isinstance(work.index, pd.Index)
        and len(work.index) > 0
        and not isinstance(work.index, pd.RangeIndex)
    ):
        return work.index.astype(str).tolist()
    if work.shape[1] > 0:
        return work.iloc[:, 0].astype(str).tolist()
    return []


def extract_gene_lists(
    de_results: pd.DataFrame | Mapping[str, pd.DataFrame],
    *,
    logfc_cutoff: float = 0.5,
    pval_cutoff: float = 0.05,
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
        Max p_adj (or p_val if no p_adj column).
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
        de_dict, logfc_cutoff=0.5, pval_cutoff=0.05, logfc_direction="up"
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

    def _get_genes_from_df(df: pd.DataFrame) -> list[str]:
        if df is None or df.empty:
            return []
        work = df
        genes = _row_gene_ids_from_df(work)

        if "p_adj" in work.columns:
            padj = pd.to_numeric(work["p_adj"], errors="coerce")
        elif "p_val" in work.columns:
            padj = pd.to_numeric(work["p_val"], errors="coerce")
        else:
            padj = pd.Series(1.0, index=work.index)

        if "logFC" in work.columns:
            lfc = pd.to_numeric(work["logFC"], errors="coerce")
        else:
            lfc = pd.Series(0.0, index=work.index)

        mask = padj < float(pval_cutoff)
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
            # (re-implement light logic to split)
            work = de_results
            if "p_adj" in work.columns:
                padj = pd.to_numeric(work["p_adj"], errors="coerce")
            elif "p_val" in work.columns:
                padj = pd.to_numeric(work["p_val"], errors="coerce")
            else:
                padj = pd.Series(1.0, index=work.index)
            lfc = (
                pd.to_numeric(work.get("logFC", 0), errors="coerce")
                if "logFC" in work.columns
                else pd.Series(0, index=work.index)
            )
            sig = padj < float(pval_cutoff)
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
            if "p_adj" in work.columns:
                padj = pd.to_numeric(work["p_adj"], errors="coerce")
            elif "p_val" in work.columns:
                padj = pd.to_numeric(work["p_val"], errors="coerce")
            else:
                padj = pd.Series(
                    1.0, index=work.index if hasattr(work, "index") else range(len(work))
                )
            lfc = (
                pd.to_numeric(work.get("logFC", 0), errors="coerce")
                if hasattr(work, "get") and "logFC" in getattr(work, "columns", [])
                else pd.Series(0, index=getattr(work, "index", range(len(work))))
            )
            sig = padj < float(pval_cutoff)
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
    ``adjust_across_clusters=False`` (the default) performs BH FDR correction
    *separately inside each cluster*. This is conservative per group but when you
    have many clusters the *overall* false discovery rate across the whole table
    can be higher than a single global correction. If you intend to compare
    significance across clusters, use ``adjust_across_clusters=True``; per-cluster
    calls then use ``return_all=True`` internally (all size-eligible terms, including
    zero-overlap) before a single global BH step, and the final table is filtered
    by ``padj_cutoff``/``pval_cutoff`` unless you also pass ``return_all=True``.

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
    eff_return_all = return_all
    if adjust_across_clusters:
        if not return_all:
            _warn_user(
                "compare_enrichment: adjust_across_clusters=True requires every tested term "
                "from each cluster (including zero-overlap and non-significant) before global "
                "BH correction. Forcing return_all=True on per-cluster calls; significance "
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
                "multiple_testing": {"scope": "per_cluster", "method": "BH", "n_tests": 0},
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

    # Multi-cluster multiple testing (new, analogous to adjust_across_all in run_go)
    multiple_testing_scope = "per_cluster"
    if adjust_across_clusters and not combined.empty and "pvalue" in combined.columns:
        combined = combined.copy()
        if "p.adjust" in combined.columns:
            combined["p.adjust.within_cluster"] = combined["p.adjust"]
        combined["p.adjust"] = _bh_p_adjust(combined["pvalue"].values)
        combined["neg_log10_padj"] = -np.log10(
            combined["p.adjust"].astype(float).clip(lower=1e-300)
        )
        combined = combined.sort_values("p.adjust").reset_index(drop=True)
        multiple_testing_scope = "all_clusters"
        if verbose:
            _log_info(
                "compare_enrichment: p.adjust re-computed across all clusters "
                f"(adjust_across_clusters=True; n_tests={len(combined)})"
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
                "method": "BH",
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
