"""scatrans.enrich._data — internal package module."""

from __future__ import annotations

import logging
import os
import sys
import warnings
from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

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
        "package_version_tag": "scATrans 0.10.x data bundle",
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
        "package_version_tag": "scATrans 0.10.x data bundle",
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
        "package_version_tag": "scATrans 0.10.x data bundle",
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
        "package_version_tag": "scATrans 0.10.x data bundle",
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


class _DeepcopyImmuneDict(dict):
    """A dict whose ``copy.deepcopy`` is an O(1) identity return.

    ``pandas.DataFrame.attrs`` is deep-copied on essentially every DataFrame
    operation (slicing, ``.head()``, ``.copy()``, column assignment) via
    ``NDFrame.__finalize__``. ``run_gsea`` stores gseapy's full per-term
    running-enrichment-score curves in ``.attrs["gsea_details"]`` — for a
    genome-wide ranked list against thousands of gene sets this is tens of
    millions of floats, so deep-copying it on every downstream DataFrame
    operation (as happens inside ``enrich_dotplot``, ``gseaplot``,
    ``filter_active_genes``, or even a plain ``.head()`` in a notebook) makes
    those calls take seconds to minutes instead of being instantaneous. This
    payload is write-once/read-only after ``run_gsea`` returns, so sharing
    the same object across copies instead of duplicating it is safe.
    """

    def __deepcopy__(self, memo):
        return self


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


def _check_gene_set_mapping_rate(
    input_genes: Iterable[Any],
    reference_genes: Iterable[Any],
    *,
    context: str = "enrichment",
    threshold: float = 0.2,
    n_examples: int = 5,
    gene_case: str | None = None,
) -> dict[str, Any]:
    """Compute input→gene-set (or universe) mapping rate; warn when low.

    Shared by ORA (:func:`run_enrichment`) and GSEA (:func:`run_gsea`) so both
    surfaces catch species / case / ID-type mismatches before silent empty results.

    Parameters
    ----------
    input_genes
        Ranked genes (GSEA) or query gene list (ORA), already case-normalized.
    reference_genes
        Gene-set members or effective universe symbols (same case as input).
    context
        Label in the warning (e.g. ``"run_gsea"``, ``"run_enrichment"``).
    threshold
        Warn when mapped/input is strictly below this fraction (default 20%).
    gene_case
        Current gene_case setting, included in the hint when mapping is low.

    Returns
    -------
    dict with n_input, n_mapped, mapping_rate, example_input, example_reference.
    """
    genes = [str(g).strip() for g in input_genes if str(g).strip()]
    ref_list = [str(g).strip() for g in reference_genes if str(g).strip()]
    ref_set = set(ref_list)
    mapped = [g for g in genes if g in ref_set]
    n_input = len(genes)
    n_mapped = len(mapped)
    rate = float(n_mapped) / float(max(n_input, 1))
    info: dict[str, Any] = {
        "n_input": int(n_input),
        "n_mapped": int(n_mapped),
        "mapping_rate": rate,
        "example_input": genes[:n_examples],
        "example_reference": ref_list[:n_examples],
        "threshold": float(threshold),
        "gene_case": gene_case,
    }
    if n_input > 0 and rate < threshold:
        case_hint = (
            f" Current gene_case={gene_case!r}."
            if gene_case is not None
            else " gene_case is None (literal match only)."
        )
        _warn_user(
            f"Low mapping rate ({rate:.1%}) in {context} "
            f"({n_mapped}/{n_input} input genes found in gene sets / universe). "
            f"Input examples: {info['example_input']}; "
            f"gene-set/universe examples: {info['example_reference']}.{case_hint} "
            "Check gene ID type, organism, and gene_case "
            "(Enrichr libraries are typically UPPERCASE symbols — try gene_case='upper' "
            "for mixed-case mouse/human symbols such as Tp53)."
        )
    return info


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


def _expand_gene_list_input(gene_list: Iterable[Any] | None) -> list[str]:
    """Expand list/Series/DataFrame inputs to raw gene identifier strings."""
    if gene_list is None:
        return []
    if isinstance(gene_list, pd.DataFrame):
        work = gene_list
        if work.empty:
            return []
        if "gene" in work.columns:
            return work["gene"].astype(str).tolist()
        if "names" in work.columns:
            return work["names"].astype(str).tolist()
        if "symbol" in work.columns:
            return work["symbol"].astype(str).tolist()
        if (
            isinstance(work.index, pd.Index)
            and len(work.index) > 0
            and not isinstance(work.index, pd.RangeIndex)
        ):
            return work.index.astype(str).tolist()
        # RangeIndex + no gene column: first column is often logFC/score, NOT genes.
        # Only treat col0 as IDs when it looks non-numeric (e.g. symbols).
        if work.shape[1] > 0:
            col0 = work.iloc[:, 0]
            as_str = col0.astype(str)
            numeric_frac = float(pd.to_numeric(as_str, errors="coerce").notna().mean())
            if numeric_frac < 0.5:
                return as_str.tolist()
            logger.warning(
                "DataFrame has a default RangeIndex and no 'gene'/'names' column; "
                "first column looks numeric (e.g. logFC), so gene IDs cannot be inferred. "
                "Use gene symbols as the index or add a 'gene' column. Returning empty list."
            )
            return []
        return []
    if isinstance(gene_list, pd.Series):
        ser = gene_list
        if (
            isinstance(ser.index, pd.Index)
            and len(ser.index) > 0
            and not isinstance(ser.index, pd.RangeIndex)
        ):
            as_str = ser.index.astype(str)
            numeric_frac = float(pd.to_numeric(as_str, errors="coerce").notna().mean())
            if numeric_frac < 0.5:
                return as_str.tolist()
        return ser.dropna().astype(str).tolist()
    if isinstance(gene_list, np.ndarray):
        return gene_list.astype(str).tolist()
    return [str(g) for g in gene_list]


def _clean_gene_list(gene_list: Iterable[Any] | None, gene_case: str | None = None) -> list[str]:
    if gene_list is None:
        return []
    s = pd.Series(_expand_gene_list_input(gene_list))
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
        load_info: dict[str, Any] = {
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

            # gseapy expects title-case species ("Mouse"/"Human"); our public API
            # uses lowercase ("mouse"/"human"). Map before calling Enrichr.
            _gseapy_org_map = {
                "mouse": "Mouse",
                "mm": "Mouse",
                "mmu": "Mouse",
                "human": "Human",
                "hs": "Human",
                "hsa": "Human",
            }
            gseapy_organism = _gseapy_org_map.get(str(organism).lower(), organism)
            try:
                gene_dict = gp.get_library(name=gene_sets_input, organism=gseapy_organism)
            except Exception as exc:
                logger.warning(
                    "gseapy get_library(name=%r, organism=%r) failed (%s); "
                    "retrying without organism (may default to Human libraries — "
                    "verify gene_set_info.actual_source and mapping rates).",
                    gene_sets_input,
                    gseapy_organism,
                    exc,
                )
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

    Non-finite p-values are left as NaN in the output and **excluded** from the
    multipletests input (statsmodels is not NaN-safe: one NaN can poison the
    whole adjusted vector).
    """
    method_norm = str(method).lower().replace("-", "_")
    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    if n == 0:
        return np.array([])
    if method_norm in ("none", "raw", "off"):
        return pvalues.copy()

    finite = np.isfinite(pvalues)
    out = np.full(n, np.nan, dtype=float)
    if not finite.any():
        return out

    def _adjust_finite(vals: np.ndarray, mt_method: str) -> np.ndarray:
        try:
            from statsmodels.stats.multitest import multipletests

            _, padj, _, _ = multipletests(vals, alpha=0.05, method=mt_method, is_sorted=False)
            return np.asarray(padj, dtype=float)
        except Exception:
            if mt_method == "fdr_bh":
                return _bh_p_adjust(vals)
            return np.minimum(vals * len(vals), 1.0)

    if method_norm in ("fdr_bh", "bh", "fdr", "padj"):
        out[finite] = _adjust_finite(pvalues[finite], "fdr_bh")
        return out

    if method_norm in ("bonferroni", "bonf"):
        out[finite] = _adjust_finite(pvalues[finite], "bonferroni")
        return out
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
