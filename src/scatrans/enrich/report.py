"""scatrans.enrich.report — internal package module."""

from __future__ import annotations

import json
import logging
import re
from typing import Any

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

    # GSEA results carry the gene list in ``leading_edge`` (renamed from gseapy
    # ``Lead_genes``) rather than the ORA ``Genes`` column. Alias it so GSEA inputs
    # expand instead of silently returning an empty table.
    if res is not None and not res.empty and "Genes" not in res.columns:
        for _alt in ("leading_edge", "Lead_genes"):
            if _alt in res.columns:
                res = res.copy()
                res["Genes"] = res[_alt]
                break
        else:
            logger.warning(
                "expand_enrichment_genes: enrichment frame has %d rows but no usable "
                "gene-list column (Genes / leading_edge / Lead_genes); the term-gene "
                "table will be empty. Available columns: %s",
                len(res),
                list(res.columns),
            )

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
        # ORA uses ';'; GSEA leading_edge / some exports use ',' — accept both.
        genes = [
            g.strip()
            for g in re.split(r"[;,]+", genes_str)
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

    Returns a dict with the written file paths, e.g.::

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
