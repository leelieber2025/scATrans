import logging
import os
import re
import warnings
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

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
            term_to_genes: Dict[str, set] = {}
            term_to_desc: Dict[str, str] = {}
            with open(gene_sets_input, encoding="utf-8") as f:
                for line in f:
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
            raise ValueError(f"Failed to load '{gene_sets_input}': {e}") from e
    raise ValueError("gene_sets must be dict, GMT path or gseapy library name")


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
    background: Optional[Iterable[Any]] = None,
    pval_cutoff: float = 0.05,
    min_size: int = 5,
    max_size: int = 500,
    restrict_background_to_gene_sets: bool = True,
    return_all: bool = False,
    verbose: bool = True,
    organism: str = "mouse",
    gene_case: Optional[str] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    Hypergeometric over-representation analysis.
    """
    genes = _clean_gene_list(gene_list, gene_case=gene_case)
    if not genes:
        if verbose:
            _log_info("gene_list is empty")
        return pd.DataFrame(columns=ORA_COLUMNS)
    if min_size < 1 or max_size < min_size or not (0 <= pval_cutoff <= 1):
        raise ValueError("Invalid min_size, max_size or pval_cutoff")
    term_to_genes, term_to_desc = _load_gene_sets(
        gene_sets, organism=organism, verbose=verbose, gene_case=gene_case
    )
    all_gs_genes = set().union(*term_to_genes.values()) if term_to_genes else set()
    if background is not None:
        if isinstance(background, str) and background.lower() == "all":
            universe = all_gs_genes
        else:
            bg = set(_clean_gene_list(background, gene_case=gene_case))
            universe = bg & all_gs_genes if restrict_background_to_gene_sets else bg
    else:
        universe = all_gs_genes
    N = len(universe)
    if N == 0:
        if verbose:
            _log_info("Universe is empty")
        return pd.DataFrame(columns=ORA_COLUMNS)
    genes_in_universe = [g for g in genes if g in universe]
    n = len(genes_in_universe)
    if verbose:
        _log_info(f"Input genes: {len(genes)}, mapped: {n}, universe: {N}")
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
    res_df = res_df.sort_values("p.adjust").reset_index(drop=True)
    res_df.attrs.update({"method": "ora", "organism": organism, "gene_case": gene_case})
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
    background: Optional[Iterable[Any]] = None,
    pval_cutoff: float = 0.05,
    min_size: int = 5,
    max_size: int = 500,
    return_all: bool = False,
    verbose: bool = True,
    gene_case: Optional[str] = None,
    kegg_library: Optional[str] = None,
    **kwargs: Any,
) -> pd.DataFrame:
    """
    KEGG pathway enrichment (wrapper around run_enrichment).
    """
    org_map = {"mouse": "Mouse", "mmu": "Mouse", "human": "Human", "hsa": "Human"}
    gseapy_org = org_map.get(str(organism).lower())
    if gseapy_org is None:
        raise ValueError(f"Unsupported organism '{organism}' for run_kegg")
    kegg_lib_map = {"Mouse": "KEGG_2026", "Human": "KEGG_2026"}
    if kegg_library is None:
        kegg_library = kegg_lib_map[gseapy_org]
    return run_enrichment(
        gene_list=gene_list,
        gene_sets=kegg_library,
        background=background,
        pval_cutoff=pval_cutoff,
        min_size=min_size,
        max_size=max_size,
        return_all=return_all,
        verbose=verbose,
        organism=gseapy_org,
        gene_case=gene_case,
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
            current = set(g.strip() for g in re.split(r"[;,]+", genes_str) if g.strip())
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
