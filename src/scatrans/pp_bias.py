"""Bias correction preprocessing for scATrans.

Includes functions to attach gene-level features (length, intron count)
needed for Huber regression bias correction in active_score().
"""

from __future__ import annotations
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Union
import warnings
import logging

logger = logging.getLogger(__name__)


def add_gene_features(
    adata,
    gene_feature_file: Optional[Union[str, Path]] = None,
    gene_length_col: str = "gene_length",
    intron_number_col: str = "intron_number",
    gene_id_col: str = None,
    fillna: float = np.nan,
    verbose: bool = True,
) -> "anndata.AnnData":
    """
    Attach precomputed gene features (length + intron number) to adata.var
    for bias correction in `active_score()`.

    The features are used by HuberRegressor to correct velocity residuals
    for technical biases (longer genes / genes with more introns tend to
    show different unspliced/spliced ratios).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix. Will be modified in-place (var columns added).
    gene_feature_file : str or Path, optional
        Path to parquet/CSV/TSV file with gene features.
        Expected columns: at least 'gene_length' and 'intron_number'.
        Index or a column should match adata.var_names (gene symbols or Ensembl IDs).
        If None, the function looks for common locations relative to package
        (data/mouse_2020A_gene_features.parquet) or does nothing.
    gene_length_col : str
        Column name for gene length in the feature file (default: "gene_length").
    intron_number_col : str
        Column name for intron count (default: "intron_number").
    gene_id_col : str, optional
        If features file uses a column for gene IDs instead of index, specify it here.
    fillna : float
        Value to fill missing genes (default NaN). Genes without features
        will have NaN and be skipped in bias correction (valid_feat mask).
    verbose : bool
        Print matching stats.

    Returns
    -------
    adata : AnnData
        The same object with added columns in .var
    """
    import anndata as ad  # local import to avoid hard dep at top

    if gene_feature_file is None:
        # Try to find bundled mouse features
        pkg_dir = Path(__file__).parent.parent.parent  # rough guess
        candidate = pkg_dir / "data" / "mouse_2020A_gene_features.parquet"
        if candidate.exists():
            gene_feature_file = candidate
            if verbose:
                logger.info(f"Using bundled gene features: {gene_feature_file}")
        else:
            warnings.warn(
                "No gene_feature_file provided and no bundled features found. "
                "Bias correction will be skipped (gene_length/intron_number = NaN). "
                "You can generate features with generate_gene_features_from_gtf() "
                "or provide a precomputed table.",
                UserWarning
            )
            adata.var[gene_length_col] = fillna
            adata.var[intron_number_col] = fillna
            return adata

    gene_feature_file = Path(gene_feature_file)
    if not gene_feature_file.exists():
        raise FileNotFoundError(f"gene_feature_file not found: {gene_feature_file}")

    # Load features
    if gene_feature_file.suffix == ".parquet":
        try:
            feat = pd.read_parquet(gene_feature_file)
        except Exception as e:
            raise ImportError(
                "pyarrow or fastparquet is required to read .parquet files. "
                "Install with: pip install pyarrow"
            ) from e
    elif gene_feature_file.suffix in {".csv", ".tsv", ".txt"}:
        sep = "\t" if gene_feature_file.suffix in {".tsv", ".txt"} else ","
        feat = pd.read_csv(gene_feature_file, sep=sep, index_col=0)
    else:
        feat = pd.read_csv(gene_feature_file, index_col=0)

    if gene_id_col is not None and gene_id_col in feat.columns:
        feat = feat.set_index(gene_id_col)

    # Standardize
    feat = feat.rename(columns={
        gene_length_col: "gene_length",
        intron_number_col: "intron_number"
    })

    if "gene_length" not in feat.columns or "intron_number" not in feat.columns:
        raise ValueError(
            f"Feature file must contain '{gene_length_col}' and '{intron_number_col}' "
            f"(after optional rename). Found columns: {feat.columns.tolist()}"
        )

    # Match to adata.var_names
    common_genes = adata.var_names.intersection(feat.index)
    n_match = len(common_genes)
    n_total = adata.n_vars

    if verbose:
        logger.info(
            f"Gene feature matching: {n_match}/{n_total} genes "
            f"({100*n_match/n_total:.1f}%) found in feature table."
        )

    if n_match == 0:
        warnings.warn(
            "No genes matched between adata.var_names and the feature file index. "
            "Check that gene identifiers are consistent (symbol vs Ensembl ID). "
            "Bias correction will use NaNs.",
            UserWarning
        )

    # Assign
    adata.var["gene_length"] = feat.reindex(adata.var_names)["gene_length"].fillna(fillna).values
    adata.var["intron_number"] = feat.reindex(adata.var_names)["intron_number"].fillna(fillna).values

    # Ensure numeric
    adata.var["gene_length"] = pd.to_numeric(adata.var["gene_length"], errors="coerce")
    adata.var["intron_number"] = pd.to_numeric(adata.var["intron_number"], errors="coerce")

    return adata


def generate_gene_features_from_gtf(
    gtf_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
    feature_type: str = "gene",
    attribute_gene_id: str = "gene_id",
    attribute_gene_name: str = "gene_name",
    attribute_gene_biotype: str = "gene_biotype",
) -> pd.DataFrame:
    """
    (Placeholder) Generate gene_length + intron_number table from a GTF file.

    This is a stub. Full implementation requires a GTF parser
    (e.g. gtfparse, pyranges, or pyensembl) and exon merging logic
    to count introns per gene and compute genomic span as length.

    For mouse, we recommend using the precomputed
    `data/mouse_2020A_gene_features.parquet` that ships with scATrans.

    Parameters
    ----------
    gtf_path : str or Path
        Path to .gtf or .gtf.gz
    output_path : str or Path, optional
        If provided, save the resulting DataFrame as parquet.
    ... other params ...

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by gene_name or gene_id with columns
        ['gene_length', 'intron_number', ...]
    """
    raise NotImplementedError(
        "generate_gene_features_from_gtf is not yet implemented in this release. "
        "Please use the pre-bundled mouse_2020A_gene_features.parquet "
        "or implement exon merging logic yourself (common in velocyto/ scVelo pipelines)."
    )
