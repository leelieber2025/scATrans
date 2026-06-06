"""
scATrans pp_bias.py

Gene feature handling and bias correction utilities:
- generate_gene_features_from_gtf() : build gene_length + intron_number table from GTF
- add_gene_features() : attach features to AnnData.var for use in active_score()
- list_available_gene_features() : helper to discover bundled parquet files
"""

import os
import logging
from pathlib import Path
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def list_available_gene_features():
    """
    List all available gene feature parquet files in the package data directory.

    Returns
    -------
    list of str
        Filenames of all .parquet files in src/scatrans/data/
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(current_dir, "data")

    if not os.path.exists(data_dir):
        logger.warning("⚠️ Data directory not found.")
        return []

    files = [f for f in os.listdir(data_dir) if f.endswith(".parquet")]
    if not files:
        logger.warning("⚠️ No gene feature files found in data/ directory.")
        return []

    logger.info("📋 Available gene feature files:")
    for f in sorted(files):
        print(f"   • {f}")
    return files


def generate_gene_features_from_gtf(
    gtf_path: str,
    output_name: str = None,
    organism: str = "mouse"
):
    """
    Generate gene features parquet from a GTF file (for developers/maintainers and CLI).

    Parameters
    ----------
    gtf_path : str
        Path to the 10X Genomics or GENCODE genes.gtf file.
    output_name : str, optional
        Output parquet filename. If None, auto-generated as {organism}_gene_features.parquet
    organism : str, default "mouse"
        Used only for default naming (you can use any name you want).
    """
    try:
        import gtfparse
    except ImportError:
        raise ImportError(
            "gtfparse is required to generate gene features. "
            "Install with: pip install 'scatrans[gene_features]' or pip install gtfparse"
        )

    if output_name is None:
        output_name = f"{organism}_gene_features.parquet"

    print(f"🚀 Parsing GTF file (may take 30-60 seconds)... {gtf_path}")
    df = gtfparse.read_gtf(gtf_path)

    if hasattr(df, "to_pandas"):
        df = df.to_pandas()
        print("✅ Converted to Pandas DataFrame")

    # 1. gene_length (sum of all exons)
    exon = df[df['feature'] == 'exon'].copy()
    exon['length'] = exon['end'] - exon['start'] + 1
    gene_length = exon.groupby('gene_id')['length'].sum().rename('gene_length')

    # 2. intron_number (max exons per transcript - 1)
    transcript_exons = exon.groupby(['gene_id', 'transcript_id']).size().rename('exon_count')
    intron_number = (transcript_exons.groupby('gene_id').max() - 1).clip(lower=0).rename('intron_number')

    # 3. Gene info - handle both GENCODE ('gene_type') and Ensembl ('gene_biotype')
    gene_cols = ['gene_id', 'gene_name']
    gene_type_col = None

    if 'gene_type' in df.columns:
        gene_type_col = 'gene_type'
    elif 'gene_biotype' in df.columns:
        gene_type_col = 'gene_biotype'
        print("ℹ️  Using 'gene_biotype' column (Ensembl-style GTF) and renaming it to 'gene_type' for consistency.")

    if gene_type_col:
        gene_cols.append(gene_type_col)

    gene_info = df[df['feature'] == 'gene'][gene_cols].drop_duplicates('gene_id')

    # Rename gene_biotype → gene_type if needed (for downstream consistency)
    if gene_type_col == 'gene_biotype':
        gene_info = gene_info.rename(columns={'gene_biotype': 'gene_type'})

    # 4. Merge
    gene_df = gene_info.set_index('gene_id').join(gene_length).join(intron_number).reset_index()

    # Ensure 'gene_type' column exists (even if empty)
    if 'gene_type' not in gene_df.columns:
        gene_df['gene_type'] = np.nan

    gene_df = gene_df[['gene_id', 'gene_name', 'gene_length', 'intron_number', 'gene_type']]
    gene_df = gene_df.dropna(subset=['gene_length'])

    print(f"✅ Processing completed! {len(gene_df):,} genes processed")
    print(gene_df.head())

    gene_df.to_parquet(output_name, index=False, compression='zstd')
    size_mb = Path(output_name).stat().st_size / (1024 * 1024)
    print(f"🎉 Parquet generated → {output_name} ({size_mb:.1f} MB)")
    return gene_df


def add_gene_features(adata, organism="mouse", gene_feature_file=None, gene_features_path=None):
    """
    Add gene features (length, intron number) to adata.var for bias correction.

    Flexible selection based on actual files in src/scatrans/data/

    Parameters
    ----------
    adata : AnnData
    organism : str, default "mouse"
        Used only as fallback when gene_feature_file is not provided.
    gene_feature_file : str, optional
        Filename in the data/ directory (e.g. "mouse_2020A_gene_features.parquet"
        or "human_gencode_v49.parquet"). Highest priority after full path.
    gene_features_path : str, optional
        Full custom path to a parquet file (highest priority).
    """
    print("🧬 Loading gene features for bias correction...")

    # Priority 1: Full custom path
    if gene_features_path is not None:
        final_path = gene_features_path
        print(f"   Using custom path: {final_path}")

    # Priority 2: Filename in package data/
    elif gene_feature_file is not None:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        final_path = os.path.join(current_dir, "data", gene_feature_file)
        print(f"   Using specified feature file: {gene_feature_file}")

    # Priority 3: Default based on organism (backward compatible)
    else:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        default_filename = f"{organism}_gene_features.parquet" if organism == "human" else "mouse_2020A_gene_features.parquet"
        final_path = os.path.join(current_dir, "data", default_filename)
        print(f"   Using default for {organism}: {default_filename}")

    # Check file exists
    if not os.path.exists(final_path):
        available = list_available_gene_features()
        raise FileNotFoundError(
            f"Gene features file not found: {final_path}\n"
            f"Available files in package data/: {available}\n\n"
            f"💡 Solutions:\n"
            f"  1. Use the CLI to generate it:\n"
            f"     generate-gene-features --gtf /path/to/genes.gtf --output {os.path.basename(final_path)}\n"
            f"  2. Provide your own file: add_gene_features(adata, gene_features_path='your_file.parquet')\n"
            f"  3. Specify filename in package data: add_gene_features(adata, gene_feature_file='mouse_2020A_gene_features.parquet')"
        )

    try:
        gf = pd.read_parquet(final_path).set_index("gene_name")
        gf = gf[~gf.index.duplicated(keep="first")]

        adata.var["gene_length"] = gf["gene_length"].reindex(adata.var_names)
        adata.var["intron_number"] = gf["intron_number"].reindex(adata.var_names)

        valid_count = adata.var["gene_length"].notna().sum()
        print(f"✅ Successfully mapped features for {valid_count} out of {adata.n_vars} genes.")
    except Exception as e:
        print(f"⚠️ Failed to load gene features ({e}). Continuing with NaN values.")
        adata.var["gene_length"] = np.nan
        adata.var["intron_number"] = np.nan

    return adata
