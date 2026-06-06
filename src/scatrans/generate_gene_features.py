#!/usr/bin/env python3
"""
scATrans Gene Features Generator

Command-line interface to generate gene feature tables (gene_length + intron_number)
from a GTF file for use with scATrans bias correction.

Usage (after pip install):
    generate-gene-features --gtf /path/to/genes.gtf --output mouse_gene_features.parquet

Or from source:
    python -m scatrans.generate_gene_features --gtf ...

This script is also available as a console entry point after installation.
"""
import argparse
from pathlib import Path
import sys

try:
    from .pp_bias import generate_gene_features_from_gtf
except ImportError:
    # Fallback for direct script execution before installation
    from pp_bias import generate_gene_features_from_gtf


def main():
    parser = argparse.ArgumentParser(
        description="Generate scATrans gene features parquet from a GTF annotation file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  generate-gene-features --gtf genes.gtf --output mouse_2020A_gene_features.parquet
  generate-gene-features --gtf gencode.vM32.primary_assembly.annotation.gtf --organism mouse --output gencode_vM32_features.parquet

After generation:
  1. (Optional) Copy the .parquet to src/scatrans/data/ to bundle it with the package
  2. Re-install the package in editable mode: pip install -e ".[gene_features]"
        """
    )
    parser.add_argument(
        "--gtf", required=True,
        help="Path to 10X Genomics or GENCODE genes.gtf file (must contain exon and gene features)"
    )
    parser.add_argument(
        "--output", default="gene_features.parquet",
        help="Output parquet filename (default: gene_features.parquet)"
    )
    parser.add_argument(
        "--organism", default="mouse",
        help="Organism name used only for default naming / metadata (default: mouse)"
    )

    args = parser.parse_args()

    gtf_path = Path(args.gtf)
    if not gtf_path.exists():
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)

    print(f"🚀 Starting gene features generation from: {gtf_path}")
    print(f"   Output will be written to: {args.output}")
    print(f"   Organism label: {args.organism}")

    try:
        df = generate_gene_features_from_gtf(
            gtf_path=str(gtf_path),
            output_name=args.output,
            organism=args.organism
        )
        print("\n✅ Gene features successfully generated!")
        print(f"   File: {args.output}")
        print(f"   Genes processed: {len(df):,}")
        print("\n💡 Next steps:")
        print("   • Use with: adata = scat.add_gene_features(adata, gene_features_path='your_file.parquet')")
        print("   • Or bundle it by copying to src/scatrans/data/ and rebuilding the package")
    except Exception as e:
        print(f"❌ Generation failed: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
