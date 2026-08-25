#!/usr/bin/env python3
"""CLI: build a gene-length / intron-number parquet from a GTF.

generate-gene-features --gtf genes.gtf --output features.parquet
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .pp_bias import generate_gene_features_from_gtf


def main():
    # Configure logging only in the CLI, never at import time.
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(
        description="Generate a gene features parquet from a GTF annotation file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  generate-gene-features --gtf genes.gtf --output features.parquet\n"
            "Then: scat.add_gene_features(adata, gene_features_path='features.parquet')"
        ),
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Path to 10X Genomics or GENCODE genes.gtf file (must contain exon and gene features)",
    )
    parser.add_argument(
        "--output",
        default="gene_features.parquet",
        help="Output parquet filename (default: gene_features.parquet)",
    )
    parser.add_argument(
        "--organism",
        default="mouse",
        help="Organism name used only for default naming / metadata (default: mouse)",
    )

    args = parser.parse_args()

    gtf_path = Path(args.gtf).expanduser()
    if not gtf_path.exists():
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)

    try:
        df = generate_gene_features_from_gtf(
            gtf_path=str(gtf_path), output_name=args.output, organism=args.organism
        )
        logging.info("Wrote %s (%s genes).", args.output, f"{len(df):,}")
    except Exception as e:
        print(f"Generation failed: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
