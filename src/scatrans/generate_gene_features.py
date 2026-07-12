#!/usr/bin/env python3
"""
scATrans Gene Features Generator (CLI)

Generates gene feature tables (gene_length + intron_number) from a GTF
for use with scat.add_gene_features(..., gene_features_path=...) and
subsequent bias correction in active_score().

Usage (after `pip install "scatrans[gene_features]"`):
    generate-gene-features --gtf /path/to/genes.gtf --output my_features.parquet

Or from source:
    python -m scatrans.generate_gene_features --gtf ...

See README for full human/custom annotation workflow.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

try:
    from .pp_bias import generate_gene_features_from_gtf
except ImportError:  # pragma: no cover - direct script execution before install
    # Fallback when run as a loose script (no package context).
    from pp_bias import generate_gene_features_from_gtf  # type: ignore[no-redef]


def main():
    # Configure logging ONLY inside the CLI entry point.
    # Never at import time: "import scatrans" must not call basicConfig or affect
    # the caller's root logger (no bare messages, no forced INFO level).
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(level=logging.INFO, format="%(message)s")
    parser = argparse.ArgumentParser(
        description="Generate scATrans gene features parquet from a GTF annotation file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  generate-gene-features --gtf genes.gtf --output mouse_2020A_gene_features.parquet
  generate-gene-features --gtf gencode.v49.primary_assembly.annotation.gtf --organism human --output human_GRCh38_gene_features.parquet

Typical workflow (end users):
  1. Generate your table from a 10x or GENCODE genes.gtf
  2. Use it directly:
       adata = scat.add_gene_features(adata, gene_features_path="human_GRCh38_gene_features.parquet")
  3. Then run: adata_res, sig, results = scat.active_score(adata, ...)

Bundling into the package (advanced):
  Copy the .parquet into src/scatrans/data/ and reinstall in editable mode.
        """,
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

    output_for_log = str(Path(args.output).expanduser())
    logging.info("Starting gene features generation from: %s", gtf_path)
    logging.info("   Output will be written to: %s", output_for_log)
    logging.info("   Organism label: %s", args.organism)

    try:
        df = generate_gene_features_from_gtf(
            gtf_path=str(gtf_path), output_name=args.output, organism=args.organism
        )
        logging.info("\nGene features successfully generated!")
        logging.info("   File: %s", args.output)
        logging.info("   Genes processed: %s", f"{len(df):,}")
        logging.info("\nNext steps:")
        logging.info(
            "   • Use in your analysis:\n"
            "       import scatrans as scat\n"
            "       adata = scat.add_gene_features(adata, gene_features_path='%s')\n"
            "       adata_res, significant, results = scat.active_score(adata, ...)",
            args.output,
        )
        logging.info(
            "   • (Advanced) To ship it inside scatrans: copy to src/scatrans/data/ and `pip install -e '.[gene_features]'`"
        )
    except Exception as e:
        print(f"Generation failed: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
