"""
Gene feature handling and bias correction utilities.

This module provides functions to attach per-gene length and intron number
information (used for bias correction) and to generate such tables from GTF
annotation files. Package data access uses importlib.resources for robustness
across installation methods.

Note on precision:
- gene_length uses full exon interval union per gene.
- intron_number uses a max-exons-per-transcript heuristic (see generate_gene_features_from_gtf).
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# --- importlib.resources compatibility (py>=3.10 stdlib, else backport)
# Using the backport on 3.9 avoids spec.origin=None issues in some install scenarios (editable/wheel in CI).
if sys.version_info >= (3, 10):
    from importlib.resources import as_file, files
else:
    from importlib_resources import as_file, files


@contextmanager
def _open_package_data(filename: str) -> Iterator[Path]:
    """Yield a real filesystem Path for a read-only file inside scatrans/data/.

    This is the recommended way to access package data resources and is
    safe for installed packages (wheel/sdist) and editable installs.
    Raises FileNotFoundError (wrapped) if the resource does not exist.
    """
    ref = files("scatrans.data") / filename
    with as_file(ref) as concrete:
        p = Path(concrete)
        if not p.exists():
            raise FileNotFoundError(f"scatrans/data/{filename} not found in package resources")
        yield p


def list_available_gene_features(verbose: bool = False) -> list[str]:
    """
    List all available gene feature parquet files shipped with the package.

    Parameters
    ----------
    verbose : bool, default False
        If True, also log the list at INFO level.

    Returns
    -------
    list of str
        Filenames of .parquet files available via the package data.
    """
    # Known files (used as fallback when discovery fails)
    known = [
        "mouse_2020A_gene_features.parquet",
        "Mus_musculus.GRCm39.115_gene_features.parquet",
        "human_GRCh38_2024A_gene_features.parquet",
    ]

    discovered: list[str] = []
    try:
        data_traversable = files("scatrans.data")
        # iterdir works on Traversable (importlib.resources)
        for item in data_traversable.iterdir():
            name = getattr(item, "name", str(item))
            if name.endswith(".parquet"):
                discovered.append(name)
    except Exception as exc:
        # Best-effort fallback for unusual environments
        logger.debug("list_available_gene_features discovery failed (best-effort): %s", exc)
        discovered = []

    files_list = sorted(set(discovered)) if discovered else known

    if verbose:
        logger.info("Available gene feature files in package data:")
        for f in files_list:
            logger.info(f"   • {f}")
    return files_list


def generate_gene_features_from_gtf(
    gtf_path: str,
    output_name: str | None = None,
    organism: str = "mouse",
):
    """
    Generate a gene features parquet (gene_length + intron_number) from a GTF.

    This is the recommended way to create custom tables for human or non-standard
    annotations so that bias correction in active_score() can be used.

    **gene_length**: proper union of all exon intervals of the gene (fixed from an
    earlier "sum all transcripts" bug).

    **intron_number**: heuristic = (max number of exons among any transcript of the gene) - 1.
    This is only an approximation:
    - Uses the "most exon-rich" transcript as a proxy for the gene.
    - Does not perform coordinate-level union of intron intervals across isoforms.
    - The main/canonical transcript may differ from the max-exon one.
    See source comments in generate_gene_features_from_gtf for details.

    Requires the `gtfparse` package (installed via `pip install "scatrans[gene_features]"`).

    After generation, attach the result with:
        adata = scat.add_gene_features(adata, gene_features_path="your_output.parquet")

    Long-running step (GTF parsing). Progress is emitted via the 'scatrans' logger.
    When used from the CLI entrypoint, the CLI configures logging so messages appear.

    Parameters
    ----------
    gtf_path : str
        Path to a 10X Genomics or GENCODE genes.gtf file (must contain "exon" and "gene" features).
    output_name : str, optional
        Output parquet filename. If None, auto-generated as {organism}_gene_features.parquet.
    organism : str, default "mouse"
        Only used for default output naming; the generated table itself is generic.
    """
    try:
        import gtfparse
    except ImportError:
        raise ImportError(
            "gtfparse is required to generate gene features. "
            "Install with: pip install 'scatrans[gene_features]' or pip install gtfparse"
        ) from None

    gtf_path = Path(gtf_path).expanduser()
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_path}")

    if output_name is None:
        output_name = f"{organism}_gene_features.parquet"
    else:
        output_name = str(Path(output_name).expanduser())

    logger.info("Parsing GTF file (may take 30-60 seconds)... %s", gtf_path)
    df = gtfparse.read_gtf(gtf_path)

    if hasattr(df, "to_pandas"):
        df = df.to_pandas()
        logger.info("Converted to Pandas DataFrame")

    # 1. gene_length: proper union of exon intervals per gene (critical!)
    # Previous code summed lengths across all transcripts → massive overcount for multi-isoform genes.
    exon = df[df["feature"] == "exon"].copy()
    if exon.empty:
        logger.warning(
            "No exon features found in GTF; gene_length will be 0 for all genes. "
            "Bias correction will fall back to median if features are attached."
        )
    exon["start"] = pd.to_numeric(exon["start"], errors="coerce")
    exon["end"] = pd.to_numeric(exon["end"], errors="coerce")
    exon = exon.dropna(subset=["gene_id", "start", "end"])
    if exon.empty:
        logger.warning(
            "No valid exon intervals (missing/NaN coordinates); gene_length will be 0 for all genes."
        )

    def _exon_union_length(grp):
        if len(grp) == 0:
            return 0
        # Collect (start, end) pairs (GTF is 1-based, inclusive)
        ivs = sorted(
            (int(s), int(e))
            for s, e in zip(grp["start"], grp["end"])
            if pd.notna(s) and pd.notna(e)
        )
        if not ivs:
            return 0
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            last = merged[-1]
            if s <= last[1] + 1:  # overlap or adjacent (inclusive)
                last[1] = max(last[1], e)
            else:
                merged.append([s, e])
        return sum(e - s + 1 for s, e in merged)

    if exon.empty:
        gene_length = pd.Series(dtype=float, name="gene_length")
    else:
        gene_length = (
            exon.groupby("gene_id")
            .apply(_exon_union_length, include_groups=False)
            .rename("gene_length")
        )

    # 2. intron_number: a *proxy* using the transcript with the largest exon count.
    #    intron_number ≈ max_exons_over_transcripts - 1
    #
    #    Limitations (documented for users of bias correction):
    #    - "Longest" by exon count is only a heuristic; the canonical/main isoform may have
    #      fewer exons or different structure.
    #    - Unlike gene_length (which now does a proper interval *union* across *all* exons
    #      of the gene), this does *not* attempt a cross-transcript "union" of intron regions
    #      (that would require coordinate-aware intron interval calculation).
    #    - Therefore gene_length and intron_number have different levels of precision.
    #    - This value is only used as a numeric feature for Huber regression bias correction
    #      of velocity delta. Small errors are usually tolerated by the robust fit.
    if exon.empty:
        intron_number = pd.Series(dtype=float, name="intron_number")
    else:
        transcript_exons = exon.groupby(["gene_id", "transcript_id"]).size().rename("exon_count")
        intron_number = (
            (transcript_exons.groupby("gene_id").max() - 1).clip(lower=0).rename("intron_number")
        )

    # 3. Gene info - handle both GENCODE ('gene_type') and Ensembl ('gene_biotype')
    gene_cols = ["gene_id", "gene_name"]
    gene_type_col = None

    if "gene_type" in df.columns:
        gene_type_col = "gene_type"
    elif "gene_biotype" in df.columns:
        gene_type_col = "gene_biotype"
        logger.info(
            "Using 'gene_biotype' column (Ensembl-style GTF) and renaming it to 'gene_type' for consistency."
        )

    if gene_type_col:
        gene_cols.append(gene_type_col)

    gene_info = df[df["feature"] == "gene"][gene_cols].drop_duplicates("gene_id")

    # Rename gene_biotype → gene_type if needed (for downstream consistency)
    if gene_type_col == "gene_biotype":
        gene_info = gene_info.rename(columns={"gene_biotype": "gene_type"})

    # 4. Merge
    gene_df = gene_info.set_index("gene_id").join(gene_length).join(intron_number).reset_index()

    # Ensure 'gene_type' column exists (even if empty)
    if "gene_type" not in gene_df.columns:
        gene_df["gene_type"] = np.nan

    gene_df = gene_df[["gene_id", "gene_name", "gene_length", "intron_number", "gene_type"]]
    gene_df["gene_length"] = pd.to_numeric(gene_df["gene_length"], errors="coerce").fillna(0)
    gene_df["intron_number"] = pd.to_numeric(gene_df["intron_number"], errors="coerce").fillna(0)
    if (gene_df["gene_length"] == 0).all():
        logger.warning(
            "All gene_length values are zero (no usable exon intervals). "
            "Downstream bias correction will use median fallback for affected genes."
        )

    n_before_names = len(gene_df)
    gene_df = gene_df.dropna(subset=["gene_name"])
    n_missing_names = n_before_names - len(gene_df)
    if n_missing_names:
        logger.warning(
            "Dropped %d genes with missing gene_name (common in some Ensembl GTF exports).",
            n_missing_names,
        )

    dup_mask = gene_df.duplicated(subset=["gene_name"], keep="first")
    n_dup_names = int(dup_mask.sum())
    if n_dup_names:
        logger.warning(
            "Dropped %d duplicate gene_name entries (multi-isoform GTFs); kept first per name. "
            "Bias correction joins on gene_name — ensure adata.var_names are unique symbols.",
            n_dup_names,
        )
        gene_df = gene_df[~dup_mask]

    logger.info("Processing completed! %d genes processed", len(gene_df))
    logger.debug("%s", gene_df.head())

    gene_df.to_parquet(output_name, index=False, compression="zstd")
    size_mb = Path(output_name).stat().st_size / (1024 * 1024)
    logger.info("Parquet generated → %s (%.1f MB)", output_name, size_mb)
    return gene_df


def add_gene_features(
    adata,
    organism: str = "mouse",
    gene_feature_file: str | None = None,
    gene_features_path: str | None = None,
):
    """
    Add gene features (length, intron number) to adata.var for bias correction.

    The function looks for files in three ways (priority order):
    1. gene_features_path (full user-provided path; ~ expanded)
    2. gene_feature_file (bare filename inside package data/; path-like values auto-treated as custom)
    3. Default filename chosen from `organism`

    Uses robust importlib.resources access so it works whether the package
    is installed from wheel, sdist, or in editable mode.

    Parameters
    ----------
    adata : AnnData
    organism : str, default "mouse"
        Used only as fallback when gene_feature_file is not provided.
    gene_feature_file : str, optional
        Bare filename inside the package data/ directory (e.g. one of the files from
        list_available_gene_features()). If you pass a path containing / or ~ here,
        it will be treated as a custom file (like gene_features_path) for convenience.
    gene_features_path : str, optional
        Full custom path to a parquet file (highest priority). ~ is expanded to $HOME.
    """
    logger.info("Loading gene features for bias correction...")

    final_path: Path | None = None
    using_package_data = False

    # Priority 1: Full custom path (user file, do not touch package resources)
    if gene_features_path is not None:
        final_path = Path(gene_features_path).expanduser()
        logger.info("Using custom path: %s", final_path)

    # Priority 2: explicit filename inside package data/
    elif gene_feature_file is not None:
        gf_val = str(gene_feature_file)
        p = Path(gf_val)
        # Detect if user passed a path (~/..., /abs, rel/dir, or with /) to the "filename" arg by mistake.
        # In that case treat it as custom path for convenience (and to avoid prepending package data dir).
        if (
            p.is_absolute()
            or gf_val.startswith("~")
            or "/" in gf_val
            or "\\" in gf_val
            or len(p.parts) > 1
        ):
            final_path = p.expanduser()
            logger.warning(
                "gene_feature_file received a path-like value (%s). "
                "Treating as custom file path (prefer the gene_features_path= parameter for this).",
                gf_val,
            )
        else:
            using_package_data = True
            pkg_filename = gf_val

    # Priority 3: Default based on organism (backward compatible)
    else:
        using_package_data = True
        avail = list_available_gene_features(verbose=False)
        organism_norm = str(organism).lower()
        if organism_norm in ("human", "hs", "hsa"):
            # Prefer any human-named file present in the package data
            human_cands = [
                f
                for f in avail
                if "human" in f.lower() or "grch" in f.lower() or "hg38" in f.lower()
            ]
            pkg_filename = (
                human_cands[0] if human_cands else "human_GRCh38_2024A_gene_features.parquet"
            )
        elif organism_norm in ("mouse", "mm", "mmu"):
            mouse_cands = [f for f in avail if "mouse" in f.lower() or f.startswith("Mus")]
            pkg_filename = mouse_cands[0] if mouse_cands else "mouse_2020A_gene_features.parquet"
        else:
            raise ValueError(
                f"Unsupported organism '{organism}' for add_gene_features. "
                "Use 'human'/'hs'/'hsa' or 'mouse'/'mm'/'mmu'."
            )
        logger.info("Using default for %s: %s", organism, pkg_filename)

    # Load the gene features table.
    # CRITICAL for wheel installs: when using_package_data, pd.read_parquet MUST be
    # performed inside the _open_package_data() context. as_file() extracts to a
    # temp path that is deleted on __exit__ of the with block.
    gf = None
    if using_package_data:
        try:
            with _open_package_data(pkg_filename) as resolved:
                logger.info("Resolved package data file: %s", resolved)
                gf = pd.read_parquet(resolved).set_index("gene_name")
        except Exception as exc:
            available = list_available_gene_features(verbose=False)
            raise FileNotFoundError(
                f"Gene features file not found inside package data: {pkg_filename}\n"
                f"Available files: {available}\n\n"
                "Solutions:\n"
                "  1. Use the CLI to generate it:\n"
                f"     generate-gene-features --gtf /path/to/genes.gtf --output {pkg_filename}\n"
                "  2. Provide your own file: add_gene_features(adata, gene_features_path='your_file.parquet')\n"
                f"  3. Specify filename in package data: add_gene_features(adata, gene_feature_file='{pkg_filename}')"
            ) from exc
    else:
        if final_path is None:
            raise RuntimeError(
                "Internal error resolving gene features path. "
                "This should not happen; please report with your add_gene_features call."
            )
        gf = pd.read_parquet(final_path).set_index("gene_name")

    gf = gf[gf.index.notna()]
    gf = gf[gf.index.astype(str).str.len() > 0]
    gf = gf[~gf.index.duplicated(keep="first")]

    adata.var["gene_length"] = gf["gene_length"].reindex(adata.var_names)
    adata.var["intron_number"] = gf["intron_number"].reindex(adata.var_names)

    valid_count = int(adata.var["gene_length"].notna().sum())
    logger.info("Successfully mapped features for %d out of %d genes.", valid_count, adata.n_vars)

    return adata
