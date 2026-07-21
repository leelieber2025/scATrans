# Installation

## Requirements

- Python ≥ 3.10
- A scientific Python environment (NumPy, SciPy, pandas, AnnData, scanpy)

## Bioconda

```bash
conda install -c conda-forge -c bioconda scatrans
```

This installs core dependencies. Optional extras listed below are distributed as
PyPI extras; install them with `pip`, or add the corresponding packages to the
conda environment.

## PyPI

```bash
pip install scatrans

# Optional extras
pip install "scatrans[advanced,gene_features]"  # scVelo mode + GTF feature CLI
pip install "scatrans[pseudobulk]"              # PyDESeq2
pip install "scatrans[memento]"                 # Memento DE backend
pip install "scatrans[gsea]"                    # GSEA (gseapy)
```

Bundled gene-feature tables (length, intron number) for mouse and human support
optional Huber bias correction. Custom GTF-derived tables:
{doc}`user_guide/gene_features`.

## Development install

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev]"
```

## Versioning

The package version is defined in `src/scatrans/_version.py` (`__version__`).
Runtime `scatrans.__version__`, packaging metadata, and documentation release
strings all read that value. For a release: bump `__version__`, update
`CHANGELOG.md`, then run `python -m build` or
`python scripts/make_release_zips.py`.

## Logging

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

## Nascent-layer quality check

Before mechanism analysis, inspect the global unspliced fraction:

```python
import scatrans as scat

ufrac = scat.qc.unspliced_global(adata)
regime = scat.qc.regime_diagnosis(adata)
print(regime["regime"], regime["reliability"], regime["message"])
# regime in {"ok", "low_unspliced", "high_unspliced"}; reliability in [0, 1]
```

`partition_de_by_mechanism` always runs this pre-flight
(`result.regime` / `meta["regime"]`, fail-soft) and scales mechanism confidence
by `reliability`. `run_default_pipeline` stores the same block and applies the
scale when `annotate_mechanism=True`.
