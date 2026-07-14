# Installation

## From Bioconda

The core package is available on [Bioconda](https://anaconda.org/bioconda/scatrans)
(noarch):

```bash
conda install -c conda-forge -c bioconda scatrans
```

This installs the core dependencies only. The optional extras below (advanced
mode, pseudobulk, Memento, GSEA) are packaged as PyPI extras — install them
with `pip` as shown, or add the corresponding packages to your conda
environment manually.

## From PyPI

```bash
# Basic installation
pip install scatrans

# With support for scVelo-based advanced mode and the gene feature generation CLI
pip install "scatrans[advanced,gene_features]"

# With support for pseudobulk differential expression using PyDESeq2
pip install "scatrans[pseudobulk]"

# Optional: Memento (Cell 2024) as an additional cell-level DE backend
pip install "scatrans[memento]"

# Optional: GSEA (pulls in gseapy)
pip install "scatrans[gsea]"
```

The package ships precomputed gene feature tables (gene length + intron
number) for both mouse and human. These are used for optional bias
correction in `active_score`. You can also supply custom tables (e.g. from
your own GTF) — see {doc}`user_guide/gene_features`.

## Install from source

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev]"
```

## Versioning (developers)

The package version is defined in a **single place**:
`src/scatrans/_version.py` (`__version__`). Packaging, runtime
`scatrans.__version__`, and documentation release metadata all read that
value. To prepare a release, bump `__version__` there, update `CHANGELOG.md`,
then run `python -m build` or `python scripts/make_release_zips.py`.
## Logging

The package logs under the name `scatrans`. You can control verbosity with:

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

## Quick data quality check

Before analysis, inspect the global unspliced fraction:

```python
import scatrans as scat
ufrac = scat.qc.unspliced_global(adata)   # logs INFO + WARNING if > 50%
```

`active_score` automatically runs this check and records the value in
diagnostics.
