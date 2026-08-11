# Installation

## Requirements

- Python 3.9+ (CI covers 3.9–3.12)
- A normal scientific stack (NumPy, SciPy, pandas, AnnData, scanpy) — pulled in
  as dependencies

## Install

**PyPI (recommended):**

```bash
pip install scatrans
```

**Bioconda:**

```bash
conda install -c conda-forge -c bioconda scatrans
```

Conda installs the core package. Optional extras below are PyPI extras; add them
with `pip` or install the matching conda packages yourself.

## Optional extras

Install only what you need:

```bash
pip install "scatrans[pseudobulk]"      # PyDESeq2 (replicate-aware DE)
pip install "scatrans[gsea]"           # GSEA via gseapy
pip install "scatrans[memento]"        # Memento DE backend
pip install "scatrans[advanced]"       # scVelo (mode="advanced")
pip install "scatrans[gene_features]"  # gtfparse for custom GTF tables
```

Combine tags as needed, e.g. `"scatrans[pseudobulk,gsea]"`.

Bundled mouse/human gene-feature tables support optional length/intron bias
correction. Custom GTF tables: {doc}`user_guide/gene_features`.

## Check the install

```python
import scatrans as scat
print(scat.__version__)
```

## After install

| Step | Page |
|------|------|
| First analysis | {doc}`quickstart` |
| Which function / backend? | {doc}`faq` |
| Real-data notebooks | {doc}`tutorials/index` |
| Full workflow knobs | {doc}`user_guide/index` |

## Development install

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev]"
```

## Versioning

The single source of truth is `src/scatrans/_version.py`. Runtime
`scatrans.__version__`, packaging metadata, and docs release strings all read
it. For a release: bump `__version__`, update `CHANGELOG.md`, then
`python -m build` or `python scripts/make_release_zips.py`.

## Logging

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

## Quick data check (before mechanism analysis)

```python
import scatrans as scat

print(scat.qc.unspliced_global(adata))
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
# regime: "ok" | "low_unspliced" | "high_unspliced"
```

`partition_de_by_mechanism` always runs this check and stores it as
`result.regime`. Low reliability does not stop the run; it down-weights
mechanism confidence.
