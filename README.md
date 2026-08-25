# scATrans

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![PyPI downloads](https://img.shields.io/pepy/dt/scatrans.svg)](https://pepy.tech/project/scatrans)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Documentation Status](https://readthedocs.org/projects/scatrans/badge/?version=latest)](https://scatrans.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

**scATrans** answers a simple follow-up after differential expression (DE):
among the genes that changed, which look *transcription-driven* and which look
*stabilization-driven*? It uses the nascent (unspliced) residual on top of a
normal DE step—something total-count fold change alone cannot sort out.

| Step | What it does |
|------|----------------|
| DE | Chooses which genes make the list |
| Mechanism | Labels those genes (transcription vs. stabilization) |
| Detection (optional) | Extra nascent-activity scores; does not rewrite mechanism labels |

Pathway- or program-level summaries (`gene_sets=`) are usually more useful than
single-gene labels. Without spliced/unspliced layers, the package still runs
ordinary DE, enrichment, and plots.

Docs: [Read the Docs](https://scatrans.readthedocs.io/en/latest/).

## What you need

- Python 3.9+ (tested 3.9–3.12)
- An AnnData object with a condition column in `.obs`
- For mechanism analysis: `spliced`/`unspliced` (or `mature`/`nascent`) layers
- Nothing yet? `scat.datasets.load_toy()` ships a synthetic example — see [First run](#first-run).

## Install

```bash
pip install scatrans
# or: conda install -c conda-forge -c bioconda scatrans
```

Optional extras (PyDESeq2, GSEA, scVelo, …):
[installation guide](https://scatrans.readthedocs.io/en/latest/installation.html).

## First run

This uses a bundled synthetic AnnData (`load_toy()`). No download. Swap in
your own object after the call works.

```python
import scatrans as scat

adata = scat.datasets.load_toy()  # Control/Disease, sample, spliced/unspliced

result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",          # "human" for human gene symbols
    sample_col="sample",       # toy has this column; pass yours when you have replicates
)
print(result.regime)           # capture quality (reliability in [0, 1])
print(result.selected.head())  # DE genes + mechanism labels
print(result.summary())        # cutoffs used, how many genes were selected
```

What to expect:

- `result.selected` is the DE gene list. Mechanism columns annotate it; they
  do not drop genes.
- Default gates are `padj < 0.05` and `logFC > 1.0`. If `selected` is empty,
  lower `logfc_cutoff` (for example `0.25`) or check the contrast.
- Omit `sample_col` only when you have no biological replicates (cell-level DE).
- Per-gene `mechanism_class` is exploratory. For claims, pass `gene_sets=`
  and/or run `program_mechanism_permutation_calibrated`.

No `spliced`/`unspliced` (or `mature`/`nascent`) layers? You can still run
DE and enrichment — see
[Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html#de-only-no-nascent-layers).
To build those layers:
[Preparing spliced/unspliced data](https://scatrans.readthedocs.io/en/latest/tutorials/t_prepare_spliced_unspliced.html).

Next: [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html) ·
[Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html) ·
[FAQ](https://scatrans.readthedocs.io/en/latest/faq.html)

## Status

**0.10.x (Beta).** Import as `import scatrans as scat` and stick to names in
`scatrans.__all__`, `scat.pl`, and `scat.qc`. Details:
[API stability](https://scatrans.readthedocs.io/en/latest/api_stability.html).

## Citation

Please cite the preprint:

> Li, Z., James, A. W. & Li, S. scATrans: annotating single-cell differential
> expression as transcription- or stabilization-weighted using unspliced RNA.
> *bioRxiv* (2026). doi:10.64898/2026.08.03.740741

For the software itself, cite the Zenodo DOI above. Pin the installed
version in Methods (this tree is `scatrans==0.10.15`). See `CITATION.cff`.

## License

Software: [Apache License 2.0](LICENSE). Bundled GO/KEGG data may carry separate
terms—see the
[license page](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial redistribution.

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)
