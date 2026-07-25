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

## Install

```bash
pip install scatrans
# or: conda install -c conda-forge -c bioconda scatrans
```

Optional extras (PyDESeq2, GSEA, scVelo, …):
[installation guide](https://scatrans.readthedocs.io/en/latest/installation.html).

## First run

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,  # needs spliced/unspliced or mature/nascent layers
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",  # or "human"
    de="builtin",
    # sample_col="sample",   # set this when you have biological replicates
    # gene_sets=my_pathways, # optional pathway / program table
    # induction_matched=True,
)
print(result.regime)           # data-quality check on unspliced capture
print(result.selected.head())  # DE genes + soft mechanism labels
print(result.summary())
# Absolute program placement (optional):
# scat.program_mechanism_permutation_calibrated(adata, gene_sets, de=frozen_de, ...)
```

That is the recommended entry point. Next:

1. [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html)
2. [Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html)
3. [FAQ](https://scatrans.readthedocs.io/en/latest/faq.html) if something looks off

## Status

**0.10.x (Beta).** Import as `import scatrans as scat` and stick to names in
`scatrans.__all__`, `scat.pl`, and `scat.qc`. Details:
[API stability](https://scatrans.readthedocs.io/en/latest/api_stability.html).

## Citation

Cite the Zenodo DOI above (and the manuscript when available). For analyses
tied to package version **0.10.9**, use `scatrans==0.10.9`. See `CITATION.cff`.

## License

Software: [Apache License 2.0](LICENSE). Bundled GO/KEGG data may carry separate
terms—see the
[license page](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial redistribution.

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)
