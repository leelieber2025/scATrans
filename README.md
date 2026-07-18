# scATrans

[!\[PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[!\[Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[!\[Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[!\[Documentation Status](https://readthedocs.org/projects/scatrans/badge/?version=latest)](https://scatrans.readthedocs.io/en/latest/?badge=latest)
[!\[CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[!\[License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[!\[DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

scATrans is a Python toolkit for single-cell differential analysis. It is
primarily designed for datasets that contain spliced/unspliced (or
mature/nascent) RNA layers. In this setting it computes a composite active
transcription score that integrates differential expression with
reference-based excess unspliced RNA to rank genes.

It also supports conventional differential expression workflows (no
velocity data required) using scanpy, PyDESeq2 pseudobulk, linear mixed
models, or optional Memento. Functional enrichment (ORA, GSEA, GO, KEGG)
uses bundled gene sets with consistent universe handling, and a set of
visualization functions is provided.

**📚 Full documentation, tutorials, and the complete API reference are on**
[**Read the Docs**](https://scatrans.readthedocs.io/en/latest/)**.**

## Installation

```bash
# From PyPI
pip install scatrans

# Or from Bioconda
conda install -c conda-forge -c bioconda scatrans

# Optional extras (PyPI): advanced (scVelo) mode, gene features, pseudobulk DE (PyDESeq2), Memento, GSEA
pip install "scatrans\[advanced,gene\_features,pseudobulk,memento,gsea]"
```

See [Installation](https://scatrans.readthedocs.io/en/latest/installation.html)
for extras, source installs, and logging setup.

## Quickstart

```python
import scatrans as scat

# One-liner pipeline: score → filter → GO enrichment
result = scat.run\_default\_pipeline(
    adata,
    groupby="condition",
    target\_group="Disease",
    reference\_group="Control",
    sample\_col="sample",   # optional; auto-selects pseudobulk when >=3 replicates/group
    organism="mouse",
)
print(result\["candidates"].head())
print(result\["enrichment"].head())
```

See the [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html)
for a complete end-to-end walkthrough, the
[Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html)
for fully worked, real-data notebooks (with and without RNA-velocity
layers), and the
[User Guide](https://scatrans.readthedocs.io/en/latest/user_guide/index.html)
for DE backends, enrichment, plotting, and advanced options.

## Before reporting results in a paper

`active\_score` is a **composite heuristic rank**, not a p-value or FDR on
its own. See
[Statistical Guidance](https://scatrans.readthedocs.io/en/latest/statistical_guidance.html)
for what each output column means, safe vs. unsafe uses, and a reporting
checklist before you cite scATrans results in a manuscript or supplement.
Domain conventions (upregulation-oriented scoring, residual vs DE, cutoff
names, GSEA ranks, within-run λ scale) are spelled out in
[Domain Assumptions](docs/domain_assumptions.md)
(also on Read the Docs after the next docs deploy).

## License

Software (Python source) is licensed under [Apache License 2.0](LICENSE).
Bundled gene-set data (GO, KEGG) carries its own licensing terms — see
[License](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial use.

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)

