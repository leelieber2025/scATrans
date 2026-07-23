# scATrans

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Documentation Status](https://readthedocs.org/projects/scatrans/badge/?version=latest)](https://scatrans.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

**scATrans** is a Python package for mechanism-aware single-cell differential
analysis. A standard differential expression (DE) step defines which genes
changed; scATrans partitions those genes into *transcription-driven* versus
*stabilization-driven* classes using the nascent (unspliced) RNA residual—a
distinction that total-count fold change alone cannot resolve.

| Component | Role |
|-----------|------|
| DE | Gene-list membership |
| Mechanism annotation | Residual-based transcription vs. stabilization labels |
| Detection (optional) | `add_nascent_score=True` adds active-transcription scores; does not drive mechanism labels |

The primary workflow requires spliced and unspliced layers (or mature and
nascent layers, e.g. from kb-python) and is most informative at the pathway or
program level. The package also supports conventional DE without nascent layers
(scanpy, PyDESeq2 pseudobulk, linear mixed models, optional Memento), enrichment
(ORA, GSEA, GO, KEGG), and plotting.

Full documentation: [Read the Docs](https://scatrans.readthedocs.io/en/latest/).

## Requirements

- Python 3.10+
- AnnData object with a condition column in `.obs`
- For mechanism analysis: `spliced`/`unspliced` or `mature`/`nascent` layers

## Installation

```bash
pip install scatrans
# or: conda install -c conda-forge -c bioconda scatrans
```

Optional extras (scVelo, gene-feature CLI, PyDESeq2, Memento, GSEA) and
development installs:
[installation guide](https://scatrans.readthedocs.io/en/latest/installation.html).

## Quickstart

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,  # AnnData with spliced/unspliced or mature/nascent layers
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    de="builtin",  # method name, kwargs dict, DataFrame, or callable
    # sample_col="sample",  # preferred when biological replicates exist
    # add_nascent_score=True,  # optional detection columns
    gene_sets=my_pathways,  # optional program-level table
    # induction_matched=True,  # induction-controlled program tests
)
result.regime    # reliability pre-flight (global unspliced fraction)
result.selected  # DE-selected genes with soft mechanism annotation
result.programs  # program-level table when gene_sets is provided
result.summary() # program-first overview
```

Further reading:

- [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html)
- [Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html)
- [User Guide](https://scatrans.readthedocs.io/en/latest/user_guide/index.html)
- [FAQ](https://scatrans.readthedocs.io/en/latest/faq.html) (scope and limitations)
- [Statistical Guidance](https://scatrans.readthedocs.io/en/latest/statistical_guidance.html)

## Status

scATrans is **0.10.x (Beta)**. Prefer `import scatrans as scat` and names in
`scatrans.__all__`, `scat.pl`, and `scat.qc`. See
[API stability](https://scatrans.readthedocs.io/en/latest/api_stability.html).

## Citation

If you use scATrans in published work, cite the software via the Zenodo DOI
above and the manuscript when available. See `CITATION.cff`.

## License

Software: [Apache License 2.0](LICENSE). Bundled gene-set data (GO, KEGG) may
carry separate terms; see the
[license page](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial redistribution.

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)
