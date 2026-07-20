# scATrans

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Documentation Status](https://readthedocs.org/projects/scatrans/badge/?version=latest)](https://scatrans.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

scATrans is a single-cell differential-analysis toolkit with an RNA-velocity
twist: standard DE **selects** the changed genes, then scATrans **partitions**
them by mechanism — *transcription-driven* vs *stabilization-driven* — from the
nascent (unspliced) signal, a call a fold-change alone cannot make.

It also supports conventional differential expression workflows (no velocity
data required) using scanpy, PyDESeq2 pseudobulk, linear mixed models, or
optional Memento. Functional enrichment (ORA, GSEA, GO, KEGG) uses bundled gene
sets with consistent universe handling, and a set of visualization functions is
provided.

**📚 Documentation, tutorials, and API reference:
[Read the Docs](https://scatrans.readthedocs.io/en/latest/).**

## Installation

```bash
pip install scatrans
# or: conda install -c conda-forge -c bioconda scatrans
```

Optional extras (scVelo, gene features, PyDESeq2, Memento, GSEA) and a source /
editable dev setup are covered in the
[Installation guide](https://scatrans.readthedocs.io/en/latest/installation.html).

## Quickstart

```python
import scatrans as scat

# DE selects the changed genes; scATrans partitions them by MECHANISM.
result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition", target_group="Disease", reference_group="Control",
    organism="mouse",
    de="builtin",            # or a DE method name / precomputed DE table / callable
    gene_sets=my_pathways,   # optional -> program-level mechanism table
)
result.regime      # reliability pre-flight
result.selected    # DE-selected genes + per-gene mechanism annotation
result.programs    # decisive program-level transcription-vs-stabilization calls
```

See the [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html),
[Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html), and
[User Guide](https://scatrans.readthedocs.io/en/latest/user_guide/index.html) for
the full workflow, DE backends, enrichment, plotting, and reporting guidance
([Statistical Guidance](https://scatrans.readthedocs.io/en/latest/statistical_guidance.html)).

## License

Software (Python source) is licensed under [Apache License 2.0](LICENSE).
Bundled gene-set data (GO, KEGG) carries its own licensing terms — see
[License](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial use.

## Author

**Zhao Li (李钊)**  
Email: [leelieber@gmail.com](mailto:leelieber@gmail.com)
