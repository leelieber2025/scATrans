# scATrans Documentation

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/leelieber2025/scATrans/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

## Overview

**scATrans** is a Python package for mechanism-aware analysis of single-cell
differential expression. Given spliced and unspliced layers (or mature and
nascent layers), it takes a DE-selected gene list and partitions genes into
*transcription-driven* versus *stabilization-driven* classes using a
reference-corrected nascent residual
({func}`~scatrans.partition_de_by_mechanism`).

| Role | Responsibility |
|------|----------------|
| **DE** | Defines gene-list membership |
| **Mechanism** | Annotates transcription vs. stabilization from the residual |
| **Detection** (optional) | `add_nascent_score=True` adds active-transcription scores; does not drive mechanism labels |

Program-level inference (`gene_sets=`) is preferred over single-gene mechanism
claims. Scope, limitations, and reporting conventions: {doc}`faq`,
{doc}`statistical_guidance`.

The package also provides conventional DE without nascent layers, enrichment
(ORA, GSEA, GO, KEGG with bundled sets), and plotting utilities.

## Installation

```bash
pip install scatrans
# or: conda install -c conda-forge -c bioconda scatrans
```

Optional extras and editable installs: {doc}`installation`.

## Minimal example

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,  # AnnData with spliced/unspliced or mature/nascent layers
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",  # or "human"
    de="builtin",  # method name, kwargs dict, DataFrame, or callable
    # add_nascent_score=True,  # optional detection columns
    gene_sets=my_pathways,  # optional program-level table
)
result.regime           # reliability pre-flight
result.selected.head()  # DE-selected genes with mechanism annotation
result.programs         # present when gene_sets is provided
```

**Reading order:** {doc}`installation` → {doc}`quickstart` →
{doc}`tutorials/index`.

:::{note}
**Beta (0.10.x).** Use the public import surface `import scatrans as scat` and
names in `scatrans.__all__`, `scat.pl`, and `scat.qc`. Leaf modules such as
`scatrans.tl.active` may move before 1.0; see {doc}`api_stability`.
:::

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

pip, Bioconda, optional extras, and logging.
:::

:::{grid-item-card} Quickstart {octicon}`rocket;1em;`
:link: quickstart
:link-type: doc

Primary workflow, lower-level scoring, enrichment, and count snapshots.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

Worked notebooks on real and synthetic data.
:::

:::{grid-item-card} User Guide {octicon}`book;1em;`
:link: user_guide/index
:link-type: doc

Workflow, backends, enrichment, plotting, and advanced options.
:::

:::{grid-item-card} Method {octicon}`beaker;1em;`
:link: method
:link-type: doc

Unspliced excess, bias correction, composite score, and permutation.
:::

:::{grid-item-card} Statistical Guidance {octicon}`alert;1em;`
:link: statistical_guidance
:link-type: doc

Output columns, safe uses, and reporting checklist.
:::

:::{grid-item-card} API Reference {octicon}`code;1em;`
:link: api/index
:link-type: doc

Public functions, parameters, and autosummary.
:::

:::{grid-item-card} API Stability {octicon}`shield-check;1em;`
:link: api_stability
:link-type: doc

Stable imports versus implementation detail (pre-1.0).
:::

:::{grid-item-card} FAQ {octicon}`question;1em;`
:link: faq
:link-type: doc

Scope, limitations, and common errors.
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/leelieber2025/scATrans

Source code, issues, and contributions.
:::
::::

```{toctree}
:hidden: true
:maxdepth: 3
:titlesonly: true

installation
quickstart
tutorials/index
user_guide/index
method
statistical_guidance
domain_assumptions
api_stability
api/index
ROADMAP
references
faq
changelog
license
GitHub <https://github.com/leelieber2025/scATrans>
```
