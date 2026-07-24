# scATrans Documentation

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![PyPI downloads](https://img.shields.io/pepy/dt/scatrans.svg)](https://pepy.tech/project/scatrans)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/scatrans.svg)](https://anaconda.org/bioconda/scatrans)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/leelieber2025/scATrans/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21365873.svg)](https://doi.org/10.5281/zenodo.21365873)

## What scATrans does

After you run differential expression, some genes go up because of
transcription and some because mRNA is stabilized. **scATrans** takes a
DE-selected gene list and soft-labels that difference using a
reference-corrected nascent (unspliced) residual.

One call does the usual path:

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,  # spliced/unspliced or mature/nascent layers
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",  # or "human"
    de="builtin",
    # sample_col="sample",   # preferred with biological replicates
    # gene_sets=my_pathways, # optional program-level table
)
result.regime           # capture-quality check
result.selected.head()  # DE genes + mechanism labels
result.summary()
```

| Role | Job |
|------|-----|
| **DE** | Builds the gene list |
| **Mechanism** | Transcription vs. stabilization labels on that list |
| **Detection** (optional) | Extra nascent scores via `add_nascent_score=True` — does not drive labels |

Program-level summaries (`gene_sets=`) are more trustworthy than single-gene
calls. No nascent layers? You can still run ordinary DE, enrichment, and
plots. Scope and caveats: {doc}`faq`.

## Get started

```bash
pip install scatrans
```

| Step | Page |
|------|------|
| 1. Install (and optional extras) | {doc}`installation` |
| 2. Run the first example | {doc}`quickstart` |
| 3. Walk a real notebook | {doc}`tutorials/index` |
| 4. Dig into options | {doc}`user_guide/index` |

:::{note}
**Beta (0.10.x).** Use `import scatrans as scat` and public names in
`scatrans.__all__`, `scat.pl`, and `scat.qc`. Private leaf modules may move
before 1.0 — see {doc}`api_stability`.
:::

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

pip, Bioconda, optional extras.
:::

:::{grid-item-card} Quickstart {octicon}`rocket;1em;`
:link: quickstart
:link-type: doc

First run, what to read back, enrichment.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

Worked notebooks on real and synthetic data.
:::

:::{grid-item-card} User Guide {octicon}`book;1em;`
:link: user_guide/index
:link-type: doc

Workflow, DE backends, enrichment, plots.
:::

:::{grid-item-card} Method {octicon}`beaker;1em;`
:link: method
:link-type: doc

Residual math and optional permutation.
:::

:::{grid-item-card} Statistical Guidance {octicon}`alert;1em;`
:link: statistical_guidance
:link-type: doc

What each column means for reporting.
:::

:::{grid-item-card} API Reference {octicon}`code;1em;`
:link: api/index
:link-type: doc

Functions and parameters.
:::

:::{grid-item-card} FAQ {octicon}`question;1em;`
:link: faq
:link-type: doc

Scope, data requirements, common errors.
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/leelieber2025/scATrans

Source, issues, contributions.
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
