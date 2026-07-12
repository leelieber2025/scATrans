# scATrans Documentation

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](https://github.com/leelieber2025/scATrans/blob/main/LICENSE)

## Single-cell active transcription analysis

**scATrans** is a Python toolkit for single-cell differential analysis. It is
primarily designed for datasets that contain spliced/unspliced (or
mature/nascent) RNA layers: in this setting it computes a composite **active
transcription score** that integrates differential expression with
reference-based excess unspliced RNA to rank genes.

It also supports conventional differential expression workflows (no velocity
data required) via scanpy, PyDESeq2 pseudobulk, linear mixed models, or
optional Memento. Functional enrichment (ORA, GSEA, GO, KEGG) uses bundled
gene sets with consistent universe handling, and a set of visualization
functions is provided.

## Try it now

```bash
pip install scatrans
```

```python
import scatrans as scat

result = scat.run_default_pipeline(
    adata,                       # AnnData with spliced/unspliced (or mature/nascent) layers
    groupby="condition", target_group="Disease", reference_group="Control",
    sample_col="sample",          # optional; auto-selects pseudobulk when >=3 replicates/group
    organism="mouse",             # or "human"
)
result["candidates"].head()       # ranked, filtered genes
result["enrichment"].head()       # GO enrichment on those genes
```

New here? Follow {doc}`installation` → {doc}`quickstart` → {doc}`tutorials/index`
(real data, fully worked) in that order.

:::{note}
**API stability.** scATrans is currently **Beta (0.10.x)**. Import the package
as `import scatrans as scat` and rely on the names in `scatrans.__all__`,
`scat.pl`, and `scat.qc`. Leaf modules such as `scatrans.tl.active` are an
implementation detail and may move between minor releases. The complete
contract is documented in {doc}`api_stability`.
:::

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} Installation {octicon}`plug;1em;`
:link: installation
:link-type: doc

Install scATrans with pip, with optional extras for pseudobulk, velocity, or
Memento backends.
:::

:::{grid-item-card} Quickstart {octicon}`rocket;1em;`
:link: quickstart
:link-type: doc

A minimal end-to-end example: load data, score genes, filter, enrich, plot.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

Worked examples on real spinal-cord-injury endothelial cell data, with and
without spliced/unspliced layers.
:::

:::{grid-item-card} User Guide {octicon}`book;1em;`
:link: user_guide/index
:link-type: doc

Core workflow, DE backends, enrichment, plotting, and advanced options.
:::

:::{grid-item-card} Statistical Guidance {octicon}`alert;1em;`
:link: statistical_guidance
:link-type: doc

What each output column means, and what it should (and should not) be used
for in a paper or supplement.
:::

:::{grid-item-card} API Reference {octicon}`code;1em;`
:link: api/index
:link-type: doc

Detailed description of every public function in scATrans.
:::

:::{grid-item-card} References {octicon}`mortar-board;1em;`
:link: references
:link-type: doc

Tutorial data source, and the methods/libraries scATrans builds on.
:::

:::{grid-item-card} FAQ / Troubleshooting {octicon}`question;1em;`
:link: faq
:link-type: doc

Common errors and how to fix them, in one place.
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/leelieber2025/scATrans

Found a bug? Want to contribute? Check out the source and open an issue.
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
