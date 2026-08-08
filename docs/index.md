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

Differential expression tells you **which** genes change. When you also have
nascent or unspliced layers, scATrans helps you ask **how** those changes look
on a transcription-versus-stabilization axis.

- Soft labels on individual genes are exploratory.
- Stronger claims usually come from **gene programs**, not single genes.
- DE still defines the gene list. scATrans does not replace DE for discovery.

No nascent layers? You can still run DE, enrichment, and plotting.

```text
 spliced + unspliced AnnData
             │
             ▼
    ┌──────────────────┐
    │ 1. DE             │   selects the gene list (Wilcoxon / PyDESeq2 / …)
    │    (defines it)   │
    └──────────────────┘
             │  DE-selected genes only
             ▼
    ┌──────────────────┐
    │ 2. Mechanism      │   unspliced-excess residual →
    │    (annotates it) │   transcription- vs stabilization-driven (soft, per gene)
    └──────────────────┘
             │  optional: gene_sets={program: [genes]}
             ▼
    ┌──────────────────┐
    │ 3. Program table  │   pooled, induction-matched / permutation-calibrated
    │    (claims here)  │   → the level to report, not single genes
    └──────────────────┘
```

Step 1 decides *which* genes; steps 2–3 decide *how* — never the reverse.

## Where to go

| Goal | Page |
|------|------|
| Install and run a first analysis | {doc}`installation` → {doc}`quickstart` |
| Pick an API or program test | {doc}`faq` and the table below |
| Full workflow, enrichment, plots | {doc}`user_guide/index` · {doc}`tutorials/index` |
| What each column means for a paper | {doc}`statistical_guidance` |
| Math and formal API | {doc}`method` · {doc}`api/index` |

### A sensible path

1. Install: `pip install "scatrans[pseudobulk]"` (add `[gsea]` if you need GSEA).
2. Follow {doc}`quickstart`.
3. Work through a real notebook
   ({doc}`tutorials/t_gse226488_partition_mechanism` or
   {doc}`tutorials/t_ga_active_transcription`).
4. Enrich the DE gene list; keep mechanism claims at the program level when you can.

### Default call

No data yet? `scat.datasets.load_toy()` returns a synthetic AnnData with
`spliced`/`unspliced` layers already in place, so the block below runs as-is,
no download needed. Swap in your own AnnData once it works; if it lacks
nascent layers, see {doc}`tutorials/t_prepare_spliced_unspliced`.

```python
import scatrans as scat

adata = scat.datasets.load_toy()  # or your own AnnData

result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",  # or "human"
    de="builtin",
    sample_col="sample",  # set when you have biological replicates
    # gene_sets=my_pathways,
    # induction_matched=True,
)
print(result.regime)           # capture quality
print(result.selected.head())  # DE genes + mechanism labels
print(result.summary())
```

### Which tool for which job

| Goal | Use |
|------|-----|
| DE + mechanism (usual case) | `partition_de_by_mechanism` |
| DE / enrichment only (no velocity layers) | `differential_expression` + `run_enrichment` |
| Program vs genome-wide background | `gene_sets=` → `result.programs` |
| Program with induction controlled | `induction_matched=True` |
| Absolute program placement | `program_mechanism_permutation_calibrated` |
| Optional nascent detection score | `add_nascent_score=True` (does not set mechanism labels) |
| Enrichment | On `result.selected`, not on `mechanism_class` splits |

Practical habits: define membership with DE, treat per-gene classes as soft,
prefer programs for claims, and pass `sample_col` when you have replicates.
See {doc}`faq`.

:::{note}
scATrans is **0.10.x (Beta)**. Prefer `import scatrans as scat` and the public
names in `scatrans.__all__`, `scat.pl`, and `scat.qc`. Details:
{doc}`api_stability`.
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

First run after install.
:::

:::{grid-item-card} Tutorials {octicon}`play;1em;`
:link: tutorials/index
:link-type: doc

Worked notebooks (read online or re-run).
:::

:::{grid-item-card} User Guide {octicon}`book;1em;`
:link: user_guide/index
:link-type: doc

Workflow, enrichment, plots, advanced options.
:::

:::{grid-item-card} FAQ {octicon}`question;1em;`
:link: faq
:link-type: doc

API choice and common errors.
:::

:::{grid-item-card} Statistical Guidance {octicon}`alert;1em;`
:link: statistical_guidance
:link-type: doc

What each column is for.
:::

:::{grid-item-card} Method {octicon}`beaker;1em;`
:link: method
:link-type: doc

Residual and calibration math.
:::

:::{grid-item-card} API Reference {octicon}`code;1em;`
:link: api/index
:link-type: doc

Functions and parameters.
:::

:::{grid-item-card} GitHub {octicon}`mark-github;1em;`
:link: https://github.com/leelieber2025/scATrans

Source and issues.
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
