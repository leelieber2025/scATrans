# Tutorials

Worked, end-to-end examples on real and synthetic data. Every notebook ships
with its outputs already rendered, so you can read it here without running
anything — or reproduce it locally with the commands at the bottom of this
page.

:::{important}
Active-transcription tutorials use spliced/unspliced layers. That scoring path
is still **experimental** and under validation. For production DE + enrichment
+ plotting without velocity layers, start with the **standalone** notebook.
:::

The two **active transcription** notebooks share a deliberately chosen theme: a
real 3-vs-3 replicate design where statistical power — not the tool —
determines how much signal survives. The **standalone** notebook covers the
common case of an ordinary count matrix with no velocity layers, and the
**visualization** gallery is a pure plotting tour on synthetic data.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} 🧬 Active transcription — real SCI vs. UN
:link: t_ec_active_transcription
:link-type: doc
:class-card: scat-meta

The full `active_score` workflow on endothelial cells from mouse spinal cord
(**with** spliced/unspliced layers): heuristic, pseudobulk, permutation, gamma
robustness, bias correction, and advanced mode — reported honestly on a
low-power 3-vs-3 design.
+++
`EC.h5ad` · velocity-aware · PyDESeq2 / scVelo
:::

:::{grid-item-card} 🔬 Active transcription — higher-power GA vs. Ctrl
:link: t_ga_active_transcription
:link-type: doc
:class-card: scat-meta

The same `active_score` pipeline on a better-powered mouse dataset (~330 cells
per replicate). A side-by-side companion to the EC notebook that shows what the
output looks like when cells-per-replicate — i.e. power — is no longer the
bottleneck.
+++
`GA_test.h5ad` · velocity-aware · pseudobulk
:::

:::{grid-item-card} 📊 Differential expression + enrichment (no velocity)
:link: t_ec_standalone_de_enrichment
:link-type: doc
:class-card: scat-meta

For the majority of users with a plain count matrix: `differential_expression`
across three backends (Wilcoxon, PyDESeq2, Memento), then the full enrichment
toolkit — ORA, KEGG, GO (all ontologies), GSEA, redundancy reduction — and the
plotting gallery.
+++
`EC.h5ad` · no layers required · ORA / GSEA
:::

:::{grid-item-card} 🔗 Gene overlap across DE methods (UpSet)
:link: t_ec_gene_upset
:link-type: doc
:class-card: scat-meta

Which genes are called by *several* DE backends, and what do the robust ones do?
The gene-level UpSet trio — `build_gene_membership` → `gene_upsetplot` →
`common_genes` — stacks Wilcoxon / *t*-test / pseudobulk into one figure, then
feeds the common-up / common-down genes into enrichment. Works equally well on
your own external DE tables (Seurat / DESeq2 / CSV).
+++
`EC.h5ad` · multi-method DE · bring-your-own tables · UpSet / ORA
:::

:::{grid-item-card} 🎨 Visualization gallery
:link: t_synthetic_visualization
:link-type: doc
:class-card: scat-meta

A scene-by-scene tour of nearly every `scat.pl.*` function — what each figure
answers, the key keyword arguments, parameter variants, multi-panel `ax=`
composition, and export helpers. Synthetic data only, fully offline.
+++
synthetic · plotting API · offline
:::
::::

```{admonition} Dataset provenance
:class: note

The EC notebooks use endothelial cells from mouse spinal cord, comparing
uninjured controls (**UN**, 3 replicates) against spinal cord injury (**SCI**,
3 replicates) — a subset of Squair et al. (2021), chosen precisely because that
paper is *about* the false-discovery risk of single-cell pseudoreplication. See
{doc}`../references` for the full citation and GEO accession (GSE165003).
`GA_test.h5ad` is a bundled mouse test dataset (GA vs. Ctrl, three individuals
per group) from Li et al. (2026); raw data GEO
[GSE266598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266598) —
see {doc}`../references`. The visualization gallery uses **synthetic** data
only.
```

```{toctree}
:hidden:
:maxdepth: 1

t_ec_active_transcription
t_ga_active_transcription
t_ec_standalone_de_enrichment
t_ec_gene_upset
t_synthetic_visualization
```

## Reproducing these notebooks locally

Notebooks ship with their outputs already rendered. To re-run them yourself:

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev,advanced,pseudobulk,gene_features,memento,gsea]"

# Real-data active-transcription tutorials (need the bundled h5ad at repo root):
jupyter lab docs/tutorials/t_ec_active_transcription.ipynb
jupyter lab docs/tutorials/t_ga_active_transcription.ipynb

# Standalone DE + enrichment (EC.h5ad, no velocity layers used):
jupyter lab docs/tutorials/t_ec_standalone_de_enrichment.ipynb

# Gene-level UpSet across DE methods (EC.h5ad):
jupyter lab docs/tutorials/t_ec_gene_upset.ipynb

# Plotting gallery (synthetic only — no h5ad required):
jupyter lab docs/tutorials/t_synthetic_visualization.ipynb
```

`EC.h5ad` and `GA_test.h5ad` are included at the repository root, which is where
the notebooks load them from (`sc.read_h5ad("../../EC.h5ad")`). The
visualization tutorial generates its own AnnData in memory.
