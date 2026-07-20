# Tutorials

Worked, end-to-end examples on real and synthetic data. Every notebook ships
with its outputs already rendered, so you can read it here without running
anything — or reproduce it locally with the commands at the bottom of this
page.

:::{important}
**Start here:** the **LPS-PBMC** notebook demonstrates the recommended primary
workflow — DE selects the changed genes, scATrans partitions them into
transcription- vs stabilization-driven, with a decisive program-level call. The
**SCI** and **GA** notebooks run the same `partition_de_by_mechanism` workflow on
velocity datasets at two power levels. For a plain count matrix with no velocity
layers, use the **standalone** DE + enrichment notebook.
:::

All three velocity notebooks use the same recommended `partition_de_by_mechanism`
workflow. The **LPS-PBMC** notebook is the lead example (10x v3, a curated
program-level mechanism call); the **SCI** and **GA** notebooks are a matched pair
where statistical power — not the tool — determines how much signal survives. The
**standalone** notebook covers the common no-velocity case, and the
**visualization** gallery is a pure plotting tour on synthetic data.

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} 🧭 Partition DE by mechanism — LPS-PBMC (10x v3)
:link: t_gse226488_partition_mechanism
:link-type: doc
:class-card: scat-meta

**The recommended workflow.** `partition_de_by_mechanism` on PBMC resting vs
LPS-4h: reliability pre-flight → DE selects → per-gene soft mechanism labels →
decisive **program-level** transcription-vs-stabilization call, validated on the
textbook NF-κB-transcription vs ARE-decay (TTP/ZFP36) biology, with an
induction-matched honesty check and a pluggable DE front-end.
+++
`GSE226488` subset · velocity-aware · program-level mechanism
:::

:::{grid-item-card} 🧬 Partition by mechanism — SCI vs UN (low power)
:link: t_ec_active_transcription
:link-type: doc
:class-card: scat-meta

`partition_de_by_mechanism` on endothelial cells from mouse spinal cord, a hard
3-vs-3 design. The honest lesson: capture regime is fine, but DE confidently
selects **nothing** at this power — so there is nothing to partition, and the
pipeline says so rather than manufacture hits.
+++
`EC.h5ad` · velocity-aware · pseudobulk · power lesson
:::

:::{grid-item-card} 🔬 Partition by mechanism — GA vs Ctrl (powered)
:link: t_ga_active_transcription
:link-type: doc
:class-card: scat-meta

The same workflow on a better-powered mouse dataset (3 individuals per group):
here pseudobulk DE selects a real program, and scATrans partitions it into
transcription-driven vs stabilization-driven with per-gene annotation, threshold
sensitivity, a mechanism plot, and enrichment.
+++
`GA_test.h5ad` · velocity-aware · pseudobulk · real partition
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
see {doc}`../references`. The mechanism-partition notebook uses human PBMC,
resting vs LPS-4h (10x 3′ v3.1), from Derbois et al. (2023), GEO
[GSE226488](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226488) — a
**downsampled subset** regenerated locally, not shipped in the repo. The
visualization gallery uses **synthetic** data only.
```

```{toctree}
:hidden:
:maxdepth: 1

t_gse226488_partition_mechanism
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

# Mechanism-partition tutorial (recommended). The GSE226488 subset is NOT shipped:
# regenerate it from GEO GSE226488 (STARsolo with spliced/unspliced -> assemble ->
# downsample) and place GSE226488_PBMC_tutorial_subset.h5ad at the repo root. The
# exact commands are in the notebook's "Reproduce" section.
jupyter lab docs/tutorials/t_gse226488_partition_mechanism.ipynb

# Velocity-aware partition tutorials, SCI (low power) and GA (powered):
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
