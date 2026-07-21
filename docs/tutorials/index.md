# Tutorials

Pre-rendered notebooks on real and synthetic data. Local reproduction notes are
at the end of this page.

| Notebook | Focus |
|----------|--------|
| LPS-PBMC (`GSE226488`) | Primary demonstration of `partition_de_by_mechanism`, program-level inference, optional `add_nascent_score` |
| SCI vs UN (EC) | Same workflow on a low-power 3-vs-3 design (empty DE list by design) |
| GA vs Ctrl | Powered mouse design with pseudobulk DE and mechanism partition |
| Standalone DE | DE and enrichment without nascent layers |
| Visualization gallery | `scat.pl` on synthetic tables |
| Gene UpSet | Multi-method gene overlap |

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Partition DE by mechanism — LPS-PBMC (10x v3)
:link: t_gse226488_partition_mechanism
:link-type: doc
:class-card: scat-meta

Reliability pre-flight, DE selection, per-gene mechanism annotation, and
program-level transcription-versus-stabilization calls (NF-κB vs ARE/ZFP36),
with an induction-matched check and pluggable DE front end.
+++
`GSE226488` subset · nascent layers · program-level mechanism
:::

:::{grid-item-card} Partition by mechanism — SCI versus UN (low power)
:link: t_ec_active_transcription
:link-type: doc
:class-card: scat-meta

`partition_de_by_mechanism` on mouse spinal cord endothelium (3 versus 3).
Capture is adequate, but DE selects no genes at this power—an illustration of
underpowered designs.
+++
`EC.h5ad` · nascent layers · pseudobulk · power illustration
:::

:::{grid-item-card} Partition by mechanism — GA versus Ctrl
:link: t_ga_active_transcription
:link-type: doc
:class-card: scat-meta

Same workflow on a better-powered mouse design (three individuals per group):
pseudobulk DE selects a program; scATrans partitions with per-gene annotation,
threshold sensitivity, and enrichment.
+++
`GA_test.h5ad` · nascent layers · pseudobulk
:::

:::{grid-item-card} Differential expression and enrichment (no nascent layers)
:link: t_ec_standalone_de_enrichment
:link-type: doc
:class-card: scat-meta

`differential_expression` across Wilcoxon, PyDESeq2, and Memento backends, then
ORA, KEGG, GO, GSEA, redundancy reduction, and plotting.
+++
`EC.h5ad` · no layers required · ORA / GSEA
:::

:::{grid-item-card} Visualization gallery (synthetic)
:link: t_synthetic_visualization
:link-type: doc
:class-card: scat-meta

Plotting helpers on synthetic results: volcano, comet, enrichment panels, and
export utilities.
+++
Synthetic · plotting only
:::

:::{grid-item-card} Gene UpSet (optional)
:link: t_ec_gene_upset
:link-type: doc
:class-card: scat-meta

Gene-set membership UpSet plots for comparing gene lists across DE methods.
+++
Gene membership visualization
:::

::::

```{toctree}
:maxdepth: 1
:hidden:

t_gse226488_partition_mechanism
t_ec_active_transcription
t_ga_active_transcription
t_ec_standalone_de_enrichment
t_synthetic_visualization
t_ec_gene_upset
```

## Reproducing notebooks locally

Install documentation extras and open notebooks under `docs/tutorials/`. Large
example AnnData files are not always shipped with the package; each notebook
documents data paths and download notes.
