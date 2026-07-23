# Tutorials

Pre-rendered notebooks on real and synthetic data. HTML builds show stored
outputs; you do not need the large `.h5ad` files to read the results online.

## Start here

| If you want… | Open |
|--------------|------|
| The main DE → mechanism story (human LPS-PBMC) | {doc}`t_gse226488_partition_mechanism` |
| The same workflow where DE actually selects genes (mouse) | {doc}`t_ga_active_transcription` |
| What underpowered DE looks like (empty gene list by design) | {doc}`t_ec_active_transcription` |
| DE + enrichment with no spliced/unspliced layers | {doc}`t_ec_standalone_de_enrichment` |
| Plotting helpers only (no real `.h5ad`) | {doc}`t_synthetic_visualization` |
| Overlap of gene lists across DE backends | {doc}`t_ec_gene_upset` |

New to scATrans? Read {doc}`../quickstart`, then the LPS-PBMC notebook. The
SCI (EC) and GA notebooks use the **same** entry point
(`partition_de_by_mechanism`) on purpose: they show how **design power**, not
API choice, decides whether there is anything to partition.

### Three partition notebooks (same API, different designs)

| Notebook | Design | What you should expect |
|----------|--------|------------------------|
| LPS-PBMC (`GSE226488`) | Human 10x PBMC, resting vs LPS 4 h | Full path: regime check, DE selection, per-gene labels, program-level calls, pluggable DE, optional detection columns |
| GA vs Ctrl | Mouse, 3 individuals per group | Pseudobulk DE selects a real program; partition + threshold sensitivity + enrichment |
| SCI vs UN (EC) | Mouse endothelium, 3 vs 3 | Capture is fine (`regime="ok"`), but DE selects **no** genes at this power — that empty list is the teaching point, not a software failure |

### Runtime and dependencies (local runs)

Times are rough wall-clock on a laptop once the data file is on disk. Online
docs already include figures and tables.

| Notebook | Approx. time | Extras | Data file (repo root unless noted) |
|----------|--------------|--------|------------------------------------|
| LPS-PBMC | 5–15 min | `[pseudobulk]` recommended | `GSE226488_PBMC_tutorial_subset.h5ad` (~4.4k cells; not shipped) |
| SCI / EC partition | 2–5 min | `[pseudobulk]` | `EC.h5ad` (not always shipped) |
| GA partition | 5–15 min | `[pseudobulk]` | `GA_test.h5ad` (not always shipped) |
| Standalone DE + enrichment | 10–30 min | `[pseudobulk,memento,gsea]` for all cells | `EC.h5ad` |
| Visualization gallery | 1–3 min | none | In-memory synthetic tables |
| Gene UpSet | 5–15 min | `[pseudobulk]` optional | `EC.h5ad` |

Install extras, for example:

```bash
pip install "scatrans[pseudobulk,memento,gsea]"
```

Dataset provenance and access notes: {doc}`../references`. Input layers and
raw-count handling: {doc}`../user_guide/workflow` (section *Input data and
layers*).

---

## Notebook cards

| Notebook | Focus |
|----------|--------|
| LPS-PBMC (`GSE226488`) | Primary demo of `partition_de_by_mechanism`, program-level inference, optional `add_nascent_score` |
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
export utilities. Useful when you only want figure recipes.
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

```bash
pip install -e ".[dev]"   # or: pip install "scatrans[docs]" if available
# plus analysis extras as needed, e.g. pseudobulk / gsea
jupyter lab docs/tutorials/
```

Large example AnnData objects are **not** always included in the git tree or
PyPI sdist (size). Place each file at the path the notebook uses (usually the
repository root, loaded as `../../....h5ad` from `docs/tutorials/`).

| File | Used by | Notes |
|------|---------|--------|
| `GSE226488_PBMC_tutorial_subset.h5ad` | LPS-PBMC | Downsampled subset (~4.4k cells). Build from GEO if you do not have a local copy — see {doc}`../references` and the notebook **Reproduce** section. |
| `EC.h5ad` | SCI partition, standalone DE, gene UpSet | Endothelium subset of GSE165003; see {doc}`../references`. |
| `GA_test.h5ad` | GA partition | Mouse GA vs Ctrl with velocity layers; see {doc}`../references`. |

If a data file is missing, you can still follow the pre-rendered HTML on
[Read the Docs](https://scatrans.readthedocs.io/). The visualization gallery
runs without external data.
