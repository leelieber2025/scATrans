# Tutorials

HTML on Read the Docs is pre-executed (tables and figures already there).
Re-run locally only if you have the `.h5ad` files at the repo root.

## Pick a notebook

| If you want… | Open |
|--------------|------|
| **No `spliced`/`unspliced` layers yet** — velocyto, kb-python, STARsolo, alevin-fry, or labeling data | {doc}`t_prepare_spliced_unspliced` |
| Full DE → mechanism story (human LPS-PBMC) | {doc}`t_gse226488_partition_mechanism` |
| Same API with real DE hits (mouse, 3 vs 3) | {doc}`t_ga_active_transcription` |
| Underpowered design (empty DE list on purpose) | {doc}`t_ec_active_transcription` |
| DE + enrichment, no nascent layers | {doc}`t_gse96583_standalone_de_enrichment` |
| Plot gallery only | {doc}`t_synthetic_visualization` |
| Gene overlap across DE backends | {doc}`t_ec_gene_upset` |

**If you are new:** read {doc}`../quickstart`, then either
{doc}`t_gse226488_partition_mechanism` or {doc}`t_ga_active_transcription`.
Open the SCI/EC notebook only if you want the empty-list lesson. If your
`AnnData` does not have nascent layers at all, start with
{doc}`t_prepare_spliced_unspliced` instead.

SCI (EC) and GA share the same entry point (`partition_de_by_mechanism`).
Whether genes are selected depends on the design, not on the function name.

### Partition notebooks

| Notebook | Design | What to expect |
|----------|--------|----------------|
| LPS-PBMC | Human 10x, resting vs LPS 4 h | Full path: DE, labels, programs, absolute placement |
| GA vs Ctrl | Mouse, 3 individuals per group | Real DE hits, programs, enrichment |
| SCI vs UN (EC) | Mouse endothelium, 3 vs 3 | Capture OK but **0 DE genes**; sample- vs cell-level contrast |

## Run locally

```bash
pip install "scatrans[pseudobulk,gsea]"
# or from a clone: pip install -e ".[dev,pseudobulk,gsea]"
jupyter lab docs/tutorials/
```

Put data files at the **repository root**. Notebooks load them with a relative
path, for example `sc.read_h5ad("../../EC.h5ad")`.

| File | Used by |
|------|---------|
| `GSE226488_PBMC_tutorial_subset.h5ad` | LPS-PBMC |
| `GA_test.h5ad` | GA |
| `EC.h5ad` | SCI partition, gene UpSet |
| `kang_ifnb_tutorial_subset.h5ad` | Standalone DE + enrichment |

These files are large and are not on PyPI. Citations: {doc}`../references`.
Without them you can still read the HTML on Read the Docs. The visualization
gallery does not need an external file.

Rough runtime with data on disk: LPS / GA about 10–30 min (calibration is
slower); EC / UpSet about 2–15 min; standalone DE + enrichment about 5–15 min;
synthetic plots about 1–3 min.

---

## Notebook cards

::::{grid} 1 1 2 2
:gutter: 3

:::{grid-item-card} Preparing spliced/unspliced data
:link: t_prepare_spliced_unspliced
:link-type: doc

velocyto / kb-python / STARsolo / alevin-fry → AnnData layers, plus merging
and sanity checks.
+++
Start here if you have no nascent layers yet
:::

:::{grid-item-card} Partition — LPS-PBMC (GSE226488)
:link: t_gse226488_partition_mechanism
:link-type: doc

Full human story: DE, mechanism labels, programs, absolute placement.
+++
Nascent layers · best first real-data tutorial
:::

:::{grid-item-card} Partition — GA vs Ctrl
:link: t_ga_active_transcription
:link-type: doc

Powered mouse design with real DE hits, GO programs, enrichment.
+++
`GA_test.h5ad` · pseudobulk
:::

:::{grid-item-card} Partition — SCI vs UN (low power)
:link: t_ec_active_transcription
:link-type: doc

Same API; DE finds nothing — that is the point. Sample- vs cell-level table.
+++
`EC.h5ad` · empty DE lists
:::

:::{grid-item-card} DE + enrichment (no nascent layers)
:link: t_gse96583_standalone_de_enrichment
:link-type: doc

Wilcoxon / PyDESeq2 / Memento, ORA, GO, KEGG, GSEA, plots.
+++
`kang_ifnb_tutorial_subset.h5ad` · counts only
:::

:::{grid-item-card} Visualization gallery
:link: t_synthetic_visualization
:link-type: doc

`scat.pl` helpers on synthetic tables, including palettes and `cmap=`.
+++
No external data
:::

:::{grid-item-card} Gene UpSet across DE methods
:link: t_ec_gene_upset
:link-type: doc

Overlap of gene lists from different DE backends.
+++
`EC.h5ad`
:::
::::

```{toctree}
:hidden: true
:maxdepth: 1

t_prepare_spliced_unspliced
t_gse226488_partition_mechanism
t_ga_active_transcription
t_ec_active_transcription
t_gse96583_standalone_de_enrichment
t_synthetic_visualization
t_ec_gene_upset
```
