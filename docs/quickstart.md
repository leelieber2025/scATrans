# Quickstart

Install the package, run one analysis, and read the main outputs. Links at the
end point to fuller tutorials and the user guide.

## 1. Install

```bash
pip install scatrans
# recommended when you have biological replicates:
pip install "scatrans[pseudobulk]"
```

```python
import scatrans as scat
print(scat.__version__)
```

### What data do you need?

| You have… | Start with |
|-----------|------------|
| AnnData with `spliced`/`unspliced` (or `mature`/`nascent`) | `partition_de_by_mechanism` below |
| Counts only (no nascent layers) | `differential_expression` — see [DE only](#de-only-no-nascent-layers) |
| Raw reads, no nascent layers yet | {doc}`tutorials/t_prepare_spliced_unspliced` — generate them with velocyto / kb-python / STARsolo / alevin-fry first |

If you plan to use PyDESeq2, Memento, or a full-gene enrichment background,
snapshot raw counts before HVG / normalize:

```python
scat.store_raw_counts(adata, layer="counts")
```

## 2. Pick the call

| Situation | Call |
|-----------|------|
| Default: DE + mechanism | `scat.partition_de_by_mechanism(...)` |
| Biological replicates | Set `sample_col="sample"` (or donor / library ID) |
| Curated gene programs | `gene_sets={...}`; use `induction_matched=True` if induction varies a lot |
| Absolute program placement | After partition: `program_mechanism_permutation_calibrated` |
| Optional detection score | `add_nascent_score=True` (does not change mechanism labels) |
| No velocity layers | `differential_expression` + enrichment |

DE builds the gene list. The residual only annotates mechanism on that list.

## 3. Run the default path

No AnnData handy? `scat.datasets.load_toy()` returns a small synthetic
object with `spliced`/`unspliced` layers already set — swap in your own
once the calls below run end to end.

```python
import scatrans as scat

adata = scat.datasets.load_toy()  # or your own AnnData

result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",       # "human" for human symbols
    sample_col="sample",    # toy has this; pass your replicate column when you have one
    # gene_sets=my_pathways,
    # induction_matched=True,
)
```

`groupby` / `target_group` / `reference_group` must match values in
`adata.obs`. `organism` selects the bundled gene-feature table and GO/KEGG
sets. The keyword default is `"mouse"`; pass `"human"` for human data.

### What to look at

```python
print(result.regime)            # capture OK? reliability in [0, 1]
print(len(result.selected))     # how many DE genes
print(result.selected.head())   # logFC, p_adj, mechanism columns
print(result.summary())         # padj_cutoff, logfc_cutoff, counts
```

| If this happens | What to do |
|-----------------|------------|
| `result.selected` is empty | Default is `padj < 0.05` and `logFC > 1.0`. Try `logfc_cutoff=0.25`, or print `result.gene_table[["logFC","p_adj"]].head()`. |
| Low gene-feature / enrichment mapping warning | Wrong `organism`, or symbol case / Ensembl vs symbol mismatch. |
| `sample_col` warning | You omitted replicates. Effect sizes are still usable; p-values are cell-level. |

Keyword defaults if you omit them: `organism="mouse"`, `logfc_cutoff=1.0`,
`padj_cutoff=0.05`, `sample_col=None`. Passing `sample_col` only switches
to pseudobulk + PyDESeq2 when there are at least 3 samples per group;
otherwise DE stays cell-level Wilcoxon. `load_toy()` ships 2 samples per
group, so the install check uses Wilcoxon. Full table: {doc}`api/index`.

Do not mix these with `active_score` / `filter_active_genes(preset="heuristic")`
defaults (`logfc_cutoff=0.35`, `de_method="t-test_overestim_var"`). Same
AnnData, different gene list.

| Field | Meaning |
|-------|---------|
| `result.regime` | Whether unspliced capture looks usable |
| `result.selected` | DE gene list plus soft mechanism labels |
| `result.gene_table` | All scored genes |
| `result.programs` | Only if you passed `gene_sets=` (relative to background) |
| `result.programs_induction_matched` | Only if `induction_matched=True` |
| `result.meta` | Version, DE source, thresholds, diagnostics |

Per-gene `mechanism_class` values are exploratory. Prefer program-level
tables when you report mechanism.

### Enrichment and plots

Enrich the DE list (`result.selected`), not genes split by `mechanism_class`:

```python
enrich = scat.run_enrichment(
    result.selected.index.tolist(),
    gene_sets="GO_Biological_Process",
    organism="mouse",
    adata=adata,
)
scat.pl.enrich_dotplot(enrich, top_n=15)
# scat.pl.volcano_plot(result.gene_table, top_n=10)
# scat.pl.comet_plot(result.gene_table, top_n=12)
```

### Absolute program placement (optional)

If you need placement against an empirical zero, not only “vs background”:

```python
de = result.gene_table[["logFC", "p_adj", "p_val"]]
cal = scat.program_mechanism_permutation_calibrated(
    adata,
    gene_sets=my_pathways,
    de=de,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    restrict_to_selected=True,
    n_perm=200,
)
print(cal[["observed_mean", "null_mean", "calibrated", "p_perm"]])
```

More detail: {doc}`user_guide/workflow`.

## 4. What next

| Goal | Page |
|------|------|
| Human LPS–PBMC example | {doc}`tutorials/t_gse226488_partition_mechanism` |
| Mouse example with real DE hits | {doc}`tutorials/t_ga_active_transcription` |
| Underpowered design (empty DE list) | {doc}`tutorials/t_ec_active_transcription` |
| DE + enrichment only | {doc}`tutorials/t_gse96583_standalone_de_enrichment` |
| Full workflow options | {doc}`user_guide/index` |
| Column meanings for reporting | {doc}`statistical_guidance` |
| Errors and API choice | {doc}`faq` |

---

## DE only (no nascent layers)

```python
adata, de_results = scat.differential_expression(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    # use_pseudobulk=True, sample_col="sample",  # still returns cell-level adata
)
candidates = scat.filter_active_genes(de_results, select_by="de")
```

`use_pseudobulk=True` aggregates internally; the returned `adata` stays
cell-level. Sample-level summary:
`adata.uns["scatrans"]["pseudobulk_obs"]`.

Enrichment and plotting work the same way:
{doc}`user_guide/standalone_de`.

## Lower-level residual table (optional)

Most users can skip this. Residual + DE without the partition wrapper:

```python
adata_res, significant, all_results = scat.active_score_simple(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
candidates = scat.filter_active_genes(all_results, select_by="de")
```

The second return value (`significant`) is a strict residual-and-DE filter and
is often empty. That is expected. Prefer `result.selected` from partition.
See {doc}`user_guide/workflow`.
