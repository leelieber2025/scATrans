# Quickstart

## Before you start

| You have… | Use |
|-----------|-----|
| AnnData with `spliced`/`unspliced` (or `mature`/`nascent`) layers | {func}`~scatrans.partition_de_by_mechanism` — start here |
| Counts only (no nascent layers) | {func}`~scatrans.differential_expression` — {doc}`user_guide/standalone_de` |

Python 3.9+, `pip install scatrans`. For replicate-aware DE backends, also
install `pip install "scatrans[pseudobulk]"`.

## One function for most users

{func}`~scatrans.partition_de_by_mechanism`:

1. Runs DE and keeps changed genes  
2. Annotates each selected gene (transcription- vs. stabilization-driven)  
3. Checks nascent-layer quality (`result.regime`)  
4. Optionally pools genes into programs if you pass `gene_sets=`  

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",  # or "human"
    de="builtin",
    sample_col="sample",  # set when biological replicates exist
    # gene_sets=my_pathways,
)
```

### What to look at

```python
print(result.regime)
# e.g. {"regime": "ok", "reliability": 0.9, "message": ...}

print(result.selected.head())
# DE genes with soft mechanism columns (transcription_support, mechanism_class, …)

print(result.summary())
# short program-first overview when gene_sets= was supplied
```

| Field | Meaning |
|-------|---------|
| `result.selected` | Your gene list (DE membership) plus mechanism labels |
| `result.gene_table` | Full gene table from the run |
| `result.programs` | Program-level table if you passed `gene_sets=` |
| `result.regime` | Whether unspliced capture looks usable |
| `result.meta` | Diagnostics and run settings |

Gene features (length / intron) are filled in when missing. With enough
replicates, the builtin DE path prefers pseudobulk + PyDESeq2; otherwise it
falls back to single-cell Wilcoxon. More options: {doc}`user_guide/index`.

### Enrichment on the selected list

Enrich the **DE list**, not mechanism-class subsets:

```python
genes = result.selected.index.tolist()
enrich = scat.run_enrichment(
    genes,
    gene_sets="GO_Biological_Process",
    organism="mouse",
    adata=adata,
)
scat.pl.enrich_dotplot(enrich, top_n=15)
# scat.pl.comet_plot(result.gene_table, top_n=12)
# scat.pl.volcano_plot(result.gene_table, top_n=10)
```

## Good habits (two minutes)

**Store raw counts early** if you will use PyDESeq2, Memento, or full-gene
enrichment backgrounds:

```python
scat.store_raw_counts(adata, layer="counts")  # before HVG / normalize / log1p
```

**Check capture quality** (also run automatically inside partition):

```python
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
```

**Prefer pathway calls** over hard single-gene mechanism claims when
induction is strong or capture is modest. See {doc}`faq`.

## Worked examples

| Goal | Go to |
|------|--------|
| Full human LPS-PBMC story | {doc}`tutorials/t_gse226488_partition_mechanism` |
| Mouse design with real DE hits | {doc}`tutorials/t_ga_active_transcription` |
| Underpowered design (empty list by design) | {doc}`tutorials/t_ec_active_transcription` |
| DE only, no nascent layers | {doc}`tutorials/t_ec_standalone_de_enrichment` |

## DE without nascent layers

```python
adata, de_results = scat.differential_expression(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
)
candidates = scat.filter_active_genes(de_results, select_by="de")
```

Same enrichment and plotting helpers apply. Details:
{doc}`user_guide/standalone_de`.

## Raw-count snapshots (optional but useful)

Call `store_raw_counts` after load/QC and **before** HVG or normalization:

- writes current `.X` to `layers["counts"]`
- with `sidecar=True` (default), keeps a full-gene snapshot under
  `adata.uns['scatrans']['raw_snapshot']` that survives subsetting and h5ad I/O

```python
adata_raw = scat.restore_raw_counts(adata, inplace=False)
adata_full = scat.restore_raw_counts(adata, full_genes=True)  # pre-HVG universe
```

Large objects:

```python
scat.store_raw_counts(adata, sidecar="ondisk", snapshot_path="raw_snapshot.h5ad")
```

```{note}
`save_raw=True` is deprecated. Prefer the sidecar snapshot and
`restore_raw_counts(..., full_genes=True)`.
```

## Lower-level scoring (optional)

Most users never need this. If you want the residual + DE table without the
partition wrapper:

```python
adata_res, significant, all_results = scat.active_score_simple(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
# Build a DE list yourself:
candidates = scat.filter_active_genes(all_results, select_by="de")
```

The second return value (`significant`) is a strict residual+DE conjunction
and is often empty — that is expected. Prefer `select_by="de"` or the primary
partition API. Details: {doc}`user_guide/workflow`, {doc}`faq`.
