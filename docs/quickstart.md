# Quickstart

## Primary workflow

{func}`~scatrans.partition_de_by_mechanism` runs DE selection, residual-based
mechanism annotation, and a reliability pre-flight from the global unspliced
fraction. Optional program-level inference requires `gene_sets=`.

| Need | Entry point |
|------|-------------|
| Mechanism partition (recommended) | {func}`~scatrans.partition_de_by_mechanism` |
| Residual / scoring engine only | {func}`~scatrans.active_score_simple` / {func}`~scatrans.active_score` |
| DE without nascent layers | {func}`~scatrans.differential_expression` ({doc}`user_guide/standalone_de`) |

Scope and limitations: {doc}`faq`.

```python
import scatrans as scat

result = scat.partition_de_by_mechanism(
    adata,  # spliced/unspliced or mature/nascent layers required
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    de="builtin",  # method name, kwargs dict, DataFrame, or callable
    # add_nascent_score=True,  # optional detection columns (not for mechanism)
    gene_sets=my_pathways,  # optional program-level table
)
print(result.regime)
print(result.selected.head())
print(result.programs)

# Residual and scoring engine only (lower-level)
adata_res, significant, all_results = scat.active_score_simple(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
```

Gene features are attached when missing. With sufficient biological replicates,
`active_score_simple` prefers pseudobulk + PyDESeq2; otherwise it uses single-cell
Wilcoxon DE. Advanced options (permutation, mixed models, Memento) are documented
in the {doc}`user_guide/index`.

## Lower-level path (`active_score`)

The example below runs scoring, filtering, enrichment, and plotting as separate
steps. The primary workflow above composes these steps and adds the mechanism
partition.

```python
import scanpy as sc
import scatrans as scat

# 1. Load data
adata = sc.read_h5ad("your_data.h5ad")

# 2. Snapshot raw counts and velocity layers before HVG/normalization.
#    sidecar=True (default) stores a full-gene snapshot in .uns that survives
#    later subsetting.
scat.store_raw_counts(adata, layer="counts")

# 3. Standard preprocessing (adjust to the analysis)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 4. Gene features for optional bias correction
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# 5. Score (residual + composite active_score)
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    show_plot=False,
)

print(all_results.head())

# 6. Filter the full table (built-in significant is often empty by design)
candidates = scat.filter_active_genes(
    all_results,
    preset="heuristic",  # or "pseudobulk" / "permissive"
    # padj_cutoff=0.05,
    # logfc_cutoff=0.3,
)

print(f"Filtered genes: {len(candidates)}")

# 7. Enrichment
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",
    organism="mouse",
    adata=adata,  # uses stored raw gene universe as background when available
    padj_cutoff=0.05,
)
print(enrich_res.head())

kegg_res = scat.run_kegg(
    gene_list=candidates.index.tolist(),
    organism="mouse",
    adata=adata,
)

# 8. Plots
scat.pl.enrich_dotplot(enrich_res, top_n=15, title="GO Enrichment")
scat.pl.enrich_dotplot(kegg_res, top_n=10, title="KEGG Pathways")
# scat.pl.comet_plot(all_results, top_n=12)
# scat.pl.volcano_plot(all_results, top_n=10)
```

Further options: {doc}`user_guide/enrichment`. For DE without spliced/unspliced
layers, replace step 5 with `scat.differential_expression(...)`
({doc}`user_guide/standalone_de`).

## Raw counts and layer snapshots

Call `store_raw_counts` after load and QC, before HVG or normalization:

- writes the current `.X` to `layers["counts"]`;
- with `sidecar=True` (default), also stores a label-indexed snapshot of full
  obs × var counts and velocity layers under
  `adata.uns['scatrans']['raw_snapshot']`.

Layers follow AnnData axis alignment (trimmed by HVG or cell subsetting). The
**snapshot** lives in `.uns` and is re-aligned by cell/gene name after
subsetting, `copy()`, and `write_h5ad()`.

```python
# Current gene set (respects subsetting / reordering)
adata_raw = scat.restore_raw_counts(adata, inplace=False)

# Full pre-HVG gene universe (counts + velocity layers)
adata_full = scat.restore_raw_counts(adata, full_genes=True)
```

Pass `adata=` to `run_enrichment` / `run_kegg` to use the stored gene list as the
background universe. For large objects:

```python
scat.store_raw_counts(adata, sidecar="ondisk", snapshot_path="raw_snapshot.h5ad")
```

```{note}
`save_raw=True` is deprecated: `adata.raw` is often reserved for log-normalized
data, and the sidecar snapshot already preserves full-gene raw counts. Prefer
`restore_raw_counts(..., full_genes=True)`.
```
