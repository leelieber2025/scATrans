# Standalone Differential Expression

Use this path when you have no spliced/unspliced layers, or you only need
ordinary DE plus enrichment and plots.

```python
import scatrans as scat

# Before HVG / normalize / log1p
scat.store_raw_counts(adata, layer="counts")

adata, de_results = scat.differential_expression(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    # de_method="wilcoxon",
    # use_pseudobulk=True, sample_col="sample",
    # use_memento_de=True,
)

candidates = scat.filter_active_genes(de_results, select_by="de")
# or: padj_cutoff=0.05, logfc_cutoff=0.3
# downregulated: logfc_direction="down"
# both directions: logfc_direction="both"

enrich = scat.run_enrichment(
    candidates.index.tolist(),
    gene_sets="GO_Biological_Process",
    adata=adata,  # uses stored full-gene background when available
)
scat.pl.volcano_plot(de_results)
scat.pl.enrich_dotplot(enrich)
```

Backends match the rest of the package (scanpy methods, PyDESeq2 pseudobulk,
mixed models, optional Memento). Do **not** set `use_mixed_model` and
`use_memento_de` together. With MixedLM, reported `logFC` is
sample-mean-of-means log2FC — see {doc}`advanced`.

`use_pseudobulk=True` aggregates internally for DE. The returned AnnData stays
cell-level (same `obs` columns, embeddings, and layers). Sample-level summary
is in `adata.uns["scatrans"]["pseudobulk_obs"]`. The same contract applies to
`active_score`, `active_score_simple`, `differential_expression_simple`,
`run_default_pipeline`, and `partition_de_by_mechanism`.

The result table plugs into `filter_active_genes`, enrichment helpers, and
`scat.pl.*`. Example script: `examples/memento_de_example.py`.

## Raw counts

Count-based backends (PyDESeq2, Memento) need raw integers. This pattern
breaks them:

```python
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# .X is log-transformed HVGs only — not usable as counts
```

Do this instead:

```python
scat.store_raw_counts(adata, layer="counts")  # first

sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata, de_res = scat.differential_expression(adata, use_memento_de=True, ...)

# DE / enrichment on the full pre-HVG gene set:
adata_full = scat.restore_raw_counts(adata, full_genes=True)
```

`store_raw_counts(adata, mode="auto")` is safe to call more than once and can
recover counts from `adata.raw` when `.X` is already normalized.

**After `anndata.concat()`:** concatenation drops `.uns`, including the
`log1p` marker. Re-set `combined.uns["log1p"] = {"base": None}` or pass
`de_preprocess="none"` if you know `.X` is already log-normalized.
