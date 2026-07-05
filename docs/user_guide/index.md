# User Guide

`active_score` performs differential expression, reference-gamma unspliced
excess calculation, optional bias correction, composite scoring, and stores
results plus diagnostics. Downstream steps commonly include gene filtering
with `filter_active_genes`, functional enrichment, and plotting.

The internal `significant` list uses strict thresholds. The complete results
table is returned as `all_results`; use `filter_active_genes` for custom
criteria. Diagnostics are available under
`adata_res.uns["scatrans"]["diagnostics"]`.

```{toctree}
:maxdepth: 2

workflow
enrichment
plotting
advanced
gene_features
standalone_de
```
