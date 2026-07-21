# User Guide

## Organization

| Page | Content |
|------|---------|
| {doc}`workflow` | Primary entry point, filtering, design diagnostics, layers |
| {doc}`standalone_de` | DE without nascent layers |
| {doc}`enrichment` | ORA, GSEA, GO, KEGG, export |
| {doc}`plotting` | `scat.pl` helpers and figure export |
| {doc}`advanced` | Permutation, mixed models, adaptive score, mechanism helpers |
| {doc}`gene_features` | Length/intron tables and GTF CLI |

## Conventions

- **Primary API:** {func}`~scatrans.partition_de_by_mechanism` →
  {class}`~scatrans.PartitionResult`
- **Lower-level tables:** `all_results` from `active_score` /
  `differential_expression`; filter with `filter_active_genes`
- **Metadata:** `result.meta` and `adata.uns["scatrans"]`
- **Scope and reporting:** {doc}`../faq`, {doc}`../statistical_guidance`

```{toctree}
:maxdepth: 2

workflow
standalone_de
enrichment
plotting
advanced
gene_features
```
