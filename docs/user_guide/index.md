# User Guide

Start with {doc}`../quickstart` if you have not run the package yet. This
guide covers the main knobs after that first call.

## Pages

| Page | When you need it |
|------|------------------|
| {doc}`workflow` | Primary API, filtering, design checks, layers |
| {doc}`standalone_de` | DE without spliced/unspliced layers |
| {doc}`enrichment` | ORA, GSEA, GO, KEGG |
| {doc}`plotting` | `scat.pl` helpers and figure export |
| {doc}`advanced` | Permutation, mixed models, adaptive score |
| {doc}`gene_features` | Length/intron tables and GTF CLI |

## Conventions

- **Default entry point:** {func}`~scatrans.partition_de_by_mechanism` →
  {class}`~scatrans.PartitionResult`
- **Gene list:** DE membership (`result.selected` or
  `filter_active_genes(..., select_by="de")`)
- **Run metadata:** `result.meta` and `adata.uns["scatrans"]`
- **Scope / reporting:** {doc}`../faq`, {doc}`../statistical_guidance`

```{toctree}
:maxdepth: 2

workflow
standalone_de
enrichment
plotting
advanced
gene_features
```
