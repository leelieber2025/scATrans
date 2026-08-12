# User Guide

Read {doc}`../quickstart` first if you have not run scATrans yet. This guide
covers the main options after that.

## Which page do I need?

| Task | Page |
|------|------|
| Default DE → mechanism workflow | {doc}`workflow` |
| DE without spliced/unspliced layers | {doc}`standalone_de` |
| ORA, GO, KEGG, GSEA | {doc}`enrichment` |
| Volcano, comet, enrichment plots | {doc}`plotting` |
| Residual permutation, mixed models, adaptive score | {doc}`advanced` |
| Length/intron tables or a custom GTF | {doc}`gene_features` |

## Default analysis steps

1. **Prepare** — AnnData with `spliced`/`unspliced` (or `mature`/`nascent`).
   Set `sample_col` if you have replicates. Optionally run
   `scat.store_raw_counts(adata, layer="counts")`.
2. **Partition** — `scat.partition_de_by_mechanism(...)`. Check
   `result.regime` and `result.selected`.
3. **Interpret** — treat per-gene labels as exploratory. For claims, pass
   `gene_sets=`, set `induction_matched=True`, and/or run
   `program_mechanism_permutation_calibrated`.
4. **Enrich** — on `result.selected` only (not on `mechanism_class` subsets).
5. **Plot and report** — `scat.pl.*`; see {doc}`../statistical_guidance` for
   what each column means.

## Defaults worth remembering

These are the installed keyword defaults, not recommendations:

| Call | Default |
|------|---------|
| `partition_de_by_mechanism` | `organism="mouse"`, `logfc_cutoff=1.0`, `padj_cutoff=0.05`, `sample_col=None` |
| `run_default_pipeline` | `select_by="composite"` (deprecated), `run_go_enrichment=True` |
| `filter_active_genes` | `select_by="composite"`, `preset=None` (permissive) |

Usual practice: use `partition_de_by_mechanism`, pass `sample_col` when you
have replicates, treat `result.selected` as the gene list, and report
mechanism at the program level. Full table: {doc}`../api/index`.

Stuck? {doc}`../faq`.

```{toctree}
:maxdepth: 2

workflow
standalone_de
enrichment
plotting
advanced
gene_features
```
