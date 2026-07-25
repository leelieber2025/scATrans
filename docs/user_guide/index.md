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
3. **Interpret** — use soft per-gene labels for exploration. For stronger
   claims, pass `gene_sets=`, turn on `induction_matched=True`, and/or run
   `program_mechanism_permutation_calibrated`.
4. **Enrich** — on `result.selected` only (not on `mechanism_class` subsets).
5. **Plot and report** — `scat.pl.*`; see {doc}`../statistical_guidance` for
   what each column means.

## Defaults worth remembering

| Topic | Usual practice |
|-------|----------------|
| Entry point | `partition_de_by_mechanism` → `PartitionResult` |
| Gene list | DE membership (`result.selected`) |
| Mechanism claims | Prefer program tables over single-gene classes |
| Absolute placement | `program_mechanism_permutation_calibrated` with a frozen `de=` table |
| Detection | Optional `add_nascent_score=True` — separate from mechanism |
| Run metadata | `result.meta` and `adata.uns["scatrans"]` |

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
