# FAQ / Troubleshooting

Also see {doc}`domain_assumptions` and {doc}`statistical_guidance`.

## What should I call?

| Goal | Function |
|------|----------|
| Default analysis (DE + mechanism) | {func}`~scatrans.partition_de_by_mechanism` |
| DE only, no nascent layers | {func}`~scatrans.differential_expression` — {doc}`user_guide/standalone_de` |
| Lower-level residual + DE table | {func}`~scatrans.active_score_simple` |

One idea to keep straight: **DE chooses the gene list; the residual labels
mechanism on that list.** The residual does not drop DE hits and is not meant
to replace DE for discovery.

| Topic | Practical rule |
|-------|----------------|
| Gene list | DE only (`result.selected`, or `filter_active_genes(..., select_by="de")`) |
| Mechanism | Soft labels on selected genes; prefer `gene_sets=` for program-level calls |
| Detection | Optional `add_nascent_score=True` — never drives mechanism labels |
| Enrichment | Run ORA on the DE list, not on `mechanism_class` splits |

Per-gene mechanism labels are modest; pathway pooling is usually more useful.
Signal tracks intron capture quality — check {func}`~scatrans.qc.regime_diagnosis`
when capture looks thin or 3′-biased.

## Do I need spliced/unspliced layers?

**For mechanism partition, yes** (`spliced`/`unspliced` or `mature`/`nascent`).
**For ordinary DE, no** — use {func}`~scatrans.differential_expression` on a
count matrix. Enrichment and plotting work the same either way. Tutorial:
{doc}`tutorials/t_ec_standalone_de_enrichment`.

## Can the residual replace DE?

**No.** Membership should come from DE. Older composite ranking
(`ranking_mode="composite"`, `run_default_pipeline(select_by="composite")`) is
not a recommended discovery path. Prefer
{func}`~scatrans.partition_de_by_mechanism` or pure DE.

## Why is `result.selected` empty?

Usually DE found no genes at your cutoffs (power, thresholds, or a quiet
contrast) — not a failed package. Check design and DE settings first. The SCI
endothelium tutorial shows an underpowered case where an empty list is the
expected teaching point: {doc}`tutorials/t_ec_active_transcription`.

## `padj_cutoff` or `pval_cutoff`?

Use **`padj_cutoff=`**. Legacy `pval_cutoff=` still works but, when adjusted
p-values exist, it filters **`p_adj`**, not raw `p_val`. The modern name
avoids warnings and makes papers clearer.

## GSEA is empty or warns about mapping / `gene_case`

- Rank genes with a **signed** metric (prefer `logFC`).
- Symbol case must match the gene set (Enrichr libraries are usually
  **UPPERCASE**). For mixed-case mouse symbols, try `gene_case="upper"`.
- Mapping rate below ~20% warns; 0% returns empty with
  `reason="no_ranked_genes_mapped"`.

See {doc}`user_guide/enrichment`.

## `ValueError` with `use_mixed_model=True`

1. **Sample size:** need ≥4 biological samples per group and ≥6 random-effect
   groups total. With 3 vs 3, use pseudobulk + PyDESeq2 instead.
2. **Clash:** do not combine `use_mixed_model=True` with `use_memento_de=True`.

See {doc}`user_guide/advanced`.

## `ImportError` for `pydeseq2` / `scvelo` / `gseapy` / `memento`

Those are optional extras:

```bash
pip install "scatrans[pseudobulk]"        # PyDESeq2
pip install "scatrans[advanced]"          # scVelo
pip install "scatrans[gene_features]"     # custom GTF tables
pip install "scatrans[memento]"           # Memento
pip install "scatrans[gsea]"              # GSEA (gseapy)
```

See {doc}`installation`.

## PyDESeq2 / Memento complain about counts

They need **raw integer counts**. If you already ran HVG + normalize + log1p,
`.X` is no longer usable. Before preprocessing:

```python
scat.store_raw_counts(adata, layer="counts")
# or: store_raw_counts(adata, mode="auto")  # also try adata.raw
```

See {doc}`user_guide/standalone_de`.

## Warnings after `anndata.concat()`

`ad.concat()` drops `.uns` by default, including the `log1p` marker. Either
restore it (`combined.uns["log1p"] = {"base": None}`) or set
`de_preprocess="none"` when you know `.X` is already log-normalized.

## Huge `logFC` values (e.g. >20)

Often a scanpy `rank_genes_groups` artifact when the reference group is near
zero for that gene — not specific to scATrans. Check raw spliced/unspliced
counts (e.g. `scat.pl.velocity_phase_portraits`) before reporting. Example:
{doc}`tutorials/t_ec_active_transcription`.

## Global unspliced fraction > 50%

Often technical (ambient RNA, swapped layers, nuclear enrichment / gDNA, or a
very immature population) rather than clean nascent signal:

```python
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
```

Reliability is high in a normal band (~10–45% unspliced) and lower at both
extremes. Partition always runs this check and scales mechanism confidence.
High reliability means the proxy is not obviously broken — not that residual
beats DE.

## Missing length / intron features (all NaN for some genes)

`add_gene_features` matches `adata.var_names`. Genes absent from the feature
table get `NaN` and skip bias correction. Custom tables need a `gene_name`
column that matches exactly ({doc}`user_guide/gene_features`).

## Design / sample-size warnings

When `sample_col` or pseudobulk is set, design notes land under
`adata.uns["scatrans"]["diagnostics"]["design"]` (and in the log). Partition
uses the same path. Open that block if notebook logs are easy to miss.

## Bundled KEGG for commercial use?

Not under Apache-2.0 alone. KEGG needs a separate commercial license from
Kanehisa Laboratories for non-academic use. To skip bundled files, pass an
Enrichr library explicitly, e.g. `run_kegg(..., kegg_library="KEGG_2021")`.
See {doc}`license`.

---

## Lower-level `active_score` only

You only need this section if you call `active_score` /
`active_score_simple` yourself. The primary path returns DE genes in
`result.selected` and does not use the built-in `significant` mask for
membership. Thresholds: {doc}`statistical_guidance`.

### Why is the built-in `significant` list empty?

It requires DE gates **and** residual gates, including residual FDR from
`use_permutation=True`. With permutation off (the default), the list is empty
on purpose. For a DE gene list from `all_results`, use
`filter_active_genes(..., select_by="de")`.

### Why does residual ranking surface genes with `p_adj ≈ 1`?

The residual is independent of DE significance. Genes filtered or filled to
`padj=1` by a DE backend can still show positive unspliced excess. That is not
a DE gene list — use `select_by="de"` or partition for membership.

---

Still stuck? Open an issue on
[GitHub](https://github.com/leelieber2025/scATrans/issues).
