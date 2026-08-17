# FAQ / Troubleshooting

First analysis: {doc}`quickstart`. Reporting columns:
{doc}`statistical_guidance`. Assumptions: {doc}`domain_assumptions`.

## What should I call?

| Goal | Function |
|------|----------|
| Default analysis (DE + mechanism) | {func}`~scatrans.partition_de_by_mechanism` |
| Absolute program placement | {func}`~scatrans.program_mechanism_permutation_calibrated` |
| DE only, no nascent layers | {func}`~scatrans.differential_expression` — {doc}`user_guide/standalone_de` |
| Lower-level residual + DE table | {func}`~scatrans.active_score_simple` |

DE chooses the gene list. The residual labels mechanism on that list; it does
not drop DE hits or replace DE for discovery.

| Topic | Practical rule |
|-------|----------------|
| Gene list | DE only (`result.selected`) |
| Mechanism | Exploratory per-gene labels; prefer `gene_sets=` for program-level calls |
| Absolute placement | `program_mechanism_permutation_calibrated` (observed − label-shuffle null) |
| Detection | Optional `add_nascent_score=True` — does not set mechanism labels |
| Enrichment | On the DE list, not on `mechanism_class` splits |
| Replicates | Set `sample_col` when you have libraries or donors |

If capture looks thin or extreme, run
{func}`~scatrans.qc.regime_diagnosis`.

## Which program test should I use?

| Question | Tool |
|----------|------|
| Is the program shifted vs genome-wide background? | `gene_sets=` / `program_mechanism` |
| Still true after matching on induction (logFC)? | `induction_matched=True` / `program_mechanism_induction_matched` |
| Absolute place on the axis (empirical zero)? | `program_mechanism_permutation_calibrated` |

For absolute placement, report `calibrated` (and usually `null_mean`). Some
gene sets keep a structural offset under shuffled labels; that is what
calibration removes. Always pass a frozen `de=` table. See
{doc}`user_guide/workflow` and {doc}`method`.

## Do I need spliced/unspliced layers?

**For mechanism partition, yes** (`spliced`/`unspliced` or `mature`/`nascent`).
**For ordinary DE, no** — use {func}`~scatrans.differential_expression` on a
count matrix. Enrichment and plotting work the same either way. Tutorial:
{doc}`tutorials/t_gse96583_standalone_de_enrichment`.

Don't have those layers yet? {doc}`tutorials/t_prepare_spliced_unspliced`
covers generating them with velocyto / kb-python / STARsolo / alevin-fry
(or renaming labeling-based `new`/`old` counts) and merging them into an
existing AnnData.

## Can the residual replace DE?

**No.** Membership should come from DE. Older composite ranking
(`ranking_mode="composite"`, `run_default_pipeline(select_by="composite")`) is
not a recommended discovery path. Prefer
{func}`~scatrans.partition_de_by_mechanism` or pure DE.

## Why is `result.selected` empty?

Usually DE found no genes at **your cutoffs**, not a crash. The partition
defaults are `padj_cutoff=0.05` and `logfc_cutoff=1.0` (strict `>`). Check
`result.summary()` for the values that were used, then try a lower
`logfc_cutoff` (for example `0.25`) or inspect
`result.gene_table[["logFC", "p_adj"]]`.

The SCI endothelium tutorial is an underpowered design where an empty list is
the expected teaching point: {doc}`tutorials/t_ec_active_transcription`.

## I got a low mapping-rate warning

- **Gene features / mechanism:** `organism` defaults to `"mouse"`. Human
  symbols (`TP53`, `GAPDH`, …) need `organism="human"`.
- **External DE table:** the DE index must match `adata.var_names` (same ID
  type and case). A match rate below 20% now warns; 0% yields an empty
  `selected` list.
- **Gene sets / GO / GSEA:** symbols must match the library. Enrichr sets are
  usually **UPPERCASE**; try `gene_case="upper"` for mixed-case mouse symbols.

## `organism` — mouse or human?

The keyword default is **`"mouse"`** on partition, enrichment, and
`add_gene_features`. That is not detected from your data. Pass
`organism="human"` for human symbols, or enrichment and bias correction will
look empty / poorly mapped.

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

## `use_pseudobulk=True` wiped my `adata.obs`

Older builds returned the intermediate sample-level object, so

```python
adata, de = scat.differential_expression(adata, use_pseudobulk=True, sample_col="...")
```

left `obs` with only `sample_col`, `groupby`, `n_cells`, `total_counts`, and
`pb_x_source`. That is a bug. Current versions keep the **cell-level** AnnData
(same `obs` columns, embeddings, layers). Cells whose `groupby` value is
not the target or reference also stay (only `subset_col` drops cells).
The DE table is the second return value. Sample-level summary is in
`adata.uns["scatrans"]["pseudobulk_obs"]`. The same contract applies to
`active_score`, `*_simple`, `run_default_pipeline`, and
`partition_de_by_mechanism`.

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
column that matches `adata.var_names` exactly ({doc}`user_guide/gene_features`).

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
