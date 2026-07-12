# FAQ / Troubleshooting

Common questions and errors, collected in one place. Each links to the
fuller explanation. Domain conventions (why top-N is not DE, why residual
can score without significant `p_adj`, cutoff names, etc.) are summarized in
{doc}`domain_assumptions`.

## "My `significant` gene list is empty" — is that a bug?

No. The built-in `significant` mask uses strict, all-of thresholds on
logFC, `p_adj`, `unspliced_excess_residual`, `active_score`, and
`unspliced_excess_fdr` (see {doc}`statistical_guidance`). On modestly
powered real designs it is often empty by design. Use the full
`all_results` table with `filter_active_genes(preset="heuristic")` or your
own cutoffs (see {doc}`user_guide/workflow`) — ranking + biological
follow-up is the intended workflow, not a `p<0.05` gate on every run.

## Top-N by `active_score` includes genes with `p_adj ≈ 1` — is that a bug?

No. The residual (nascent excess) leg of `active_score` is **independent of
DE significance**. Genes that DE backends did not test or filled as neutral
(e.g. PyDESeq2 independent filtering → `padj` set to 1) can still score from
positive unspliced excess. Top-N ranking is therefore **not** a
DE-significant gene list. For DE-gated candidates use
`filter_active_genes(..., padj_cutoff=0.05, logfc_cutoff=...)` (preferred;
legacy `pval_cutoff=` still works but maps to **adjusted** p) or the
built-in `significant` conjunction. See {doc}`statistical_guidance` for the
full default cutoff table.

## Can I compare `active_score` across cell types or datasets?

**Not as absolute numbers.** Soft-scale λ is estimated from the genes present
in **that** run (`median(positive)/ln(2)`). The same gene with the same raw
logFC can score ~40 in one subset and ~70 in another if the background gene
pool (or HVG filter) differs. Safe uses: rank genes **within one**
`active_score` call. For cross-group claims use `logFC` / `p_adj` / residual
on a shared analysis design, or re-run with a **shared gene universe** knowing
scores still are not absolute units. Inspect the within-run λ scale via:

```python
print(adata.uns["scatrans"]["diagnostics"]["scoring"])
# typically: lambda_fc, lambda_res, lambda_pval, …
```

See {doc}`domain_assumptions` and {doc}`statistical_guidance`.

## Should I use `padj_cutoff` or `pval_cutoff`?

**Prefer `padj_cutoff=`** on `filter_active_genes`, `extract_gene_lists`, and
enrichment helpers. Legacy `pval_cutoff=` still works but is a misleading
name: when adjusted p-values exist it filters **`p_adj` / `p.adjust`**, not
raw `p_val`. Using the modern name avoids deprecation warnings and makes
reporting clearer.

## `run_gsea` returns empty / warns about mapping rate or `gene_case`

Preranked GSEA needs **signed** ranks (prefer `logFC`; `active_score` is
not auto-selected). It also checks symbol overlap with gene sets (same
`_check_gene_set_mapping_rate` gate as ORA): below **20%** mapping rate you
get a warning with input vs gene-set examples; at **0%** the result is empty
with `reason="no_ranked_genes_mapped"`. Duplicate gene IDs (e.g. after
case-folding) keep the entry with **max |score|**. Enrichr libraries are
typically **UPPERCASE** — for mixed-case mouse symbols (`Tp53`) pass
`gene_case="upper"`. See {doc}`user_guide/enrichment`.

## `ValueError` from `use_mixed_model=True`

Two common causes:

1. **Sample-size gate:** the mixed-model path requires **≥4 biological
   samples per group** and **≥6 total random-effect groups**. With fewer
   replicates (e.g. 3 vs. 3), use `use_pseudobulk=True` +
   `pseudobulk_de_backend="pydeseq2"` instead.
2. **Backend clash:** `use_mixed_model=True` and `use_memento_de=True` are
   mutually exclusive (both are cell-level DE backends). Choose one.

See {doc}`user_guide/advanced`.

## `ImportError` for `pydeseq2` / `scvelo` / `gseapy` / `memento`

These are optional extras, not installed by the base `pip install scatrans`:

```bash
pip install "scatrans[pseudobulk]"        # PyDESeq2
pip install "scatrans[advanced]"          # scVelo (mode="advanced")
pip install "scatrans[gene_features]"     # gtfparse (custom gene-feature tables)
pip install "scatrans[memento]"           # Memento (use_memento_de=True)
pip install "scatrans[gsea]"              # GSEA (run_gsea), pulls in gseapy
```

See {doc}`installation`.

## `differential_expression(..., use_memento_de=True)` or PyDESeq2 raises a data-shape / non-integer-count error

Count-based backends need **raw integer counts**. A common mistake is
running HVG selection + `normalize_total` + `log1p` first, which leaves
`.X` log-transformed. Call `scat.store_raw_counts(adata, layer="counts")`
(or `ensure_raw_counts`) **before** any preprocessing. See
{doc}`user_guide/standalone_de`.

## After `anndata.concat()`, I get double-log1p / preprocessing warnings

`ad.concat()` drops `.uns` by default, including the `uns["log1p"]` marker
scATrans uses to detect already-log-normalized data. `de_preprocess="auto"`
still guards against double-log1p via heuristics on `.X`, but for
certainty either re-set the marker after concatenating
(`combined.uns["log1p"] = {"base": None}`) or pass `de_preprocess="none"`
explicitly. See {doc}`user_guide/standalone_de`.

## Some genes show implausibly large `logFC` (e.g. >20)

This is a known artifact of scanpy's `rank_genes_groups` log-fold-change
calculation when a gene's expression is near-zero in the reference group
(the denominator approaches zero). It is not specific to scATrans. Cross-
check any such gene against raw spliced/unspliced counts (e.g.
`scat.pl.velocity_phase_portraits`) before reporting it — see the
{doc}`tutorials/t_ec_active_transcription` tutorial for a real example, and
{doc}`statistical_guidance` for the general reporting checklist.

## I see a warning that the global unspliced fraction is > 50%

That usually indicates a technical issue (ambient RNA, mismatched
spliced/unspliced layers, or a very immature cell population) rather than a
real biological signal. `scat.qc.unspliced_global(adata)` reports this
value directly; `active_score` runs it automatically and stores it in
`adata.uns["scatrans"]["diagnostics"]`. See {doc}`installation`.

## `add_gene_features` silently produced all-`NaN` length/intron columns for some genes

`add_gene_features` does a `reindex` against `adata.var_names`, so any gene
not present in the bundled (or custom) feature table gets `NaN` and falls
back to no bias correction for that gene. If you are using a custom table,
make sure it has a `gene_name` column that matches your `var_names`
exactly. See {doc}`user_guide/gene_features`.

Huber bias correction only uses genes with **`gene_length > 0`** (and finite
intron). Missing or non-positive lengths from GTF generation are stored as
`NaN` / excluded from the fit so they cannot act as `log1p(0)` leverage
points. Partial feature tables on the simple path may be completed from the
bundled organism table without overwriting existing positive lengths.

## Design / sample-size warnings from `active_score`

When `sample_col` or `use_pseudobulk` is set, `active_score` runs
`diagnose_design` automatically, logs each design warning, and stores the
full payload under `adata.uns["scatrans"]["diagnostics"]["design"]` (including
`warnings` and `recommendations`). Inspect that block if logs are hard to
see in notebooks.

## Do I need spliced/unspliced layers at all?

No. If your data is a standard count matrix with no RNA-velocity layers,
use `differential_expression(...)` instead of `active_score(...)` — same
downstream tooling (`filter_active_genes`, enrichment, `scat.pl.*`), no
unspliced-excess term. See {doc}`user_guide/standalone_de` and the
{doc}`tutorials/t_ec_standalone_de_enrichment` tutorial.

## Can I use the bundled KEGG gene sets commercially?

Not under Apache-2.0. KEGG pathway data requires a separate commercial
license from Kanehisa Laboratories for non-academic use. To avoid the
bundled files entirely, pass an Enrichr/gseapy version explicitly, e.g.
`run_kegg(..., kegg_library="KEGG_2021")`. See {doc}`license`.

Still stuck? Open an issue on
[GitHub](https://github.com/leelieber2025/scATrans/issues).
