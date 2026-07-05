# FAQ / Troubleshooting

Common questions and errors, collected in one place. Each links to the
fuller explanation.

## "My `significant` gene list is empty" — is that a bug?

No. The built-in `significant` mask uses strict, all-of thresholds on
logFC, `p_adj`, `unspliced_excess_residual`, `active_score`, and
`unspliced_excess_fdr` (see {doc}`statistical_guidance`). On modestly
powered real designs it is often empty by design. Use the full
`all_results` table with `filter_active_genes(preset="heuristic")` or your
own cutoffs (see {doc}`user_guide/workflow`) — ranking + biological
follow-up is the intended workflow, not a `p<0.05` gate on every run.

## `ValueError` from `use_mixed_model=True`

The mixed-model path requires **≥4 biological samples per group** and **≥6
total random-effect groups**. With fewer replicates (e.g. 3 vs. 3), use
`use_pseudobulk=True` + `pseudobulk_de_backend="pydeseq2"` instead. See
{doc}`user_guide/advanced`.

## `ImportError` for `pydeseq2` / `scvelo` / `gseapy` / `memento`

These are optional extras, not installed by the base `pip install scatrans`:

```bash
pip install "scatrans[pseudobulk]"        # PyDESeq2
pip install "scatrans[advanced]"          # scVelo (mode="advanced")
pip install "scatrans[gene_features]"     # gtfparse (custom gene-feature tables)
pip install "scatrans[memento]"           # Memento (use_memento_de=True)
pip install "scatrans[gsea]" gseapy       # GSEA (run_gsea)
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
