# FAQ / Troubleshooting

Related pages: {doc}`domain_assumptions` (conventions),
{doc}`statistical_guidance` (reporting).

## Scope and entry points

scATrans does not replace differential expression for gene discovery. Empirical
checks on steady-state and labeling-style designs indicate that the nascent
residual does not systematically recover more true positives than DE. Its role
is to annotate **mechanism** among DE genes (*transcription-driven* versus
*stabilization-driven*). Per-gene labels are modest in accuracy; prefer
program-level pooling. Check nascent-signal reliability with
{func}`~scatrans.qc.regime_diagnosis` (often uninformative under low capture or
strong 3′ bias).

| Topic | Practice |
|-------|----------|
| Gene-list membership | DE only (`partition_de_by_mechanism`, `filter_active_genes(..., select_by="de")`, or `differential_expression`). The residual does not remove DE hits. |
| Mechanism | Residual-based `transcription_support` / `mechanism_class`; soft per-gene labels. Use `gene_sets=` / `program_mechanism` for program-level inference. |
| Detection | Optional `add_nascent_score=True` → `nascent_poisson_z`, `de_reproducible`. Does not drive mechanism labels. |
| Primary entry point | {func}`~scatrans.partition_de_by_mechanism` |
| DE without mechanism | {doc}`user_guide/standalone_de` or `run_default_pipeline(..., select_by="de")` |

## Is nascent-layer scoring suitable for production gene discovery?

No. Composite ranking that depends on spliced and unspliced layers remains
experimental. Prefer DE-defined membership (table above). Differential
expression, enrichment, and plotting without nascent layers are suitable for
routine use.

## Why is the built-in `significant` list empty?

The built-in `significant` mask requires conjunction of thresholds on logFC,
`p_adj`, `unspliced_excess_residual`, and `unspliced_excess_fdr`
({doc}`statistical_guidance`). On modestly powered designs it is often empty by
design. Filter the full `all_results` table with
`filter_active_genes(preset="heuristic")` or explicit cutoffs
({doc}`user_guide/workflow`).

## Why does the nascent residual rank genes with `p_adj ≈ 1`?

The `unspliced_excess_residual` is **independent of DE significance**. Genes
that DE backends did not test or filled as neutral (e.g. PyDESeq2 independent
filtering → `padj` set to 1) can still show positive unspliced excess. Ranking
by the residual is therefore not a DE-significant gene list. Use
`filter_active_genes(..., select_by="de")`, explicit `padj_cutoff` /
`logfc_cutoff`, or the built-in `significant` conjunction. Default cutoffs:
{doc}`statistical_guidance`.

## Should I use `padj_cutoff` or `pval_cutoff`?

**Prefer `padj_cutoff=`** on `filter_active_genes`, `extract_gene_lists`, and
enrichment helpers. Legacy `pval_cutoff=` still works but is a misleading
name: when adjusted p-values exist it filters **`p_adj` / `p.adjust`**, not
raw `p_val`. Using the modern name avoids deprecation warnings and makes
reporting clearer.

## `run_gsea` returns empty / warns about mapping rate or `gene_case`

Preranked GSEA needs **signed** ranks (prefer `logFC`; one-sided score
columns are not auto-selected). It also checks symbol overlap with gene sets (same
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
(or `store_raw_counts(adata, mode="auto")` to also recover counts from
`adata.raw` when `.X` is already normalized) **before** any preprocessing. See
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
spliced/unspliced layers, nuclear enrichment / gDNA contamination, or a very
immature cell population) rather than a clean nascent-transcription signal.
`scat.qc.unspliced_global(adata)` reports the raw fraction; prefer
`scat.qc.regime_diagnosis(adata)` for a structured verdict:

```python
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
# "high_unspliced" / "low_unspliced" / "ok"; reliability in [0, 1]
```

Reliability is **U-shaped**: near 1 in a normal unspliced band (~10–45%), and
lower when the fraction is too low (weak nascent signal) or too high (gamma /
proxy may mis-fit). `partition_de_by_mechanism` always runs this pre-flight
(`result.regime` / `meta["regime"]`) and scales `mechanism_confidence`.
`run_default_pipeline` records `meta["regime"]` and applies the scale when
`annotate_mechanism=True`.

This check is data-quality / gamma reliability only. It does not classify
dynamic versus steady-state regimes. High reliability means the proxy is not
clearly corrupted; it does not imply that the residual outperforms DE.
See {doc}`installation` and {doc}`user_guide/advanced`.

## `add_gene_features` silently produced all-`NaN` length/intron columns for some genes

`add_gene_features` reindexes against `adata.var_names`. Genes absent from the
feature table receive `NaN` and skip bias correction. Custom tables require a
`gene_name` column that matches `var_names` exactly
({doc}`user_guide/gene_features`).

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

For unresolved issues, open a ticket on
[GitHub](https://github.com/leelieber2025/scATrans/issues).
