# Statistical Guidance & Reporting Checklist

**Read this before writing a paper or supplement.** scATrans combines
several heuristics; not every output column is a calibrated statistical
claim.

| Output | Safe use | Do **not** use it for |
|--------|----------|------------------------|
| `active_score` (0–100) | **Ranking** and visualization within one analysis | p-values, FDR, or "statistically significant activation" on its own |
| `unspliced_excess_delta` / `unspliced_excess_residual` | Exploratory signal for **group-contrast** nascent excess (after reference γ) | Literal transcription rates, causal claims, or equivalence to dynamical RNA velocity |
| `logFC`, `p_adj` (DE leg) | Standard DE reporting (with usual pseudoreplication caveats). Under **`use_mixed_model=True`**, `logFC` is **sample-mean-of-means log2FC**, not the LMM fixed-effect coefficient — see `diagnostics["mixed_model"]["logFC_method"]` | Treating MixedLM `logFC` as the LMM coef, or ignoring high `n_genes_logFC_mixedlm_sign_discordant` |
| `unspliced_excess_fdr` (with `use_permutation=True`) | **Primary** active-gene significance filter (one-sided, conditional null) | Claims without inspecting diagnostics and replicate structure |

## Reporting checklist

1. Rank genes with `active_score`; state clearly that it is a **composite
   heuristic**, not a test statistic.
2. For significance, use DE `p_adj` and/or `unspliced_excess_fdr`
   (permutation). The built-in `significant` list is intentionally strict
   and often empty.
3. Describe the unspliced excess term as a **reference-gamma group
   contrast**, not full stochastic velocity inference.
4. When `use_permutation=True`, note the **conditional permutation** (labels
   shuffled; layers and γ fixed) in methods — see diagnostics
   `permutation_approximation_note`.
5. Cross-check top hits with raw spliced/unspliced counts, phase portraits,
   and (when possible) orthogonal DE or replicate-aware models.

The simple wrappers (`active_score_simple`, `run_default_pipeline`) keep
permutation off by default so new users explore ranked tables first; enable
`use_permutation=True` explicitly when you need FDR on unspliced excess.

## Quick reference (one page)

| Step | Function | Key outputs |
|------|----------|-------------|
| 0. Pre-flight | `recommend_workflow(adata, groupby, target, ref, sample_col=...)` | `workflow_preset`, `suggested_kwargs`, `filter_preset`, `power_summary` |
| 1. Score | `active_score(..., **rec["suggested_kwargs"])` or `active_score_simple(...)` | `all_results` (rank by `active_score`), `adata.uns["scatrans"]` |
| 2. Filter | `filter_active_genes(all_results, preset=rec["filter_preset"])` | candidate gene list for plots / enrichment |
| 3. Enrich | `run_enrichment(candidates, gene_sets="GO_Biological_Process", adata=adata)` | ORA table; cite `attrs["gene_set_info"]["provenance"]` |
| 4. Plot | `scat.pl.comet_plot(...)`, `volcano_plot(..., label_repel=True)` | `(fig, ax)`; batch export via `scat.pl.figure_export_context` or `save_all_figures` |

**Workflow presets** (via `recommend_workflow` → `WORKFLOW_PRESETS`):

- `explore` — ranking only, no permutation (fast)
- `report` — `use_permutation=True`, `n_perm=500`, `perm_de_backend="same"`
- `pseudobulk_report` — multi-replicate pseudobulk + permutation
- `nascent_focus` — `ranking_mode="nascent_excess"` (residual-only ranking;
  DE-backend / DE-null mismatches cannot affect `active_score` FDR because
  `weight_fc`/`weight_pval` are forced to 0 — only residual terms matter)

**Paper checklist (minimal):**

1. State that `active_score` is a **composite heuristic rank**, not a p-value.
2. Report DE with `p_adj`; report nascent excess significance with
   `unspliced_excess_fdr` when `use_permutation=True`.
3. Describe unspliced excess as a **reference-gamma group contrast** (not
   full dynamical velocity).
4. Enrichment: name bundled library (`Hs/Mm_*_2026`) and `p_adjust_method`;
   see `src/scatrans/data/README.md`.
5. Plots: note `adjustText` is optional (`label_repel=False` to skip); label
   density via `min_label_score` / `label_fontsize`.

## Result interpretation

### Column naming (v0.9+)

Primary result columns use **unspliced / nascent excess** terminology (not
RNA velocity):

| Primary column | Legacy alias (deprecated) | Meaning |
|----------------|---------------------------|---------|
| `unspliced_excess_delta` | `velocity_delta_raw` | Raw U − γ_ref·S in target group |
| `unspliced_excess_residual` | `velocity_residual` | Bias-corrected excess residual |
| `unspliced_excess_pval` | — | One-sided permutation p-value on residual |
| `unspliced_excess_fdr` | — | BH-FDR on `unspliced_excess_pval` |

`active_score` (0–100) is a **heuristic ranking score** (weighted
soft-scaled composite of logFC + unspliced excess residual + -log p_adj).
The logFC and significance legs are **upregulation-gated** (`logFC > 0`, or
`mixedlm_coef > 0` when MixedLM is used so the gate matches what `p_adj`
tests): strongly downregulated genes do not earn composite score from
p-values alone. The p-value soft-scale λ is estimated on direction-positive
genes only. The residual leg remains one-sided on positive unspliced excess
(independent of mature-RNA DE direction). The score is intended **for ranking
and visualization only** and should **not** be interpreted or reported as a
p-value or statistical significance measure. Use the permutation-derived
`unspliced_excess_fdr` (when enabled) or your own post-hoc statistics for claims.

### Built-in `significant` gene list

When `use_permutation=True`, the built-in mask uses the same default
thresholds as `filter_active_genes(..., preset="heuristic")` (see
`HEURISTIC_FILTER_DEFAULTS` in `tl.py`). To recover that exact list later
from `all_results`, use `filter_active_genes(all_results,
preset="significant")` — it reads the stored filter context rather than
re-guessing cutoffs.

Under default parameters the built-in mask requires **all** of:

- `logFC > logfc_cutoff` (default **0.35**)
- `p_adj < pval_cutoff` (default 0.05)
- `unspliced_excess_residual > 1.0` (default residual magnitude cutoff)
- `active_score >= 55.0`
- `active_score_fdr < 0.25` (when permutation computed composite-score FDR)
- `unspliced_excess_fdr < unspliced_excess_fdr_cutoff` (default 0.05)

Without `use_permutation=True`, the built-in `significant` list is **empty**
(FDR on unspliced excess cannot be computed). Use `all_results` +
`filter_active_genes` for custom thresholds.

On low-signal data the built-in list may still be small. Use the full table
in `all_results`, sorted by `active_score` descending. If you need different
cutoffs, pass explicit arguments to `filter_active_genes` rather than
assuming the built-in list matches a custom `logfc_cutoff` override on
`active_score()`.

After each run inspect the diagnostics:

```python
meta = adata_res.uns["scatrans"]
print(meta["diagnostics"]["unspliced_global_fraction"])
print(meta["diagnostics"]["bias_correction"])
print(meta.get("permutation_approximation_note"))
```

Global unspliced fractions above ~50% frequently indicate technical issues.
Bias-correction diagnostics report the number of genes used and any
fallback behavior. The permutation note records that unspliced/spliced
layers and the reference gamma were fixed for speed while labels were
shuffled.

## Limitations

The unspliced excess term (used by the primary `active_score` workflow) is
a group-contrast proxy derived from a reference-group gamma calculation. It
is not a full stochastic or dynamical model.

The unspliced excess term is most directly applicable to binary group
contrasts. Within-group heterogeneity can reduce observed signal. When
`use_permutation=True`, labels are shuffled while unspliced/spliced layers
and the reference gamma remain fixed; this is noted in the results. Global
unspliced fractions above ~50% are reported in diagnostics. Bias correction
effectiveness depends on annotation coverage. Small replicate numbers limit
power for the unspliced excess term and FDR estimates. Mixed-model results
tend to be conservative with large between-sample variation.

When used purely as a differential expression + enrichment toolkit (via
`differential_expression`, `run_enrichment`, etc.), scATrans relies on
established backends (scanpy, PyDESeq2, etc.) whose standard statistical
caveats apply.

Always examine diagnostics, score distributions, and (when available) the
original spliced/unspliced counts before biological interpretation.
