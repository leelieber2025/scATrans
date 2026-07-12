# Statistical Guidance & Reporting Checklist

**Read this before writing a paper or supplement.** scATrans combines
several heuristics; not every output column is a calibrated statistical
claim.

For the full list of **implicit domain rules made explicit** (what “active”
means, p-value cutoffs, residual vs DE, GSEA ranks, gene length), see
{doc}`domain_assumptions`.

| Output | Safe use | Do **not** use it for |
|--------|----------|------------------------|
| `active_score` (0–100) | **Within-run ranking** and visualization (λ is data-adaptive; not absolute) | Cross-dataset / cross-subset numeric comparison; p-values, FDR, or "statistically significant activation" on its own |
| `unspliced_excess_delta` / `unspliced_excess_residual` | Exploratory signal for **group-contrast** nascent excess (after reference γ) | Literal transcription rates, causal claims, or equivalence to dynamical RNA velocity |
| `logFC`, `p_adj` (DE leg) | Standard DE reporting (with usual pseudoreplication caveats). Under **`use_mixed_model=True`**, `logFC` is **sample-mean-of-means log2FC**, not the LMM fixed-effect coefficient — see `diagnostics["mixed_model"]["logFC_method"]`. Sign discordance vs `mixedlm_coef` triggers a **warning** and is counted in `n_genes_logFC_mixedlm_sign_discordant` | Treating MixedLM `logFC` as the LMM coef, or ignoring high `n_genes_logFC_mixedlm_sign_discordant` |
| `unspliced_excess_fdr` (with `use_permutation=True`) | **Primary** active-gene significance filter (one-sided, conditional null) | Claims without inspecting diagnostics and replicate structure |

## Reporting checklist

1. Rank genes with `active_score`; state clearly that it is a **composite
   heuristic**, not a test statistic, and that 0–100 scores are **within-run
   relative** (data-adaptive λ — see {doc}`domain_assumptions`).
2. For significance, use DE `p_adj` and/or `unspliced_excess_fdr`
   (permutation). The built-in `significant` list is intentionally strict
   and often empty.
3. Describe the unspliced excess term as a **reference-gamma group
   contrast**, not full stochastic velocity inference (do not equate the
   default track with scVelo dynamical modeling).
4. When `use_permutation=True`, note the **conditional permutation** (labels
   shuffled; layers and γ fixed) and the **(1+exceed)/(n+1)** permutation
   *p*-value convention (Phipson & Smyth — {doc}`references`).
5. Cite backends and databases you used (scanpy / PyDESeq2 / GSEApy / GO /
   KEGG, etc.) from {doc}`references`.
6. Cross-check top hits with raw spliced/unspliced counts, phase portraits,
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
p-values alone. Soft-scale λ for each leg is **`median(positive values)/ln(2)`
from this run's gene vector** (p-value λ uses direction-positive genes only).
**Scores are therefore within-analysis relative ranks**, not absolute units:
changing the gene set (HVG), cell subset, or dataset re-estimates λ and can
move a gene from ~40 to ~70 with unchanged raw statistics. Lambdas are
recorded in `diagnostics["scoring"]`.

The residual leg remains one-sided on positive unspliced excess
(**independent of DE significance** — genes with `p_adj` filled to 1 after
backend filters, or weak DE, can still rank highly on nascent excess alone).
**Top-N by `active_score` is therefore not a DE-significant gene list**; use
`filter_active_genes` or the built-in `significant` conjunction when you need
DE gates. The score is intended **for ranking and visualization only** and
should **not** be interpreted or reported as a p-value or statistical
significance measure. Use the permutation-derived `unspliced_excess_fdr`
(when enabled) or your own post-hoc statistics for claims.

### Default filter thresholds (`preset="heuristic"`)

Single source of truth in code: `scatrans.tl.HEURISTIC_FILTER_DEFAULTS`
(`src/scatrans/tl/_common.py`). Used by `filter_active_genes(...,
preset="heuristic")`, the built-in `significant` conjunction (when
permutation ran and the run is not pseudobulk), and the default
`active_score(logfc_cutoff=...)`. Values may be tuned in a future minor
release based on scientific feedback — **always treat the installed
code dict as authoritative**.

| Key | Default | Applied as |
|-----|---------|------------|
| `logfc_cutoff` | **0.35** | `logFC >` cutoff (magnitude gate; direction via `logfc_direction`) |
| `pval_cutoff` / prefer **`padj_cutoff=`** | **0.05** | `p_adj <` cutoff (legacy name `pval_cutoff` still accepted) |
| `unspliced_excess_residual_cutoff` | **1.0** | residual `>` cutoff |
| `active_score_cutoff` | **55.0** | `active_score >=` cutoff |
| `active_score_fdr_cutoff` | **0.25** | composite-score FDR (only if permutation ran) |
| `unspliced_excess_fdr_cutoff` | **0.05** | residual FDR (only if permutation ran) |
| `effective_gamma_min` / `max` | `None` | optional γ bounds (off by default) |

After **pseudobulk** aggregation residual and score scales shrink; use
`preset="pseudobulk"` / `PSEUDOBULK_FILTER_DEFAULTS` instead
(`active_score_cutoff=5.0`, `unspliced_excess_residual_cutoff=0.05`,
`logfc_cutoff=0.2`, same FDR/p_adj defaults).

```python
from scatrans.tl import HEURISTIC_FILTER_DEFAULTS, PSEUDOBULK_FILTER_DEFAULTS
print(HEURISTIC_FILTER_DEFAULTS)
```

### Built-in `significant` gene list

When `use_permutation=True`, the built-in mask uses the same defaults as
`filter_active_genes(..., preset="heuristic")` (or pseudobulk defaults when
`is_pseudobulk`). To recover that exact list later from `all_results`, use
`filter_active_genes(all_results, preset="significant")` — it reads the
stored filter context rather than re-guessing cutoffs.

Under default **heuristic** parameters the built-in mask requires **all** of
the gates in the table above, plus:

- When MixedLM was used: also `mixedlm_coef > 0` (direction aligned with the
  tested effect / `p_adj`, not merely sample-mean `logFC`)
- `unspliced_excess_fdr` and (when present) `active_score_fdr` gates as in
  the table

Without `use_permutation=True`, the built-in `significant` list is **empty**
(FDR on unspliced excess cannot be computed). Use `all_results` +
`filter_active_genes` for custom thresholds.

On low-signal data the built-in list may still be small. Use the full table
in `all_results`, sorted by `active_score` descending. If you need different
cutoffs, pass explicit arguments to `filter_active_genes` rather than
assuming the built-in list matches a custom `logfc_cutoff` override on
`active_score()`.

### MixedLM: `logFC` vs tested coefficient

With `use_mixed_model=True`, reported `logFC` is **sample-mean-of-means
log2FC** (scanpy-style, not cell-count-weighted). Inference (`p_val` /
`p_adj`) tests the LMM fixed effect `mixedlm_coef`. When the two disagree in
sign, the fit emits a **warning** and records
`n_genes_logFC_mixedlm_sign_discordant` under
`diagnostics["mixed_model"]` (and on `de_df.attrs`). Built-in `significant`
and the composite significance leg **exclude** those discordant genes
(`mixedlm_coef > 0` required). Always inspect before interpreting
borderline genes — see {doc}`user_guide/advanced`.

After each run inspect the diagnostics:

```python
meta = adata_res.uns["scatrans"]
diag = meta["diagnostics"]
print(diag["unspliced_global_fraction"])
print(diag["bias_correction"])
print(meta.get("permutation_approximation_note"))
# Within-run soft-scale λ (why 0–100 is not absolute across runs)
print(diag.get("scoring"))  # lambda_fc, lambda_res, lambda_pval, …
# MixedLM only:
print(diag.get("mixed_model"))  # logFC_method, n_genes_logFC_mixedlm_sign_discordant, …
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
