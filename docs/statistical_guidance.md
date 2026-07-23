# Statistical Guidance and Reporting Checklist

Not every output column is a calibrated statistical claim. Domain conventions
are listed in {doc}`domain_assumptions`; product scope in {doc}`faq`.

## Capability roles

| Capability | Primary columns / API | Intended use | Not intended for |
|------------|----------------------|--------------|------------------|
| DE membership | `logFC`, `p_adj`; `select_by="de"` / partition `selected` | Gene-list definition | — |
| Mechanism | residual → `transcription_support`, `mechanism_class`; `program_mechanism` | Annotate DE genes (prefer program level) | Replacing DE; driving membership with `nascent_poisson_z` |
| Detection | `nascent_poisson_z`, `de_reproducible` via `add_nascent_score=True` | Optional active-transcription annotation | Mechanism labels; sole production filter |

## Output columns

| Output | Safe use | Do **not** use it for |
|--------|----------|------------------------|
| `adaptive_score` / `adaptive_score_pct` | **Optional** post-hoc rank that reweights the nascent leg by a data-driven reliability AUC (see `add_adaptive_score`). Report `diagnostics["reliability_auc"]`, `w_proxy`, `verdict`, and `anchor` | Treating it as FDR, calibrated probability, or a replacement for DE/`unspliced_excess_fdr`; cross-dataset numeric comparison of the score itself |
| `unspliced_excess_delta` / `unspliced_excess_residual` | Exploratory signal for **group-contrast** nascent excess (after reference γ); **mechanism** residual for `annotate_mechanism_class` | Literal transcription rates, causal claims, or equivalence to dynamical RNA velocity |
| `unspliced_excess_residual_abnorm` | Interpretable residual ranking after abundance (and optional length) normalization — demotes nuclear-retained / extreme-abundance outliers | A significance test, or assuming it restores residual reliability on steady-state velocity snapshots (kinetic limitation remains) |
| `nascent_poisson_z` / `de_reproducible` | Optional **detection** annotations (`add_nascent_score=True` / `nascent_activity_score`) | Mechanism labels / program pooling; sole production gene filter |
| `transcription_support` / `mechanism_class` / `mechanism_confidence` / `induction_confounded` | Soft **annotation** of DE-selected genes. Prefer `program_mechanism` / `program_mechanism_induction_matched` for program-level calls. Confidence is scaled by regime reliability and by the induction-confound flag | High-confidence per-gene claims; ORA on `mechanism_class` subsets; gating DE membership by the residual |
| `meta["regime"]` / `result.regime` / `qc.regime_diagnosis` | Pre-flight data-quality reliability of the nascent proxy from global unspliced fraction (U-shaped map) | Dynamic-vs-steady-state claims (not yet implemented); sole justification for production gene lists |
| `logFC`, `p_adj` (DE leg) | Standard DE reporting (with usual pseudoreplication caveats). Under **`use_mixed_model=True`**, `logFC` is **sample-mean-of-means log2FC**, not the LMM fixed-effect coefficient — see `diagnostics["mixed_model"]["logFC_method"]`. Sign discordance vs `mixedlm_coef` triggers a **warning** and is counted in `n_genes_logFC_mixedlm_sign_discordant` | Treating MixedLM `logFC` as the LMM coef, or ignoring high `n_genes_logFC_mixedlm_sign_discordant` |
| `unspliced_excess_fdr` (with `use_permutation=True`) | Exploratory significance on residual under conditional permutation | Sole production filter without DE; claims without inspecting diagnostics and replicate structure |

## Reporting checklist

1. Prefer DE-defined gene lists (`partition_de_by_mechanism` / `select_by="de"`).
   Rank within a list by DE (`p_adj`, `logFC`) or by the nascent residual
   (`unspliced_excess_residual`), not by a cross-run absolute scale.
2. For significance, use DE `p_adj` and/or `unspliced_excess_fdr`
   (permutation). The built-in `significant` list is intentionally strict
   and often empty.
3. Describe the unspliced excess term as a **reference-gamma group
   contrast**, not full stochastic velocity inference (do not equate the
   default track with scVelo dynamical modeling). Separate **detection**
   (`nascent_poisson_z`) from **mechanism** residual if both are reported.
4. When `use_permutation=True`, note the **conditional permutation** (labels
   shuffled; layers and γ fixed) and the **(1+exceed)/(n+1)** permutation
   *p*-value convention (Phipson & Smyth — {doc}`references`).
5. Cite backends and databases you used (scanpy / PyDESeq2 / GSEApy / GO /
   KEGG, etc.) from {doc}`references`.
6. Cross-check top hits with raw spliced/unspliced counts, phase portraits,
   and (when possible) orthogonal DE or replicate-aware models.

`active_score_simple` and `run_default_pipeline` leave permutation off by
default. Enable `use_permutation=True` when residual FDR is required.

## Quick reference (one page)

**Recommended primary path** (DE selects → mechanism partition):

| Step | Function | Key outputs |
|------|----------|-------------|
| 0. Pre-flight | `scat.qc.regime_diagnosis(adata)` (also inside partition) | `regime`, `reliability`, `message` |
| 1. Primary | `partition_de_by_mechanism(...)` | `PartitionResult`: `selected`, `gene_table`, `programs`, `regime`, `meta`; optional `programs_induction_matched` |
| 1b. Optional detection | `add_nascent_score=True` on partition (or standalone `nascent_activity_score`) | `nascent_poisson_z`, `de_reproducible` / `de_repro_frac` — **detection only**, not mechanism |
| 1c. Optional programs | `gene_sets=` / `induction_matched=True` | competitive and induction-controlled program tables |
| 2. Enrich / plot | ORA on `result.selected` (not on `mechanism_class` subsets); plots on residual / DE columns | enrichment table; figures |

**Lower-level / pure DE path:**

| Step | Function | Key outputs |
|------|----------|-------------|
| 0. Pre-flight | `recommend_workflow(...)`; with velocity layers also `scat.qc.regime_diagnosis(adata)` | workflow presets; `regime` / `reliability` / message |
| 1. Score | `active_score(...)` / `active_score_simple(...)` **or** pure DE via `differential_expression` | `all_results` / `de_results`, `adata.uns["scatrans"]` |
| 1b. Optional | `add_adaptive_score` / `add_abundance_normalized_residual` / pipeline `bias_method` & `adaptive_weighting` | `adaptive_score`, `unspliced_excess_residual_abnorm` + diagnostics |
| 1c. Optional | `annotate_mechanism_class` (pass `reliability=` from regime) / `program_mechanism` / `program_mechanism_induction_matched` / `threshold_sensitivity` | soft mechanism labels; program tables; threshold grid |
| 2. Filter | `filter_active_genes(..., select_by="de")` for production DE lists, or `preset=...` for exploratory thresholds | candidate gene list for plots / enrichment |
| 3. Enrich | `run_enrichment(candidates, gene_sets="GO_Biological_Process", adata=adata)` on DE (or detection-filtered) lists — not `mechanism_class` partitions | ORA table; cite `attrs["gene_set_info"]["provenance"]` |
| 4. Plot | `scat.pl.comet_plot(...)`, `volcano_plot(..., label_repel=True)` | `(fig, ax)`; batch export via `scat.pl.figure_export_context` or `save_all_figures` |

**Workflow presets** (via `recommend_workflow` → `WORKFLOW_PRESETS`):

- `explore` — ranking only, no permutation (fast)
- `report` — `use_permutation=True`, `n_perm=500`, `perm_de_backend="same"`
- `pseudobulk_report` — multi-replicate pseudobulk + permutation
- `nascent_focus` — `ranking_mode="nascent_excess"` (residual-only ranking;
  the residual FDR is DE-backend-independent, so DE-null mismatches cannot
  affect it — only residual terms matter)

**Paper checklist (minimal):** DE membership; residual vs detection vs score
roles as in the table above; cite backends/libraries; regime pre-flight if
mechanism labels used. Product scope: {doc}`faq`.

## Result interpretation

### Column naming (v0.9+)

Primary result columns use **unspliced / nascent excess** terminology (not
RNA velocity):

| Primary column | Legacy alias (deprecated) | Meaning |
|----------------|---------------------------|---------|
| `unspliced_excess_delta` | `velocity_delta_raw` | Raw U − γ_ref·S in target group |
| `unspliced_excess_residual` | `velocity_residual` | Bias-corrected excess residual |
| `unspliced_excess_residual_abnorm` | — | Optional post-hoc abundance-/length-normalized residual (`add_abundance_normalized_residual` / `bias_method=`) |
| `nascent_poisson_z` | — | Pseudobulk variance-stabilized nascent **detection** score (`nascent_activity_score` / `add_nascent_score=True`); **not** the mechanism residual |
| `dlog_unspliced` / `dlog_spliced` | — | CPM log fold-changes from the same pseudobulk contrast (diagnostic) |
| `de_reproducible` / `de_repro_frac` | — | Spliced-side DE-reproducibility flag / fold agreement (**annotation only** — never gates membership; flat genes are not flagged) |
| `adaptive_score` / `adaptive_score_pct` | — | Optional reliability-weighted combined score (`add_adaptive_score` / `adaptive_weighting=`) |
| `transcription_support` / `mechanism_class` / `mechanism_confidence` | — | Optional mechanism annotation (`annotate_mechanism_class`) |
| `unspliced_excess_pval` | — | One-sided permutation p-value on residual |
| `unspliced_excess_fdr` | — | BH-FDR on `unspliced_excess_pval` |

The `unspliced_excess_residual` is one-sided on positive unspliced excess and
**independent of DE significance** — genes with `p_adj` filled to 1 after
backend filters, or with weak DE, can still show positive nascent excess.
Ranking by the residual is therefore **not** a DE-significant gene list; use
`filter_active_genes` or the built-in `significant` conjunction when you need
DE gates. The residual is intended for **ranking and visualization**, not as a
p-value; use the permutation-derived `unspliced_excess_fdr` (when enabled) or
your own post-hoc statistics for claims.

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
| `unspliced_excess_fdr_cutoff` | **0.05** | residual FDR (only if permutation ran) |
| `effective_gamma_min` / `max` | `None` | optional γ bounds (off by default) |

After **pseudobulk** aggregation residual scales shrink; use
`preset="pseudobulk"` / `PSEUDOBULK_FILTER_DEFAULTS` instead
(`unspliced_excess_residual_cutoff=0.05`, `logfc_cutoff=0.2`, same FDR/p_adj
defaults).

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
- `unspliced_excess_fdr` gate as in the table

Without `use_permutation=True`, the built-in `significant` list is **empty**
(FDR on unspliced excess cannot be computed). Use `all_results` +
`filter_active_genes` for custom thresholds.

On low-signal data the built-in list may still be small. Use the full table
in `all_results`, sorted by `p_adj` then `logFC`. If you need different
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
**excludes** those discordant genes
(`mixedlm_coef > 0` required). Always inspect before interpreting
borderline genes — see {doc}`user_guide/advanced`.

After each run inspect the diagnostics:

```python
meta = adata_res.uns["scatrans"]
diag = meta["diagnostics"]
print(diag["unspliced_global_fraction"])
print(diag["bias_correction"])
print(meta.get("permutation_approximation_note"))
# Within-run soft-scale λ diagnostics (data-adaptive per run)
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

The unspliced excess term (the core nascent residual computed by
`active_score`) is a group-contrast proxy derived from a reference-group gamma
calculation. It is not a full stochastic or dynamical model.

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
