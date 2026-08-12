# Statistical Guidance and Reporting Checklist

Use this once you have results and need to decide what is safe to report. If
you are new to scATrans, start with {doc}`quickstart` and
{doc}`user_guide/index`. Scope: {doc}`faq`. Assumptions:
{doc}`domain_assumptions`.

### Short reporting checklist

1. Gene list from DE (`result.selected`, with padj and logFC).
2. Note capture quality (`result.regime` / `reliability`).
3. Mechanism: per-gene labels for exploration; program tables for claims;
   `calibrated` if you ran absolute placement.
4. Enrich only the DE list (not `mechanism_class` splits).
5. Record `scatrans.__version__`, DE backend, and cutoffs.

## Capability roles

| Capability | Primary columns / API | Intended use | Not intended for |
|------------|----------------------|--------------|------------------|
| DE membership | `logFC`, `p_adj`; `select_by="de"` / partition `selected` | Gene-list definition | — |
| Mechanism | residual → `transcription_support`, `mechanism_class`; `program_mechanism` / `program_mechanism_induction_matched` | Annotate DE genes (prefer program level) | Replacing DE; driving membership with `nascent_poisson_z` |
| Absolute program placement | `program_mechanism_permutation_calibrated` → `calibrated`, `null_mean`, `null_sd`, `p_perm` (+ optional `z`) | Program mean vs label-shuffle null; report `calibrated` for absolute claims (other columns optional) | Per-gene labels; reading raw means as absolute; treating `null_sd` as sampling SE of $m$ |
| Detection | `nascent_poisson_z`, `de_reproducible` via `add_nascent_score=True` | Optional active-transcription annotation | Mechanism labels; sole production filter |

## Output columns

| Output | Safe use | Do **not** use it for |
|--------|----------|------------------------|
| `adaptive_score` / `adaptive_score_pct` | **Optional** post-hoc rank that reweights the nascent leg by a data-driven reliability AUC (see `add_adaptive_score`). Report `diagnostics["reliability_auc"]`, `w_proxy`, `verdict`, and `anchor` | Treating it as FDR, calibrated probability, or a replacement for DE/`unspliced_excess_fdr`; cross-dataset numeric comparison of the score itself |
| `unspliced_excess_delta` / `unspliced_excess_residual` | Exploratory signal for **group-contrast** nascent excess (after reference γ); **mechanism** residual for `annotate_mechanism_class` | Literal transcription rates, causal claims, or equivalence to dynamical RNA velocity |
| `unspliced_excess_residual_abnorm` | Interpretable residual ranking after abundance (and optional length) normalization — demotes nuclear-retained / extreme-abundance outliers | A significance test, or assuming it restores residual reliability on steady-state velocity snapshots (kinetic limitation remains) |
| `nascent_poisson_z` / `de_reproducible` | Optional **detection** annotations (`add_nascent_score=True` / `nascent_activity_score`) | Mechanism labels / program pooling; sole production gene filter |
| `transcription_support` / `mechanism_class` / `mechanism_confidence` / `induction_confounded` | Soft **annotation** of DE-selected genes. Prefer `program_mechanism` / `program_mechanism_induction_matched` for program-level calls. Confidence is scaled by regime reliability and by the induction-confound flag | High-confidence per-gene claims; ORA on `mechanism_class` subsets; gating DE membership by the residual |
| `calibrated` / `null_mean` / `null_sd` / `p_perm` / `z` (`program_mechanism_permutation_calibrated`) | Absolute program displacement (`calibrated = observed − null_mean`). `null_sd` = SD of null means (offset stability). `p_perm` = Phipson–Smyth (floor $1/(n+1)$). `z` optional. Report frozen `de=` and `n_perm` | Reading raw program means as absolute; treating `null_sd` as sampling SE of the observed mean; gene-wise residual FDR as DE membership |
| `meta["regime"]` / `result.regime` / `qc.regime_diagnosis` | Pre-flight data-quality reliability of the nascent proxy from global unspliced fraction (U-shaped map) | Dynamic-vs-steady-state claims (not yet implemented); sole justification for production gene lists |
| `logFC`, `p_adj` (DE leg) | Standard DE reporting (with usual pseudoreplication caveats). Under **`use_mixed_model=True`**, `logFC` is **sample-mean-of-means log2FC**, not the LMM fixed-effect coefficient — see `diagnostics["mixed_model"]["logFC_method"]`. Sign discordance vs `mixedlm_coef` triggers a **warning** and is counted in `n_genes_logFC_mixedlm_sign_discordant` | Treating MixedLM `logFC` as the LMM coef, or ignoring high `n_genes_logFC_mixedlm_sign_discordant` |
| `unspliced_excess_fdr` (with `use_permutation=True`) | Exploratory significance on residual under conditional permutation | Sole production filter without DE; claims without inspecting diagnostics and replicate structure |

## Reporting checklist

1. Define gene lists with DE (`partition_de_by_mechanism` / `select_by="de"`).
   You may sort within that list by `p_adj`, `logFC`, or residual for plots —
   residual rank is not membership.
2. Report membership with DE `p_adj`. Residual FDR (`unspliced_excess_fdr`) is
   optional and not required for partition. An empty lower-level `significant`
   list from `active_score` is often expected, not a crash.
3. Call the unspliced excess a **reference-gamma group contrast**, not full
   RNA-velocity inference. Keep **detection** (`nascent_poisson_z`) separate
   from **mechanism** residual if both appear.
4. For program claims: prefer induction-matched tests when induction varies;
   for **absolute** placement report `calibrated` (and `null_mean`) from
   `program_mechanism_permutation_calibrated` (frozen `de=`, `n_perm` often 50
   in the manuscript, optional `block_col`), not the raw program mean alone.
   A transcriptional control should approach zero after calibration if the
   offset correction is working.
5. If `use_permutation=True` (gene-wise residual FDR), note the **conditional**
   shuffle and the **(1+exceed)/(n+1)** p-value convention (Phipson & Smyth —
   {doc}`references`). Program-level `p_perm` uses the same convention.
6. Cite backends and gene-set sources (scanpy / PyDESeq2 / GSEApy / GO / KEGG)
   from {doc}`references`.
7. Spot-check top genes against raw spliced/unspliced counts or phase
   portraits when you can.

Gene-wise residual permutation is off by default. Turn it on only when you need
`unspliced_excess_fdr` — not for ordinary DE membership. Program-level
calibration is a separate API.

## Quick reference (one page)

**Default path** (DE selects → mechanism partition):

| Step | Function | Key outputs |
|------|----------|-------------|
| 0. Pre-flight | `scat.qc.regime_diagnosis(adata)` (also inside partition) | `regime`, `reliability`, `message` |
| 1. Primary | `partition_de_by_mechanism(...)` | `PartitionResult`: `selected`, `gene_table`, `programs`, `regime`, `meta`; optional `programs_induction_matched` |
| 1b. Optional detection | `add_nascent_score=True` on partition (or standalone `nascent_activity_score`) | `nascent_poisson_z`, `de_reproducible` / `de_repro_frac` — **detection only**, not mechanism |
| 1c. Optional programs | `gene_sets=` / `induction_matched=True` | competitive and induction-controlled program tables (relative) |
| 1d. Absolute placement | `program_mechanism_permutation_calibrated(..., de=frozen_de)` | `calibrated`, `null_mean`, `p_perm` (program-level) |
| 2. Enrich / plot | ORA on `result.selected` (not on `mechanism_class` subsets); plots on residual / DE columns | enrichment table; figures |

**Lower-level / pure DE path:**

| Step | Function | Key outputs |
|------|----------|-------------|
| 0. Pre-flight | `recommend_workflow(...)`; with velocity layers also `scat.qc.regime_diagnosis(adata)` | workflow presets; `regime` / `reliability` / message |
| 1. Score | `active_score(...)` / `active_score_simple(...)` **or** pure DE via `differential_expression` | `all_results` / `de_results`, `adata.uns["scatrans"]` |
| 1b. Optional | `add_adaptive_score` / `add_abundance_normalized_residual` / pipeline `bias_method` & `adaptive_weighting` | `adaptive_score`, `unspliced_excess_residual_abnorm` + diagnostics |
| 1c. Optional | `annotate_mechanism_class` (pass `reliability=` from regime) / `program_mechanism` / `program_mechanism_induction_matched` / `program_mechanism_permutation_calibrated` / `threshold_sensitivity` | soft mechanism labels; program tables; absolute placement; threshold grid |
| 2. Filter | `filter_active_genes(..., select_by="de")` for production DE lists, or `preset=...` for exploratory thresholds | candidate gene list for plots / enrichment |
| 3. Enrich | `run_enrichment(candidates, gene_sets="GO_Biological_Process", adata=adata)` on DE (or detection-filtered) lists — not `mechanism_class` partitions | ORA table; cite `attrs["gene_set_info"]["provenance"]` |
| 4. Plot | `scat.pl.comet_plot(...)`, `volcano_plot(..., label_repel=True)` | `(fig, ax)`; batch export via `scat.pl.figure_export_context` or `save_all_figures` |

**Workflow presets** (via `recommend_workflow` → `WORKFLOW_PRESETS`; these
tune the lower-level `active_score` kwargs — the primary path is still
`partition_de_by_mechanism`):

- `explore` — ranking only, no permutation (fast)
- `report` — `use_permutation=True`, `n_perm=500`, `perm_de_backend="same"`
- `pseudobulk_report` — multi-replicate pseudobulk + permutation
- `nascent_focus` — `ranking_mode="nascent_excess"` (residual-only **display
  ranking**; exploratory — not a production DE list)

**Minimal paper note:** say membership is DE-defined; distinguish residual
(mechanism) from detection if both appear; cite backends and libraries; report
regime when you interpret mechanism labels. Scope: {doc}`faq`.

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
| `de_reproducible` / `de_repro_frac` | — | Spliced-side DE-reproducibility flag / fold agreement (does not change gene-list membership; genes with zero spliced fold-change are not flagged) |
| `adaptive_score` / `adaptive_score_pct` | — | Optional reliability-weighted combined score (`add_adaptive_score` / `adaptive_weighting=`) |
| `transcription_support` / `mechanism_class` / `mechanism_confidence` | — | Optional mechanism annotation (`annotate_mechanism_class`) |
| `unspliced_excess_pval` | — | One-sided permutation p-value on residual |
| `unspliced_excess_fdr` | — | BH-FDR on `unspliced_excess_pval` |

The `unspliced_excess_residual` is one-sided on positive unspliced excess and
**independent of DE significance** — genes with `p_adj` filled to 1 after
backend filters, or with weak DE, can still show positive nascent excess.
Ranking by the residual is therefore **not** a DE gene list. Prefer
`partition_de_by_mechanism` / `filter_active_genes(..., select_by="de")` for
membership. The residual is for **mechanism annotation and visualization**
within a DE-selected set, not a p-value; optional
`unspliced_excess_fdr` (when permutation is enabled) calibrates residual
magnitude under a conditional null.

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

### Built-in `significant` gene list (lower-level `active_score` only)

This subsection applies to the second return value of `active_score` /
`active_score_simple`. The primary API
{func}`~scatrans.partition_de_by_mechanism` exposes DE membership as
`result.selected` and does **not** use this conjunction.

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
(FDR on unspliced excess cannot be computed). That is expected — use
`filter_active_genes(..., select_by="de")` for production DE lists.

On low-signal data the built-in list may still be small even with
permutation. Prefer the full `all_results` table with explicit DE cutoffs
rather than assuming the built-in list matches a custom `logfc_cutoff`
override on `active_score()`.

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
