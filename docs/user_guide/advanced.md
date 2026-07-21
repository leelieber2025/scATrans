# Optional Advanced Features

These options are off by default; enable them only when the design requires them:

- `use_permutation=True`
- `bias_correction="none"`
- `show_effective_gamma=True`
- `gamma_method="robust_median"` (or `"raw"`)
- `use_mixed_model=True`
- `prioritize_velocity=True` (**deprecated** — prefer `ranking_mode="nascent_excess"`)

`diagnose_design` summarizes cell and sample counts and the global unspliced
fraction, and returns warnings plus a suggested `filter_active_genes` preset.
It runs automatically when `sample_col` or `use_pseudobulk=True` is set.
Inspect diagnostics after enabling any advanced option.

## `use_permutation=True`

**Required for the built-in `significant` list** (via `unspliced_excess_fdr`).

Adds:

- `unspliced_excess_pval` / `unspliced_excess_fdr` — permutation significance
  on the bias-corrected unspliced excess residual (one-sided, positive
  direction). **Use these for active-gene calls.**
- `active_score_pval` / `active_score_fdr` — permutation on the composite
  heuristic score (**ranking aid only; do not report as primary
  significance**).

The permutation shuffles only group labels; unspliced/spliced layers and the
reference gamma are fixed from the original labeling for speed. **This is a
conditional permutation** (conditioned on the observed velocity structure
and gamma). It is a speed/tractability tradeoff and **not an unconditional
permutation of the full data**. In small reference groups or strong batch
effects, interpret the resulting FDR with extra caution; always inspect
diagnostics and consider biological replicates.

**`perm_de_backend` (default: `"same"`)** — controls which DE method builds
the permutation null:

| Value | Behavior | When to use |
|-------|----------|-------------|
| `"same"` (**default**) | Each permutation uses the **same** DE backend as the main analysis: scanpy `de_method` / pseudobulk backend, **and MixedLM** when `use_mixed_model=True` (same `sample_col`, `mixed_model_pval`, `paired_replicates`) | **Recommended for manuscripts** — null and observed statistics match |
| `"fast"` | Always uses scanpy `t-test_overestim_var` inside permutations (faster) | Large screens / exploration only; **may bias FDR** if main analysis uses Wilcoxon, Memento, PyDESeq2, or **MixedLM**. With MixedLM, composite `active_score_pval`/`active_score_fdr` are **not valid** under `'fast'` (residual/`unspliced_excess_fdr` still are). Under `ranking_mode="nascent_excess"` the composite active_score FDR is residual-only and is **not** affected by DE-null mismatch |

```python
adata_res, significant, all_results = scat.active_score(
    adata,
    use_permutation=True,
    n_perm=500,
    perm_de_backend="same",   # default; matches main de_method
    unspliced_excess_fdr_cutoff=0.05,
)

# Faster exploration (not recommended for final FDR claims):
# ..., perm_de_backend="fast"
```

See `diagnostics["velocity"]` for the actual `gamma_method` and
`prior_weight` used.

**Realistic runtimes (heuristic mode, rough guide):** `diagnose_design` /
`recommend_workflow` return `power_summary` with an estimated duration. Rule
of thumb on an 8-core workstation:

| Genes | `n_perm` | ~Time (heuristic, parallel) |
|-------|----------|----------------------------|
| ~5k | 100 | 2–8 min |
| ~20k | 100 | 5–20 min |
| ~20k | 500 | 25–90 min |

Pseudobulk designs with few samples cap exact permutations
(`auto_adjust_n_perm=True`). `perm_de_backend="same"` with PyDESeq2 or
Memento, and `mode="advanced"`, can be **several times slower**. Use
`n_perm=100` for exploration; reserve `n_perm≥500` for final FDR claims.

## `bias_correction`

By default the package applies a
[Huber regression](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.HuberRegressor.html)
(see also {doc}`../references`) of the raw unspliced
excess delta on log(gene length) and log(intron number) and uses the
residuals as `unspliced_excess_residual`. This step can be disabled by
setting `bias_correction="none"`, in which case the raw (reference-gamma
corrected) delta is used directly.

The correction is intended to reduce technical contributions from gene
length and intron number to the unspliced excess term. Whether length or
intron number carry biological signal of interest in a given dataset is a
scientific judgment that the user must make; the correction is therefore
optional. The `bias_diagnostic_plot` function can be used to inspect the
relationship before and after correction.

### Abundance / nuclear-retention artifacts (post-hoc)

Built-in Huber length/intron correction does **not** remove extreme
abundance or nuclear-retention confounds (e.g. *MALAT1*). After
`active_score`, you can add a scale-free residual:

```python
import scatrans as scat

all_results, diag = scat.add_abundance_normalized_residual(
    all_results, method="abundance"  # or "abundance_length"
)
# column: unspliced_excess_residual_abnorm
print(diag["abundance_floor"], diag["method"])
```

Or fold the same step into the one-liner (fail-soft; diagnostics in
`meta["bias"]`):

```python
result = scat.run_default_pipeline(
    adata, groupby="condition", target_group="Disease",
    reference_group="Control", organism="mouse",
    bias_method="abundance",  # or "abundance_length"
)
```

This improves **interpretability** of residual rankings; it does not fix a
kinetically uninformative nascent proxy on steady-state snapshots.

## Reliability-adaptive nascent weight (post-hoc)

When the residual anti-correlates with induction (common on late /
steady-state velocity snapshots), a fixed nascent weight can pull the
composite below plain DE. `add_adaptive_score` estimates reliability as the
AUC of `unspliced_excess_residual` recovering an **anchor** gene set and
builds `adaptive_score` with weight
`w = clip(k * (reliability - 0.5), 0, w_max)` (defaults `k=4`, `w_max=2`):

```python
all_results, diag = scat.add_adaptive_score(all_results)  # anchor="de" default
# or end-to-end:
# all_results, diag = scat.adaptive_active_score(
#     adata, groupby="condition", target_group="Disease",
#     reference_group="Control", organism="mouse",
# )
print(diag["reliability_auc"], diag["w_proxy"], diag["verdict"], diag.get("anchor"))
```

### Reliability anchor (`anchor=` / `adaptive_anchor=`)

| Anchor | Meaning |
|--------|---------|
| `"de"` (default) | Strong DE-induced genes (`logFC >= 1` & `p_adj < 0.05`) |
| `scat.labeling_anchor(column="new_log2fc", threshold=1.0)` | Metabolic-labeling truth column |
| callable / boolean array / Series | Custom induced set |

On metabolic-labeling time courses the DE anchor can under-estimate reliability
(e.g. fast IEGs with depleted unspliced excess) and force `w_proxy=0`; a labeling
anchor recovers graded down-weighting as the proxy becomes less informative.
Pipeline form:

```python
result = scat.run_default_pipeline(
    adata, ..., adaptive_weighting=True,
    adaptive_anchor=scat.labeling_anchor("new_log2fc"),  # or "de"
)
print(result.meta.get("adaptive"))
```

`adaptive_score` is still a **heuristic rank**, not FDR. Report the
diagnostics when you use it. Core `active_score` columns and defaults are
unchanged.

## DE selects, proxy annotates

**Preferred:** {func}`~scatrans.partition_de_by_mechanism` already selects by DE
and annotates mechanism (optionally `add_nascent_score=True` for detection
columns). For the legacy pipeline path, use `select_by="de"`:

```python
result = scat.run_default_pipeline(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    select_by="de",              # candidates = DE gates only
    bias_method="abundance",     # optional residual cleanup
    annotate_mechanism=True,     # optional per-gene mechanism labels
)
# result.candidates: padj/logFC-selected
# result.all_results: still carries residual / adaptive / mechanism columns
print(result.meta["select_by"], result.meta.get("mechanism"))
```

Standalone filter:

```python
candidates = scat.filter_active_genes(all_results, select_by="de")
```

## Regime / proxy-reliability pre-flight

Before trusting nascent annotations, map the global unspliced fraction to a
dataset-level reliability scalar:

```python
r = scat.qc.regime_diagnosis(adata)
# keys: unspliced_fraction, reliability [0, 1], regime, basis, message
# regime: "ok" | "low_unspliced" | "high_unspliced"
```

| Regime | Typical meaning |
|--------|-----------------|
| `ok` | Fraction in the normal band (~10–45%); reliability ≈ 1 |
| `low_unspliced` | Little nascent signal; reliability ramps down toward 0 |
| `high_unspliced` | Possible nuclear/gDNA contamination; gamma / proxy may mis-fit |

`run_default_pipeline` always stores this in `meta["regime"]` (fail-soft if
layers are missing). Scope: **data-quality / gamma** only — not yet
dynamic-vs-steady-state (that needs a velocity-magnitude signal, pending
validation). High reliability means the proxy is not clearly corrupted; it is
not evidence that the residual outperforms DE.

## Mechanism annotation

Annotation only — does not change gene-list membership. Scale confidence with
regime reliability when velocity layers are present. Product rules: {doc}`../faq`.

### Preferred: primary workflow (+ optional detection columns)

```python
# Mechanism always uses the induction-normalized residual
res = scat.partition_de_by_mechanism(
    adata, groupby="condition", target_group="Disease", reference_group="Control",
    de="builtin", gene_sets=my_pathways,
)

# Same mechanism path, plus additive DETECTION columns (decoupled from mechanism)
res = scat.partition_de_by_mechanism(
    adata, groupby="condition", target_group="Disease", reference_group="Control",
    de="builtin",
    add_nascent_score=True,
    gene_sets=my_pathways,
)
# gene_table gains: nascent_poisson_z, de_reproducible, de_repro_frac, …
# meta["nascent_score"] records enabled / status / n_reproducible
# transcription_support / program directions are unchanged vs residual-only run
```

### Manual building blocks

```python
r = scat.qc.regime_diagnosis(adata)
# Optional: detection score (do NOT pass as residual_col for mechanism)
nz = scat.nascent_activity_score(
    adata, groupby="condition", target_group="Disease", reference_group="Control",
    sample_col="sample",
)
all_results = all_results.join(nz, how="left")

# Per-gene mechanism (low confidence by design; prefer program-level pooling)
# omit residual_col to auto-pick unspliced_excess_residual / abnorm residual
all_results, mdiag = scat.annotate_mechanism_class(
    all_results,
    reliability=r["reliability"],
)
# columns: transcription_support, mechanism_class, mechanism_confidence
# classes: transcription-driven | stabilization-driven | ambiguous |
#          unclassified_down | unknown

# Or via the pipeline (uses meta["regime"]["reliability"] automatically):
# result = scat.run_default_pipeline(..., select_by="de", annotate_mechanism=True)
# result.meta["regime"], result.meta["mechanism"]

# Program-level (threshold-free competitive Mann–Whitney on support)
prog = scat.program_mechanism(all_results, gene_sets={"IEG": ieg_list, "inflam": inflam_list})

# Threshold robustness of a DE-selected list
sens = scat.threshold_sensitivity(all_results)  # padj × logFC grid + Jaccard vs reference
```

Per-gene labels are exploratory. Prefer `program_mechanism` for
program-level transcription-versus-stabilization calls, and report
`threshold_sensitivity` rather than relying on a single DE cutoff.

## `gamma_method` and reference gamma robustness

The core unspliced excess uses a per-gene reference gamma = U_ref / S_ref
(shrunk).

- Default: `gamma_method="heuristic_shrink"` + `prior_weight=5.0` (additive
  pseudo-count shrinkage toward a global ratio).
- For small reference groups, try `gamma_method="robust_median"`: a
  **heuristic variant** of the above that uses the *median of per-gene* U/S
  ratios (instead of the global sum ratio) as the shrinkage anchor. It is
  *not* an empirical-Bayes or hierarchical method; see source/docstring for
  details.
- **`gamma_method="empirical_bayes"`** (optional, recommended for small
  reference): **hierarchical gamma estimation** using robust log-ratio
  empirical Bayes shrinkage. Prior hyperparameters are estimated once from
  the reference group (trimmed median + MAD); per-gene gammas are shrunk
  toward the shared prior on the log-ratio scale (hierarchical model across
  genes). During permutation, the **same fixed prior** is reused while
  observed ratios are recomputed from shuffled labels (conditional
  permutation preserved).
- `gamma_method="raw"` disables most shrinkage (exploratory only).

```python
adata_res, _, all_results = scat.active_score(
    adata,
    gamma_method="empirical_bayes",
    show_effective_gamma=True,  # optional: expose per-gene gamma
)
v = adata_res.uns["scatrans"]["diagnostics"]["velocity"]
print(v["gamma_prior_mean"], v["shrinkage_summary"], v["effective_gamma_stats"])
scat.pl.gamma_shrinkage_plot(all_results)  # needs gamma_shrinkage_weight column
```

`diagnose_design` recommends `empirical_bayes` (the hierarchical gamma
estimator) when the reference group is small (<80 cells).

## `show_effective_gamma=True`

Adds the column `effective_gamma` (reference-group shrunk U/S ratio) to
`adata.var` and to the results tables. Many genes will have similar values
in pure heuristic mode; advanced (moments) mode usually shows more per-gene
variation.

Example filter using the column (when present):

```python
final = scat.filter_active_genes(
    all_results,
    effective_gamma_min=0.05,   # removes genes whose gamma is dominated by the prior
    effective_gamma_max=1.0,    # optional
)
```

## `use_mixed_model=True` + `delta_variance`

Requires `sample_col` (the column identifying biological replicates/individuals).

- Replaces the simple DE statistics with LMM estimates (cell-level with
  sample as random intercept). Inference (`p_adj`, LRT, random effects)
  comes from the MixedLM fit.
- **`logFC` is not the LMM fixed-effect coefficient.** It is a
  **sample-mean-of-means log2 fold-change** (mean of per-sample means in
  target vs reference on the same scale as the fit). This avoids
  cell-count-weighted logFC when samples have unequal n_cells.
- Adds `delta_variance` (fraction of total modeled variance explained by
  condition) and `delta_var_pval` (LRT).
- `delta_variance` is always available in `all_results` when the flag is
  on; you can use it post-hoc as an additional filter.
- Use `use_delta_variance_pval=True` only if you want the LRT p-value to
  participate in the built-in `significant` mask.
- Diagnostics under `adata.uns["scatrans"]["diagnostics"]["mixed_model"]`
  include `logFC_method`, `n_genes_logFC_mixedlm_sign_discordant` (genes
  where sample-mean logFC sign disagrees with the MixedLM coef), and
  `failed_fit_rate`. When the discordant count is &gt; 0, a **logger
  warning** is emitted and those genes are excluded from the built-in
  `significant` path / active_score significance leg (`mixedlm_coef > 0`
  required). Inspect these before publication.
- **Incompatible with `use_memento_de=True`** — both are cell-level DE
  backends; enabling both raises `ValueError` (choose one).
- **With `use_permutation=True`:** default `perm_de_backend="same"` refits
  MixedLM under each shuffle so `active_score_pval`/`active_score_fdr` match
  the observed estimator. Labels are shuffled at the **sample /
  random-effect cluster** level (not per cell): a biological sample is never
  split across conditions under the null, which is required for
  hierarchical exchangeability of `(1|sample)`. Paired designs shuffle
  within subject. `perm_de_backend="fast"` deliberately uses a t-test null
  and **invalidates** those active_score FDRs (a warning is logged);
  `unspliced_excess_fdr` is unaffected. MixedLM×n_perm can be slow;
  `auto_adjust_n_perm` caps `n_perm` by the sample-level permutation space.

**Small-sample guidance:** The mixed-model path requires **≥4 biological
samples per group** and **≥6 total random-effect groups**; otherwise
`active_score(..., use_mixed_model=True)` raises `ValueError`. With fewer
replicates, use **`use_pseudobulk=True`** + `pseudobulk_de_backend="pydeseq2"`
instead (and prefer `filter_active_genes(preset="pseudobulk")` or DE
`p_adj` for significance). `recommend_workflow()` and `diagnose_design()`
surface this automatically when `sample_col` is provided.

**Paired replicate designs:** When the same `sample_col` IDs appear in both
conditions (e.g. `rep1`/`rep2` reused as labels in Disease and Control), the
default mixed-model grouping uses composite `{condition}::{sample}` random
effects so unpaired samples are not merged. For true paired/blocking designs
(same individual measured in both conditions), pass
**`paired_replicates=True`** so the raw `sample_col` IDs define the random
intercept.

```python
adata_res, significant, all_results = scat.active_score(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="mouse_id",
    use_mixed_model=True,       # only when >=4 samples/group
    paired_replicates=True,     # paired/blocking: same ID in both conditions
)
```

The mixed-model settings and median `delta_variance` are recorded in
diagnostics.

## `ranking_mode="nascent_excess"`

When ranking is residual-only, custom `weight_fc` / `weight_pval` /
`weight_unspliced` are **ignored** and forced to `(0, 0, 1)` (with a
warning). This keeps the mode name accurate: composite DE weights cannot
silently alter a residual-only ranking.

**Impact on DE / permutation mismatches:** because `active_score` then
depends only on `unspliced_excess_residual`, **DE-null mismatches do not
affect** `active_score_pval` / `active_score_fdr` under this mode — including
historical cell-level vs sample-level MixedLM shuffle issues, or
`perm_de_backend="fast"` with MixedLM/Memento/PyDESeq2 observed DE. Those
mismatches matter only under default `ranking_mode="composite"` (non-zero
`weight_fc` / `weight_pval`). Residual-only columns
(`unspliced_excess_pval` / `unspliced_excess_fdr`) never use the DE backend.

## `mode="advanced"`

Uses scVelo moments for local smoothing before computing the group-wise
gamma delta. It is still a simple reference-gamma excess calculation on the
smoothed moments, not a full stochastic or dynamical model.

Use when you have sufficient cells and want local smoothing. The function
falls back to heuristic mode on failure (`advanced_fallback=True` by
default).
