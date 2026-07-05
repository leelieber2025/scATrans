# Optional Advanced Features

The following flags are disabled by default and should be enabled only when
required by the experimental design:

- `use_permutation=True`
- `bias_correction="none"`
- `show_effective_gamma=True`
- `gamma_method="robust_median"` (or `"raw"`)
- `use_mixed_model=True`
- `prioritize_velocity=True`

`diagnose_design` summarizes cell and sample counts plus global unspliced
fraction and returns warnings and a suggested `filter_active_genes` preset.
It runs automatically when `sample_col` or `use_pseudobulk=True` is
supplied.

Inspect the corresponding diagnostics after enabling any advanced option.

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
| `"same"` (**default**) | Each permutation uses the **same** DE backend and `de_method` as the main analysis | **Recommended for manuscripts** — null and observed statistics match |
| `"fast"` | Always uses scanpy `t-test_overestim_var` inside permutations (faster) | Large screens / exploration only; **may bias FDR** if main analysis uses Wilcoxon, Memento, or PyDESeq2 |

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

By default the package applies a Huber regression of the raw unspliced
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
  sample as random intercept).
- Adds `delta_variance` (fraction of total modeled variance explained by
  condition) and `delta_var_pval` (LRT).
- `delta_variance` is always available in `all_results` when the flag is
  on; you can use it post-hoc as an additional filter.
- Use `use_delta_variance_pval=True` only if you want the LRT p-value to
  participate in the built-in `significant` mask.

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

## `mode="advanced"`

Uses scVelo moments for local smoothing before computing the group-wise
gamma delta. It is still a simple reference-gamma excess calculation on the
smoothed moments, not a full stochastic or dynamical model.

Use when you have sufficient cells and want local smoothing. The function
falls back to heuristic mode on failure (`advanced_fallback=True` by
default).
