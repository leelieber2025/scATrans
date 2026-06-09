# scATrans

scATrans is a Python package for the analysis of active transcription in single-cell RNA sequencing data. It computes a composite score for genes based on the unspliced (nascent) RNA fraction, differential expression between experimental groups, and correction for technical biases related to gene length and intron number.

The package implements a dual-track approach (heuristic and advanced) with optional permutation testing for assessing the statistical significance of the resulting scores.

## Installation

```bash
# Basic installation
pip install scatrans

# With support for scVelo-based advanced mode and the gene feature generation CLI
pip install "scatrans[advanced,gene_features]" gseapy

# With support for pseudobulk differential expression using PyDESeq2
pip install "scatrans[pseudobulk]"
```

The package includes precomputed gene feature tables (gene length and intron number) for mouse. These are used by default for bias correction.

To install from source for development:

```bash
git clone https://github.com/scATrans/scatrans.git
cd scatrans
pip install -e ".[dev]"
```

If the editable install fails outside a git repository (a common issue with dynamic versioning), force a version with:

```bash
SETUPTOOLS_SCM_PRETEND_VERSION_FOR_SCATRANS=0.7.0-dev pip install -e ".[dev]"
```

**Logging.** The package uses the Python `logging` module under the name `scatrans`. The command-line tools configure basic INFO-level output. In scripts or notebooks, logging can be configured as follows:

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

**Quick data quality check (strongly recommended).** Before (or during) analysis it is useful to inspect the global unspliced fraction:

```python
import scatrans as scat
ufrac = scat.qc.unspliced_global(adata)   # logs INFO + WARNING if > 50%
```

`active_score` now automatically calls this and records the value in diagnostics (see the mode selection guide below for interpretation). A very high fraction can indicate technical problems with the velocity layers.

## Quick Start

```python
import scanpy as sc
import scatrans as scat

# Load data containing spliced/unspliced or mature/nascent layers
adata = sc.read_h5ad("your_data.h5ad")

# Attach gene features for bias correction (uses bundled mouse data by default)
adata = scat.add_gene_features(adata)

# Compute active transcription scores
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="heuristic",
    use_permutation=True,
    n_perm=200,
    show_plot=True,
)

print(significant.head())
```

> **Important**: The default `significant` list is deliberately strict. On real data it is very common to get 0 or very few genes.  
> The primary deliverable is `all_results` (full ranked table). Most users then apply custom filters on `all_results` to obtain their final biological gene list (see the detailed workflow section below).

## The active_score Function

The main analysis function is `active_score`. It accepts a large number of parameters to control data preprocessing, differential expression, velocity estimation, bias correction, statistical testing, and output.

### Core Parameters

| Parameter                  | Default          | Description |
|----------------------------|------------------|-------------|
| `adata_input`              | (required)       | AnnData object with layers for spliced/unspliced (or mature/nascent). |
| `groupby`                  | `"condition"`    | Column in `adata.obs` defining the groups to compare. |
| `target_group`             | `"GA"`           | Name of the group of interest. |
| `reference_group`          | `"Ctrl"`         | Name of the reference group. |
| `subset_col`, `subset_values` | `None`        | Optional column and values to subset cells before analysis. |
| `weight_fc`, `weight_unspliced`, `weight_pval` | `1.0` | Relative weights used when combining fold change, unspliced residual, and p-value into the final active score. |
| `pval_cutoff`, `logfc_cutoff` | `0.05`, `0.5` | Thresholds used for the **default** `significant` list (combined with `velocity_residual > 0`). These do **not** limit what is present in `all_results`. |
| `active_fdr_cutoff`        | `0.05`           | FDR threshold for the default significant list when `use_permutation=True`. |

### Differential Expression and Pseudobulk

| Parameter                  | Default              | Description |
|----------------------------|----------------------|-------------|
| `de_method`                | `"t-test_overestim_var"` | Method passed to `scanpy.tl.rank_genes_groups`. |
| `pseudobulk_de_backend`    | `"pydeseq2"`         | Backend for differential expression when `use_pseudobulk=True`. Options: `"pydeseq2"`, `"scanpy"`. |
| `use_pseudobulk`           | `False`              | Whether to aggregate to pseudobulk before analysis. |
| `sample_col`               | `None`               | Column identifying biological samples (required when `use_pseudobulk=True`, and also when `use_mixed_model=True` for the random intercept). |
| `min_cells`, `min_counts`  | `10`, `1000`         | Minimum cells and counts required to retain a pseudobulk sample. |
| `pb_x_layer`, `pb_use_total_for_x` | `"spliced"`, `True` | Controls the expression matrix used for pseudobulk aggregation. |
| `strict_pydeseq2_counts`   | `True`               | If True, raises an error when input does not appear to be raw counts for PyDESeq2. |

### Mixed Models and Delta Variance (Replicate-Aware DE)

scATrans supports an optional mixed linear model (LMM) path for the differential component of the active score. This is particularly useful for datasets with multiple biological samples (replicates) per condition/group, where treating individual cells as independent observations can lead to inflated significance (pseudoreplication).

- Enabled via `use_mixed_model=True` (requires `sample_col`).
- Uses `statsmodels.formula.api.mixedlm("y ~ C(condition)", ..., groups=sample)` on log1p library-size normalized expression (LMM on log-scale is the standard lightweight approach also used by dream/dreamlet/variancePartition).
- Returns the usual `logFC` (fixed-effect coef), `p_val`/`p_adj` (Wald test on condition term) **and**:
  - `delta_variance`: fraction of total modeled variance (fixed-effect contribution + random intercept var + residual) attributable to the condition of interest. This is directly inspired by `variancePartition` "fraction of variation explained".
  - `delta_var_pval`: likelihood-ratio test (LRT) p-value comparing the full model vs. a reduced model without the condition term.
- `delta_variance` is **always** added to `all_results` (and `adata.var`) when `use_mixed_model=True`. You can use it post-hoc to filter (e.g. `all_results[all_results["delta_variance"] > 0.05]` for genes where condition explains a non-trivial fraction of variance).
- Use `use_delta_variance_pval=True` (and `delta_var_pval_cutoff=0.05`) to make the LRT p-value part of the **significant gene mask** (supplementary filter alongside p_adj, logFC, velocity residual, etc.).
- Relation to referenced packages (implementation choice):
  - Follows the spirit of **dreamlet / variancePartition** (pseudobulk LMM + explicit variance fractions) and its 2026 Python port **dreampy** — we use a native statsmodels LMM + explicit delta variance computation without requiring heavy limma/voom reimplementation.
  - **Libra** (R) provides a menu of mixed models (LMM, NB-GLMM etc.); we expose the fast Gaussian LMM path.
  - For true count-based **negative binomial mixed models** optimized for scRNA, **NEBULA** (R primary, with Python interest) is excellent but heavier; our LMM is a practical, fast, widely-validated approximation on the log scale. For production count NBMM you can run NEBULA/dreampy separately on pseudobulks and combine results.
- Recommendation: Use `use_pseudobulk=True + pydeseq2` (or scanpy) for count-driven pseudobulk DE when you have a reasonable number of samples per group (ideally ≥5–6 after filtering). Use `use_mixed_model=True` (cell-level with sample random effect) when you want to retain single-cell resolution for the velocity signal while still properly accounting for sample structure (avoids pseudoreplication).

**Important practical guidance for small numbers of samples** (e.g. only 3 biological samples, as seen in some real datasets like the EC.h5ad example):

- With `use_pseudobulk=True` on very few samples, aggregation often causes `velocity_residual` to collapse to near-zero for the large majority of genes (most signal is averaged away). Combined with any non-trivial weight on the unspliced term, the composite `active_score` becomes very low for most genes. Permutation power is also extremely limited.
- In such low-n regimes, running with `use_pseudobulk=False` + `use_mixed_model=True` (sample as random intercept) frequently yields more informative results: the velocity excess signal remains visible at cell level for genes that truly differ, while the LMM still respects sample correlation and `delta_variance` provides an additional robust filter.
- Always inspect the actual distributions in `all_results` (`active_score`, `velocity_residual`, `active_score_fdr`, etc.) before choosing cutoffs. The `filter_active_genes(..., preset="pseudobulk")` or `preset="heuristic"` helpers can help pick starting thresholds appropriate to your run mode.

In short: pseudobulk is powerful and statistically cleaner **when you have enough samples**. With very few samples it can make the velocity component of the active score almost disappear. The mixed-model path at cell level is often the better compromise.

| Parameter                    | Default | Description |
|------------------------------|---------|-------------|
| `use_mixed_model`            | `False` | Enable LMM DE + delta variance computation (requires `sample_col`). |
| `use_delta_variance_pval`    | `False` | If True, `delta_var_pval < delta_var_pval_cutoff` is added to the significant gene criteria. |
| `delta_var_pval_cutoff`      | `0.05`  | Threshold for the optional delta variance LRT p-value filter. |
| `sample_col`                 | `None`  | Also used as the grouping factor for `(1 \| sample)` random intercept. |

Example:

```python
import scatrans as scat
adata_res, significant, all_results = scat.active_score(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="patient_id",      # or "sample", "donor" etc.
    use_mixed_model=True,
    use_delta_variance_pval=True, # include delta_var_pval in sig filter
    delta_var_pval_cutoff=0.05,
    # other params...
)
# Post-filter example using delta_variance (supplementary)
top_drivers = all_results[all_results["delta_variance"] > 0.1].head(20)
```

The mixed model settings, number of samples used for RE, and median delta_variance are recorded in `adata.uns["scatrans"]["diagnostics"]["mixed_model"]` and the top-level metadata.

### Permutation Testing

| Parameter                  | Default     | Description |
|----------------------------|-------------|-------------|
| `use_permutation`          | `False`     | Enable permutation testing to obtain gene-level p-values and FDR. |
| `perm_de_backend`          | `"fast"`    | DE backend used inside permutations (`"fast"` = scanpy, `"same"` = same as main analysis). |
| `n_perm`                   | `100`       | Number of permutations. |
| `auto_adjust_n_perm`       | `True`      | Reduce `n_perm` automatically for very small pseudobulk designs. |
| `random_seed`              | `42`        | Random seed for reproducibility of permutations and neighbor graphs. |

### Dual-Track Velocity Estimation

| Parameter                      | Default     | Description |
|--------------------------------|-------------|-------------|
| `mode`                         | `"heuristic"` | `"heuristic"` uses a global ratio method. `"advanced"` uses scVelo moments for local smoothing (requires the `advanced` extra). |
| `advanced_fallback`            | `True`      | Fall back to heuristic mode if the advanced track fails. |
| `advanced_n_neighbors`, `advanced_n_pcs` | `30`, `30` | Neighbor and PCA parameters for the advanced (scVelo moments) track. |
| `advanced_use_precomputed`     | `False`     | Use existing Mu/Ms layers instead of recomputing moments. |
| `allow_advanced_pseudobulk`    | `False`     | Permit advanced mode on pseudobulk data (experimental). |
| `advanced_recompute_neighbors` | `True`      | Force recomputation of the neighbor graph in advanced mode. |
| `prior_weight`                 | `5.0`       | Strength of the prior used when computing the reference gamma (velocity delta). |

### Choosing `mode`: heuristic vs advanced (and common pitfalls)

**Recommendation (quick decision guide):**

- Use the default `mode="heuristic"` for most analyses, especially:
  - Small cell numbers (< ~50-100 total in the two groups being compared)
  - When the reference group is small or has low unspliced counts
  - When you want the most direct "target excess relative to reference baseline" interpretation
  - For pseudobulk data (advanced is experimental and can over-smooth)

- Try `mode="advanced"` (requires `pip install "scatrans[advanced]"`) when:
  - You have a reasonable number of cells (ideally > 100-200 in the contrast)
  - You want local neighborhood smoothing (via scVelo `pp.moments`) to reduce Poisson noise in U/S counts
  - The data is not pseudobulk (or you explicitly set `allow_advanced_pseudobulk=True`)

**Key conceptual difference**
- Both modes ultimately compute a group-level "excess unspliced" delta using a reference-derived gamma.
- `heuristic` uses raw (or size-factor normalized) counts + a simple reference-group gamma (with shrinkage prior).
- `advanced` first computes locally-smoothed Mu/Ms moments, then applies the **same** reference-gamma delta formula on the smoothed values. The main practical benefit is noise reduction from moments, **not** a fundamentally different velocity model.

**Common pitfalls & limitations**
- Advanced can be unstable or slow with very few cells or when neighbor graph construction fails — it falls back gracefully when `advanced_fallback=True` (default).
- The "velocity" here is **not** scVelo's full stochastic or dynamical model (we deliberately kept it simple and group-contrast focused). It is a lightweight proxy for differential active transcription.
- Permutation testing (when enabled) fixes the velocity/Mu/Ms layers from the original labeling for speed; only the group labels are shuffled. This is documented in `adata.uns["scatrans"]["permutation_approximation_note"]`.
- Global unspliced fraction > ~50% often indicates technical issues (nuclear enrichment, gDNA contamination). The package now automatically computes and logs this (see `qc.unspliced_global`).
- Bias correction quality depends on having enough genes with reliable length/intron annotations. Check `uns["scatrans"]["diagnostics"]["bias_correction"]` after a run.
- Results are most interpretable for **clear binary biological contrasts**. Mixed cell states within a "target" group can dilute the signal.

A rich set of diagnostics (including the global unspliced fraction, bias regression coefficients, number of genes used for fitting, effective gamma per gene, etc.) is now stored in `adata.uns["scatrans"]["diagnostics"]` and a concise summary is printed to the log at the end of every `active_score` run. Always inspect these before interpreting significant genes.

### Layer Handling

| Parameter             | Default      | Description |
|-----------------------|--------------|-------------|
| `spliced_layer`       | `"spliced"`  | Name of the spliced (mature) layer. |
| `unspliced_layer`     | `"unspliced"`| Name of the unspliced (nascent) layer. |
| `de_preprocess`       | `"auto"`     | Preprocessing applied before DE (`"auto"`, `"normalize_log1p"`, or `"none"`). |

### Bias Correction and Filtering

| Parameter             | Default | Description |
|-----------------------|---------|-------------|
| `min_total_counts`    | `50`    | Minimum total (spliced + unspliced) counts required for a gene to be considered expressed. |
| `gene_type_filter`    | `None`  | If provided, restricts analysis to genes where `adata.var["gene_type"]` equals this value. |

### Output Control

| Parameter     | Default | Description |
|---------------|---------|-------------|
| `show_plot`   | `True`  | Whether to display a summary plot at the end of the analysis (delegates to `scat.pl.comet_plot`). |
| `n_jobs`      | `-1`    | Number of parallel jobs (used for permutation testing and certain DE backends). `-1` uses all available cores. |

**Return value.** The function returns a tuple `(adata_res, significant, all_results)`:
- `adata_res`: The processed AnnData object (subsetted and/or aggregated) with results written to `.var` and `.uns["scatrans"]`.
- `significant`: DataFrame of genes passing the **built-in** significance criteria (see below), sorted by active score. This list is intentionally conservative and can easily be empty (or very small) on real data.
- `all_results`: The most important output for most users — a DataFrame containing **all** genes with every computed score and intermediate value (`active_score`, `logFC`, `p_val`, `p_adj`, `velocity_residual`, `delta_variance`, etc.), sorted by active score descending. This table is what you should inspect and filter to obtain your final biological gene list.

Columns added to `adata.var` typically include `active_score`, `velocity_residual`, `logFC`, `p_val`, `p_adj`, and (when permutation is used) `active_score_pval` and `active_score_fdr`. When `use_mixed_model=True`, `delta_variance` and `delta_var_pval` are also added (and included in `all_results`).

### Recommended practical workflow: default `significant` vs. custom post-filtering on `all_results`

The built-in significance mask used for `significant` is a **strict conjunction** of several conditions:

```text
p_adj < pval_cutoff
AND logFC > logfc_cutoff
AND velocity_residual > 0
AND valid_expr (sufficient total counts)
AND active_score > 0
AND (if use_permutation) active_score_fdr < active_fdr_cutoff
AND (if use_delta_variance_pval) delta_var_pval < delta_var_pval_cutoff
```

In practice, especially when `use_permutation=True` (recommended for rigor) together with a non-trivial `velocity_residual` requirement, it is **very common** for the default `significant` to contain 0 or only a handful of genes. This does **not** mean the analysis failed — the composite `active_score` and the full `all_results` table are still highly informative.

**The recommended real-world usage pattern is therefore two-step:**

1. Run `active_score` once (with your preferred parameters, `show_plot=True` is useful) to obtain the complete ranking in `all_results`.
2. Inspect `all_results` (it is already sorted by `active_score` descending) and derive a custom `final_significant` list using thresholds that make biological sense for your experiment.

Example:

```python
import scatrans as scat

# 1. Run the software normally to obtain the full ranked results
#    (the default `sig_targets` is frequently empty or very small)
adata_res, sig_targets, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="SCI",
    reference_group="UN",
    use_permutation=True,
    n_perm=300,
    active_fdr_cutoff=0.2,
    pval_cutoff=0.05,
    logfc_cutoff=0.3,
    weight_fc=0.5,
    weight_unspliced=2.0,
    weight_pval=1.0,
    show_plot=True,
)

print(f"Number of genes in the default sig_targets: {len(sig_targets)}")   # often 0

# 2. Use the built-in helper (recommended) or your own logic on all_results
#    to obtain a final, biologically meaningful gene list.
#    The helper only applies filters for columns that are actually present.

# Calling with no arguments returns everything (sorted). To apply
# the commonly useful "recommended" filters, pass explicit values:
final_significant = scat.filter_active_genes(
    all_results,
    active_score_cutoff=55,
    pval_cutoff=0.05,
    velocity_residual_cutoff=1.0,
    logfc_cutoff=0.35,
    active_score_fdr_cutoff=0.25,      # only used if use_permutation=True
    effective_gamma_min=0.05,          # recommended default (see effective_gamma section)
    effective_gamma_max=1.0,           # optional upper bound
    # delta_variance_min=0.05,         # uncomment if you used use_mixed_model=True
)

print(f"Number of genes after custom filtering: {len(final_significant)}")

# Inspect the most relevant columns
display_cols = [
    "active_score", "logFC", "p_val", "p_adj",
    "active_score_pval", "active_score_fdr",
    "velocity_residual", "effective_gamma"
]
display_cols = [c for c in display_cols if c in final_significant.columns]
print(final_significant.head(15)[display_cols])

# Save for downstream analysis (enrichment, plotting, validation, etc.)
final_significant.to_csv("final_significant_genes_SCI.csv")
```

**Tips for choosing custom thresholds**
- Start by inspecting the distributions in `all_results.head(50)` or `all_results.describe()` for `active_score`, `velocity_residual`, `logFC`, and (if available) `active_score_fdr`.
- `velocity_residual > 0` is already required by the default significance mask. Raising the threshold (e.g. > 0.5 or > 1.0) is an effective way to retain only genes with a clear excess of nascent RNA.
- When `use_permutation=True`, `active_score_fdr` is usually the most relevant statistical filter for the final composite score. A cutoff in the 0.2–0.3 range is common for exploratory analyses (stricter than the default 0.05 used internally for `significant`).
- The condition on `effective_gamma > 0` (and not NaN) is a light guard against genes where the velocity reference gamma was unreliable (very low counts in the reference group). In most cases the existing `valid_expr` column already handles the worst offenders.
- When `use_mixed_model=True`, consider also requiring `delta_variance > X` (higher values mean the experimental condition explains a larger fraction of the gene's variance after modeling sample-level random effects).
- You can combine any columns available in `all_results`. Always guard for column existence in your script if you sometimes run with/without `use_permutation` or `use_mixed_model`.
- After obtaining `final_significant`, pass its index (or `.index.tolist()`) to `scat.run_enrichment(...)` or the plotting functions.

This two-step approach (trust the software for computation + rich ranking, apply domain knowledge for final selection) is the most common and productive way to use scATrans on real single-cell velocity datasets.

### The `filter_active_genes` helper function

Because the default `significant` list returned by `active_score` is intentionally conservative and frequently empty, most users derive their final gene list from `all_results`.

`scat.filter_active_genes` provides a convenient, documented way to do this.

Calling it with only the DataFrame (no extra keyword arguments) and no `preset` returns the **full** `all_results` table sorted by active_score descending (fully permissive, no filtering applied).

The function now supports a `preset` argument that automatically chooses sensible default thresholds for different analysis styles:

- `preset="heuristic"` (or `"single_cell"`): stricter defaults suitable for typical single-cell heuristic runs with default weights (active_score >= 55, velocity_residual > 1.0, logFC > 0.35, etc.).
- `preset="pseudobulk"`: more lenient defaults that account for the much smaller scale of `active_score` and `velocity_residual` after sample-level aggregation (active_score >= 5, velocity_residual > 0.05, logFC > 0.2, etc.).
- `preset="permissive"` (or `"none"`): no filtering.

Explicitly passed cutoff arguments always override the preset.

Example using a preset for a pseudobulk analysis:

```python
import scatrans as scat

adata_res, sig_targets, all_results = scat.active_score(...)

final_significant = scat.filter_active_genes(
    all_results,
    preset="pseudobulk",                 # chooses appropriate lenient defaults
    # You can still override individual values:
    # active_score_cutoff=8,
    # active_score_fdr_cutoff=0.25,
)
```

See the function docstring for complete details. The helper safely ignores columns that do not exist (e.g. `active_score_fdr` when `use_permutation=False`).

**Important note on scale differences**: With `use_pseudobulk=True` (or heavy re-weighting such as `weight_unspliced=2.0`), the composite `active_score` values are often much lower and `velocity_residual` values are typically near zero. The "single-cell heuristic" numbers will frequently return almost nothing in these cases. Always inspect the actual distributions in your `all_results` first:

```python
print(all_results[["active_score", "velocity_residual", "logFC", "active_score_fdr"]].describe())
```

Then choose (or let the preset choose) cutoffs appropriate for your run.

### Understanding and filtering on `effective_gamma`

`effective_gamma` is the per-gene **reference-group unspliced-to-spliced ratio** (U/S) computed with a small shrinkage prior (`prior_weight=5` by default). It is used internally to calculate the velocity delta:

```
delta = U_target - (effective_gamma × S_target)
```

It therefore reflects the transcriptional “baseline” in the reference (control) group.

**Typical range (heuristic mode, realistic data):**
- 5th percentile: ~0.07
- Median: ~0.40–0.45
- 95th percentile: ~0.74
- Overall observed range in well-powered data: roughly 0.02 – 0.85

In pure heuristic mode the values for many genes can be quite similar (close to the global reference gamma). In advanced mode (scVelo moments) they show more per-gene variation.

**Recommended filtering practice**

- **Lower bound (most important):** `effective_gamma > 0.05` (default in `filter_active_genes`).
  - Removes genes whose gamma estimate in the reference group is dominated by the prior because unspliced counts were extremely low. These produce noisy or unreliable velocity deltas.
  - On many datasets this keeps ~90–95% of genes while discarding the noisiest tail.
  - If you have very deep data you may safely use `0.03` or even `0.02`. If your reference group is small or low-quality, consider `0.08`.

- **Upper bound (optional):** `effective_gamma < 1.0` (or 0.8–0.9).
  - Genes with very high effective_gamma were already producing a lot of nascent RNA relative to mature RNA in the reference group.
  - Use this if you want to focus on genes that were relatively “quiet” in the control and became active in the target condition.

- **Best practice:** Always inspect the actual distribution in your data:

  ```python
  print(all_results["effective_gamma"].describe())
  print(all_results["effective_gamma"].quantile([0.05, 0.10, 0.90, 0.95]))
  ```

  Then choose thresholds that make sense for your biological question rather than using universal magic numbers.

The `filter_active_genes` helper exposes `effective_gamma_min` and `effective_gamma_max` with the defaults above so that the most common safe choice is easy to apply while still allowing full customization.

## Gene Feature Attachment

Bias correction requires per-gene length and intron number information.

```python
adata = scat.add_gene_features(adata, organism="mouse")
```

Priority order for supplying features:
1. `gene_features_path`: Full path to a user-provided parquet file.
2. `gene_feature_file`: Filename of a file present inside the package data directory.
3. `organism`: Selects a default bundled table (`"mouse"` or `"human"`).

Use `scat.list_available_gene_features()` to inspect the tables shipped with the installation.

The command-line tool `generate-gene-features` (installed with the `gene_features` extra) can be used to create custom tables from GTF files:

```bash
generate-gene-features --gtf genes.gtf --output mouse_features.parquet --organism mouse
```

## Plotting

The `scat.pl` submodule provides functions for visualizing analysis results. A consistent, clean style suitable for scientific publication (vector output, minimal non-data ink, readable sizes) is used by default. The style is inspired by professional single-cell visualization libraries such as OmicVerse.

```python
import scatrans as scat
scat.pl.set_style()   # Call once for good defaults (or pass parameters)
```

A temporary style can be applied with:

```python
with scat.pl.style_context(linewidth=0.8, labelsize=10):
    scat.pl.comet_plot(...)
```

### Main Plot Functions and External Axes (`ax=` / `axes=`)

All primary plotting functions accept an `ax` (or `axes`) parameter. This enables embedding plots into custom multi-panel figures — a common requirement for publication-quality composite figures.

- `comet_plot(df, top_n=12, save_path=None, title=..., point_scale=1.0, figsize=(8, 6), dpi=300, fontsize=12, cmap="coolwarm", ax=None)`
  - Recommended: log fold change vs. bias-corrected unspliced residual.

- `volcano_plot(df, top_n=10, save_path=None, logfc_cutoff=0.5, pval_cutoff=0.05, color_by="active_score", ax=None)`
  - 2D volcano plot.

- `volcano_3d(df, top_n=8, save_path=None, ..., ax=None)`
  - 3D visualization.

- `bias_diagnostic_plot(results_df, save_path=None, ..., show_regression=True, axes=None)`
  - Before/after bias correction (pass `axes=(ax1, ax2)` for two panels).

- `enrich_dotplot(enrich_df, top_n=15, save_path=None, x="GeneRatio", color_by=..., size_by=..., cmap=..., ax=None)`
  - Enrichment dot plot (highly customizable).

Additional functions (`enrich_barplot`, `active_score_rankplot`, `active_genes_heatmap`) are available for convenience. `active_score_rankplot` now has a lightweight real implementation; `active_genes_heatmap` delegates to `scanpy.pl.heatmap`.

**Example – multi-panel figure with external axes**

```python
import matplotlib.pyplot as plt
import scatrans as scat

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
scat.pl.comet_plot(all_results, top_n=8, ax=axes[0], title="Comet")
scat.pl.volcano_plot(all_results, top_n=6, ax=axes[1], title="Volcano")
plt.tight_layout()
fig.savefig("multi_panel.pdf", dpi=300, bbox_inches="tight")
```

All functions write high-quality output when `save_path` is provided (300 dpi default, tight bbox, vector-friendly fonts via `set_style`).

See the `examples/` directory:
- `examples/synthetic_active_transcription.py` — fully runnable synthetic demo (good for quick testing of the API and plotting).
- `examples/real_data_template.py` — heavily commented template showing the **recommended real-data workflow**, including pre-flight QC, diagnostics inspection, bias diagnostic plots, and enrichment. Copy and adapt it.

All plotting functions support `ax=` / `axes=` for embedding in multi-panel publication figures.

All plotting functions accept `save_path`. When supplied, the figure is saved at the requested DPI with `bbox_inches="tight"`.

## Functional Enrichment

Over-representation analysis is available via `run_enrichment`:

```python
res = scat.run_enrichment(
    gene_list=significant.index.tolist(),
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=adata.var_names.tolist(),
    pval_cutoff=0.05,
    min_size=5,
    max_size=500,
)
```

Convenience wrappers `run_kegg` and `simplify_enrichment` (Jaccard-based redundancy reduction) are also provided.

## kb_python and Custom Layer Names

When input data use `'mature'` and `'nascent'` layers (kb_python output), the package automatically detects and remaps them to the internal `'spliced'` / `'unspliced'` names. Custom layer names can also be supplied explicitly:

```python
scat.active_score(
    adata,
    spliced_layer="mature",
    unspliced_layer="nascent",
    ...
)
```

## Command-Line Interface

After installing with the `gene_features` extra, the `generate-gene-features` command is available for creating custom gene feature tables from GTF annotation files.

## License

MIT License.