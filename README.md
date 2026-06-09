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
- Recommendation: Use `use_pseudobulk=True + pydeseq2` (or scanpy) for count-driven pseudobulk DE. Use `use_mixed_model=True` when you specifically want cell-level modeling + sample random effect for the active transcription score.

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

# 2. Manually (or semi-automatically) derive the final gene list from all_results
#    using experiment-specific thresholds
final_significant = all_results[
    (all_results['active_score'] >= 55) & 
    (all_results['p_val'] < 0.05) & 
    (all_results['velocity_residual'] > 1.0) & 
    (all_results['logFC'] > 0.35)
].sort_values('active_score', ascending=False)

print(f"Number of genes after custom filtering: {len(final_significant)}")
print(final_significant.head(12)[['active_score', 'logFC', 'p_val', 'velocity_residual']])

# Optional: when use_mixed_model=True, you can also incorporate delta_variance
# final_significant = all_results[
#     (all_results['active_score'] >= 55) &
#     (all_results['delta_variance'] > 0.08) &   # condition explains substantial variance
#     (all_results['logFC'] > 0.35)
# ].sort_values('active_score', ascending=False)

# Save for downstream analysis (enrichment, plotting, validation, etc.)
final_significant.to_csv("final_significant_genes_SCI.csv")
```

**Tips for choosing custom thresholds**
- Start by looking at the distribution of `active_score`, `velocity_residual`, and `logFC` in `all_results.head(50)`.
- `velocity_residual > 0` is already enforced internally for the default list; raising it (e.g. > 0.5 or > 1.0) is a common way to focus on genes with clearer nascent RNA excess.
- When you enabled `use_mixed_model=True`, also consider `delta_variance > X` (higher = condition explains more of the gene's variance after accounting for sample random effect) and/or `delta_var_pval < 0.05`.
- You can combine any columns present in `all_results` (including `active_score_pval`, `active_score_fdr`, `p_adj`, etc.).
- After obtaining `final_significant`, you can feed its index into `scat.run_enrichment(...)` or plotting functions.

This two-step approach (trust the software for computation + rich ranking, apply domain knowledge for final selection) is the most common and productive way to use scATrans on real single-cell velocity datasets.

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