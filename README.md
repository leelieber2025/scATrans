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
| `pval_cutoff`, `logfc_cutoff` | `0.05`, `0.5` | Thresholds applied when selecting significant genes (used together with velocity residual > 0). |
| `active_fdr_cutoff`        | `0.05`           | FDR threshold applied to permutation-derived p-values when `use_permutation=True`. |

### Differential Expression and Pseudobulk

| Parameter                  | Default              | Description |
|----------------------------|----------------------|-------------|
| `de_method`                | `"t-test_overestim_var"` | Method passed to `scanpy.tl.rank_genes_groups`. |
| `pseudobulk_de_backend`    | `"pydeseq2"`         | Backend for differential expression when `use_pseudobulk=True`. Options: `"pydeseq2"`, `"scanpy"`. |
| `use_pseudobulk`           | `False`              | Whether to aggregate to pseudobulk before analysis. |
| `sample_col`               | `None`               | Column identifying biological samples (required when `use_pseudobulk=True`). |
| `min_cells`, `min_counts`  | `10`, `1000`         | Minimum cells and counts required to retain a pseudobulk sample. |
| `pb_x_layer`, `pb_use_total_for_x` | `"spliced"`, `True` | Controls the expression matrix used for pseudobulk aggregation. |
| `strict_pydeseq2_counts`   | `True`               | If True, raises an error when input does not appear to be raw counts for PyDESeq2. |

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
- `significant`: DataFrame of genes passing the significance criteria, sorted by active score.
- `all_results`: DataFrame containing all genes with computed scores and intermediate values, sorted by active score.

Columns added to `adata.var` typically include `active_score`, `velocity_residual`, `logFC`, `p_val`, `p_adj`, and (when permutation is used) `active_score_pval` and `active_score_fdr`.

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

See the `examples/` directory (in particular `examples/synthetic_active_transcription.py`) for a complete runnable demonstration that generates synthetic data, runs the full workflow, and produces multi-panel figures using `ax=`.

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