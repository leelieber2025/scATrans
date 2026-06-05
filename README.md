# scATrans 🚀

**scATrans** — Single-cell Active Transcription toolkit

A lightweight, beginner-friendly Python package that helps you **quickly identify genes that are currently being actively transcribed** in scRNA-seq data (high unspliced fraction + differential expression), with built-in quality control and splicing bias correction.

Made with ❤️ by [@leelieber2025](https://github.com/leelieber2025)

---

## ✨ Key Features

- **Quality Control**: Automatically warns when the global unspliced fraction exceeds 50-60%.
- **Bias Correction**: Smart correction for gene length and intron number bias using Huber regression.
- **Advanced Scoring Engine**: Composite Active Score combining DE (logFC + p_adj) and bias-corrected velocity residual, with fully customizable weights.
- **Functional Enrichment**: Built-in `gseapy` wrapper for GO/KEGG analysis.
- **Publication-Ready Plotting**: High-quality editable vector graphics (Comet plots, 3D Volcano, Dotplots, etc.).
- **New in v0.7+**: **Dual-Track Design** (`mode="heuristic"` vs `mode="advanced"` using scVelo moments).

---

## 📦 Installation

```bash
pip install scatrans

# Recommended: install with advanced mode support
pip install "scatrans[advanced]"
```

**Gene Features** (for bias correction):

```python
import scatrans as scat
scat.add_gene_features(adata, gene_feature_file="mouse_2020A_gene_features.parquet")
```

---

## 🚀 Quick Start (Standard Workflow)

```python
import scanpy as sc
import scatrans as scat

adata = sc.read_h5ad("your_data.h5ad")

# Add gene features for bias correction
adata = scat.add_gene_features(adata, gene_feature_file="mouse_2020A_gene_features.parquet")

# Run analysis
adata_res, sig_targets, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="heuristic",           # or "advanced"
    use_permutation=True,
    n_jobs=-1,
    show_plot=True
)

# Visualization
scat.pl.comet_plot(all_results, top_n=12, save_path="Comet_Plot.pdf")
scat.pl.volcano_3d(all_results, top_n=10, save_path="3D_Volcano.pdf")

# Enrichment
my_active_genes = sig_targets.index.tolist()
go_res = scat.run_enrichment(gene_list=my_active_genes, gene_sets="GO_Biological_Process_2023", organism="mouse")
scat.pl.enrich_dotplot(go_res, title="GO Enrichment of Active Drivers", save_path="GO_Dotplot.pdf")
```

---

## Dual-Track Design (New in v0.7+)

scATrans now supports two analysis modes:

| Mode          | Description                                      | Speed     | Robustness | Recommendation                     |
|---------------|--------------------------------------------------|-----------|------------|------------------------------------|
| `heuristic`   | Original fast group-wise method                  | Very Fast | Good       | Default for exploration            |
| `advanced`    | scVelo moments + Huber bias correction           | Slower    | Higher     | Final figures, noisy data          |

**Usage**:

```python
# Fast exploration
res = scat.active_score(adata, mode="heuristic")

# More robust analysis (experimental)
res = scat.active_score(adata, mode="advanced", advanced_fallback=True)
```

> **Note**: `mode="advanced"` is currently **experimental**. It uses neighborhood smoothing via `scv.pp.moments()` on the target vs reference subset, then applies the same Huber correction as the heuristic mode.

---

## Three Usage Modes (Progressive Enhancement)

### 1. Default Mode (Exploration)
```python
adata_res, sig, all_res = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    use_permutation=False,
    show_plot=True
)
```

### 2. Fast Validation Mode
```python
adata_res, sig, all_res = scat.active_score(
    ..., 
    use_permutation=True,
    perm_de_backend="fast",
    n_perm=100
)
```

### 3. Strict Publication Mode
```python
adata_res, sig, all_res = scat.active_score(
    ..., 
    use_permutation=True,
    perm_de_backend="same",
    n_perm=1000
)
```

---

## Cell-type / Cluster-specific Analysis

```python
# Subset to specific cell type first
adata_t = adata[adata.obs["cell_type"] == "T_cells"].copy()

adata_res, sig_targets, all_results = scat.active_score(
    adata_input=adata_t,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="advanced",
    use_permutation=True,
    show_plot=True
)
```

---

## Functional Enrichment

```python
go_res = scat.run_enrichment(
    gene_list=my_active_genes,
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=adata.var_names.tolist()
)

scat.pl.enrich_dotplot(go_res, title="GO Enrichment", save_path="GO_Dotplot.pdf")
```

---

## Visualization

All plots support `save_path`, `figsize`, `dpi`, and return `(fig, ax)`.

```python
scat.pl.set_style()

scat.pl.comet_plot(all_results, top_n=12, save_path="Comet_Plot.pdf")
scat.pl.bias_diagnostic_plot(all_results, save_path="Bias_Diagnostic.pdf")
scat.pl.enrich_dotplot(go_res, save_path="GO_Dotplot.pdf")
```

**Core Plotting Functions**:
- `comet_plot()` — Signature plot (recommended)
- `volcano_plot()` / `volcano_3d()`
- `bias_diagnostic_plot()`
- `enrich_dotplot()` / `enrich_barplot()`

---

## Metadata

All settings are recorded in:

```python
adata.uns["scatrans"]
```

Includes: `mode`, `velocity_source`, `bias_correction_method`, `moments_info`, `weights`, `cutoffs`, etc.

---

## Recommended Weight Settings

- **Default (Balanced)**: `weight_fc=1.0, weight_unspliced=1.0, weight_pval=1.0`
- **Early Response**: `weight_fc=0.5, weight_unspliced=2.0, weight_pval=1.0`
- **Steady-state / Late**: `weight_fc=2.0, weight_unspliced=0.5, weight_pval=1.0`

---

## License

MIT License

Happy analyzing! 🚀
