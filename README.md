# scATrans 🚀

**scATrans** — Single-cell Active Transcription Analysis Toolkit

A lightweight, beginner-friendly Python package to **identify genes that are currently being actively transcribed** in scRNA-seq data. It combines:

- High unspliced fraction (velocity signal)
- Differential expression between conditions
- Smart bias correction (gene length + intron number via Huber regression)
- Optional permutation-based significance testing

The result is a composite **Active Score** (0–100) that highlights "active driver" genes.

> **Current Status (v0.7+)**: The core `active_score()` engine with **Dual-Track** (`heuristic` / `advanced`) is fully functional and well-tested.  
> Full plotting submodule (`pl.comet_plot`, `volcano_3d`, etc.) and some advanced enrichment helpers are still under active development.  
> The package gracefully informs you when a feature is not yet available.

Made with ❤️ by the scATrans team (original concept by [@leelieber2025](https://github.com/leelieber2025))

---

## ✨ Key Features (Implemented)

| Feature                        | Status          | Notes |
|--------------------------------|-----------------|-------|
| **Dual-Track Engine**          | ✅ Full         | `mode="heuristic"` (fast) or `mode="advanced"` (scVelo moments + robust) |
| **Huber Bias Correction**      | ✅ Full         | Corrects for gene length & intron number bias |
| **Composite Active Scoring**   | ✅ Full         | Weighted combination of logFC + velocity residual + p-value |
| **Permutation Testing**        | ✅ Full         | Empirical p-values / FDR for the composite score |
| **PyDESeq2 + Scanpy DE**       | ✅ Full         | Supports both pseudobulk and single-cell DE backends |
| **Built-in Diagnostic Plot**   | ✅ Full         | `active_score(..., show_plot=True)` shows logFC vs velocity residual scatter |
| **Gene Feature Attachment**    | ✅ Full         | `add_gene_features()` (uses bundled mouse data or your own parquet/CSV) |
| **Metadata Recording**         | ✅ Full         | Everything stored in `adata.uns["scatrans"]` for reproducibility |
| **Functional Enrichment**      | 🟡 Partial      | `run_enrichment()` thin wrapper around `gseapy` (recommended to use gseapy directly for now) |
| **Publication Plots** (Comet, 3D Volcano, etc.) | 🔴 Planned | Use the returned DataFrame + seaborn/matplotlib for now; full module coming soon |

---

## 📦 Installation

```bash
# Basic
pip install scatrans

# With advanced mode (scVelo moments) + optional enrichment
pip install "scatrans[advanced]" gseapy
```

**From source (recommended while under active development)**:

```bash
git clone https://github.com/scATrans/scatrans.git
cd scatrans
pip install -e ".[advanced]"
```

---

## 🚀 Quick Start (Working Example)

```python
import scanpy as sc
import scatrans as scat
import pandas as pd

# 1. Load your data (must have 'spliced' and 'unspliced' layers)
adata = sc.read_h5ad("your_data.h5ad")

# 2. (Recommended) Attach gene features for bias correction
#    Uses the bundled mouse table if you don't provide one
adata = scat.add_gene_features(
    adata, 
    gene_feature_file=None   # or path to your own parquet/CSV
)

# 3. Run Active Transcription Analysis
adata_res, sig_targets, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="heuristic",           # fast default; use "advanced" for noisy data
    use_permutation=True,
    n_perm=200,                 # increase to 1000 for publication
    n_jobs=-1,
    show_plot=True              # shows nice diagnostic scatter
)

print(f"Found {len(sig_targets)} significant active driver genes")
print(sig_targets.head(10))

# 4. Inspect everything that was recorded
print(adata_res.uns["scatrans"].keys())

# 5. Functional enrichment (thin wrapper — gseapy is powerful)
my_active_genes = sig_targets.index.tolist()
go_res = scat.run_enrichment(
    gene_list=my_active_genes,
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=adata.var_names.tolist()   # highly recommended
)
print(go_res.res2d.head())   # or use gseapy's plotting functions directly
```

---

## Dual-Track Design (New in v0.7)

| Mode          | Description                                      | Speed     | Robustness | When to use                     |
|---------------|--------------------------------------------------|-----------|------------|---------------------------------|
| `heuristic`   | Fast global group-wise U/S ratio + Huber       | ⚡ Very Fast | Good      | Exploration, large datasets    |
| `advanced`    | scVelo `pp.moments()` neighborhood smoothing + same Huber | 🐢 Slower   | Higher    | Noisy data, final figures      |

```python
# Fast exploration
res = scat.active_score(adata, mode="heuristic", show_plot=True)

# More robust (uses local neighborhood information)
res = scat.active_score(adata, mode="advanced", advanced_fallback=True)
```

> **Note**: `mode="advanced"` is **experimental** but often more stable. It falls back automatically unless you set `advanced_fallback=False`.

---

## Three Usage Modes (Progressive Rigor)

1. **Default / Exploration**
   ```python
   adata_res, sig, all_res = scat.active_score(..., use_permutation=False)
   ```

2. **Fast Validation** (recommended starting point)
   ```python
   ..., use_permutation=True, perm_de_backend="fast", n_perm=100
   ```

3. **Strict Publication Mode**
   ```python
   ..., use_permutation=True, perm_de_backend="same", n_perm=1000
   ```

---

## Cell-type Specific Analysis

```python
adata_t = adata[adata.obs["cell_type"] == "T_cells"].copy()
adata_res, sig, all_res = scat.active_score(
    adata_input=adata_t,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="advanced"
)
```

---

## Working with Results

All key outputs are in `adata.var`:

- `active_score` — main composite score (0-100)
- `velocity_residual` — bias-corrected velocity signal
- `logFC`, `p_adj` — from DE
- `active_score_pval`, `active_score_fdr` — if permutation used

Significant genes are also returned as a sorted DataFrame.

You can create your own beautiful plots easily:

```python
import seaborn as sns
import matplotlib.pyplot as plt

top = all_results.head(15)
sns.scatterplot(data=top, x="logFC", y="velocity_residual", 
                size="active_score", hue="active_score", palette="viridis")
plt.title("Top Active Driver Genes")
plt.show()
```

---

## Recommended Weight Settings

```python
# Balanced (default)
weight_fc=1.0, weight_unspliced=1.0, weight_pval=1.0

# Early / bursty response (emphasize unspliced)
weight_fc=0.5, weight_unspliced=2.0, weight_pval=1.0

# Steady-state / late response (emphasize fold change)
weight_fc=2.0, weight_unspliced=0.5, weight_pval=1.0
```

Pass them directly to `active_score()`.

---

## Current Limitations & Roadmap

- Plotting module (`pl.comet_plot`, `volcano_3d`, etc.) → coming in next minor release
- Full `generate_gene_features_from_gtf()` → stub (use bundled mouse file or external tools)
- QC warning for high global unspliced fraction → not yet wired (easy to add)
- More tutorials and example notebooks → planned

We welcome contributions and feedback!

---

## Citation

If you use scATrans in your research, please cite the original method paper / preprint (to be added) and this package.

---

## License

MIT License

Happy active transcription hunting! 🚀
