# scATrans 🚀

**scATrans** — Single-cell Active Transcription Analysis Toolkit

A lightweight, beginner-friendly Python package to **identify genes that are currently being actively transcribed** in scRNA-seq data. It combines:

- High unspliced fraction (velocity signal)
- Differential expression between conditions
- Smart bias correction (gene length + intron number via Huber regression)
- Optional permutation-based significance testing

The result is a composite **Active Score** (0–100) that highlights "active driver" genes.

> **Current Status**: Fully production-ready Dual-Track engine + complete publication-quality plotting (`pl.*`) and enrichment modules.

Made with ❤️ by the scATrans team.

---

## ✨ Installation

```bash
# Basic
pip install scatrans

# With advanced mode (scVelo) + gene feature generation CLI + enrichment
pip install "scatrans[advanced,gene_features]" gseapy
```

> **📦 Bundled data**: The package ships with `src/scatrans/data/mouse_2020A_gene_features.parquet`  
> (precomputed gene_length + intron_number for mouse). It is automatically included during installation  
> and used as the default by `add_gene_features(adata)`.

**From source (recommended for development)**:

```bash
git clone https://github.com/scATrans/scatrans.git
cd scatrans
pip install -e ".[advanced,gene_features]"
```

---

## 🚀 Quick Start

```python
import scanpy as sc
import scatrans as scat

# 1. Load your data (must have spliced/unspliced or mature/nascent layers)
adata = sc.read_h5ad("your_data.h5ad")

# 2. (Recommended) Attach gene features for bias correction
adata = scat.add_gene_features(adata)   # uses bundled mouse table by default

# 3. Run Active Transcription Analysis
adata_res, sig_targets, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="heuristic",           # or "advanced" for noisy data
    use_permutation=True,
    n_perm=200,
    show_plot=True
)

print(f"Found {len(sig_targets)} significant active driver genes")
print(sig_targets.head(10))
```

---

## 🧬 CLI: Generate Gene Features (New!)

After installing with the `gene_features` extra, you get a convenient command-line tool:

```bash
# Basic usage
generate-gene-features --gtf /path/to/genes.gtf \
                       --output mouse_gene_features.parquet \
                       --organism mouse

# For human GENCODE
generate-gene-features --gtf gencode.v49.primary_assembly.annotation.gtf \
                       --output human_gencode_v49_gene_features.parquet \
                       --organism human
```

This generates the `gene_length` + `intron_number` table needed for Huber bias correction.

You can then use it with:
```python
adata = scat.add_gene_features(adata, gene_features_path="mouse_gene_features.parquet")
```

---

## 🔄 kb_python Compatibility (Important!)

Many users generate velocity data with **kb_python** (kallisto | bustools). In these objects the layers are named:

- `'mature'`   → spliced / mature mRNA
- `'nascent'`  → unspliced / nascent pre-mRNA

**scATrans now fully supports this automatically!**

```python
# Just load your kb_python output — no extra work needed
adata = sc.read_h5ad("kb_python_velocity.h5ad")
# adata.layers will contain: 'mature', 'nascent', ...

adata_res, sig, all_res = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    # No need to specify anything — auto-detection happens
    mode="heuristic",
    show_plot=True
)
```

### What happens internally:
1. If `'spliced'` and `'unspliced'` are missing but `'mature'` + `'nascent'` exist → warning is shown and layers are automatically remapped.
2. All analysis (DE, velocity delta, bias correction, permutation, plotting) proceeds normally.
3. You can also manually control it:

```python
scat.active_score(
    adata,
    spliced_layer="mature",
    unspliced_layer="nascent",
    ...
)
```

This makes scATrans work seamlessly with:
- Standard velocyto / scVelo outputs
- kb_python lamanno/velocity mode outputs
- Any custom layer naming (just pass the parameter)

---

## 📊 Publication-Quality Plots

```python
import scatrans as scat
scat.pl.set_style()   # Call once for Nature/Cell style

# Recommended signature plot
scat.pl.comet_plot(all_results, top_n=12, save_path="Comet_Plot.pdf")

# 3D impact view
scat.pl.volcano_3d(all_results, top_n=8, save_path="Active_Volcano_3D.pdf")

# Bias correction diagnostic (unique to scATrans!)
scat.pl.bias_diagnostic_plot(all_results, save_path="Bias_Diagnostic.pdf")

# Enrichment
go_res = scat.run_enrichment(my_active_genes, gene_sets="GO_Biological_Process_2023")
scat.pl.enrich_dotplot(go_res, top_n=15, save_path="GO_Dotplot.pdf")
```

---

## Three Usage Modes (Progressive Rigor)

1. **Exploration** (fast)
   ```python
   ..., use_permutation=False
   ```

2. **Recommended starting point**
   ```python
   ..., use_permutation=True, perm_de_backend="fast", n_perm=100
   ```

3. **Publication mode**
   ```python
   ..., use_permutation=True, perm_de_backend="same", n_perm=1000
   ```

---

## Cell-type Specific Analysis

```python
adata_t = adata[adata.obs["cell_type"] == "T_cells"].copy()
adata_res, sig, all_res = scat.active_score(adata_t, groupby="condition", ...)
```

---

## Working with Results

All important columns are added to `adata.var`:
- `active_score`
- `velocity_residual` (bias-corrected)
- `logFC`, `p_adj`
- `active_score_pval`, `active_score_fdr` (when permutation used)

Significant genes are returned sorted by Active Score.

---

## Recommended Weight Settings

```python
# Early/bursty response (emphasize unspliced signal)
scat.active_score(..., weight_fc=0.5, weight_unspliced=2.0, weight_pval=1.0)

# Steady-state / late response (emphasize fold change)
scat.active_score(..., weight_fc=2.0, weight_unspliced=0.5, weight_pval=1.0)
```

---

## Functional Enrichment

```python
my_genes = sig_targets.index.tolist()
go_res = scat.run_enrichment(
    gene_list=my_genes,
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=adata.var_names.tolist()
)
scat.pl.enrich_dotplot(go_res)
```

Also available: `run_kegg()`, `simplify_enrichment()` (Jaccard redundancy reduction).

---

## Citation

If you use scATrans in your research, please cite the original method paper (to be added) and this package.

---

## License

MIT License

Happy active transcription hunting! 🚀
