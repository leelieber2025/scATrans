# scATrans

scATrans is a tool that helps users mine condition-wise nascent RNA relative excess signals from single-cell spliced/unspliced (or mature/nascent) data. It provides basic analysis pipelines (differential expression + velocity-layer signal integration + gene filtering + enrichment + visualization) in a concise, transparent way, and treats advanced features such as `effective_gamma`, mixed models, and `delta_variance` as opt-in exploration options for users who need them.

We explicitly acknowledge the method's limitations. We do not claim to be a gold standard. We offer one interpretable, comparable analysis perspective.

## Installation

```bash
# Basic installation
pip install scatrans

# With support for scVelo-based advanced mode and the gene feature generation CLI
pip install "scatrans[advanced,gene_features]" gseapy

# With support for pseudobulk differential expression using PyDESeq2
pip install "scatrans[pseudobulk]"
```

The package ships precomputed gene feature tables (gene length and intron number) for mouse. These are used for bias correction when available.

To install from source:

```bash
git clone https://github.com/scATrans/scatrans.git
cd scatrans
pip install -e ".[dev]"
```

**Logging.** The package logs under the name `scatrans`. You can control verbosity with:

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

**Quick data quality check (recommended).** Before analysis it is useful to inspect the global unspliced fraction:

```python
import scatrans as scat
ufrac = scat.qc.unspliced_global(adata)   # logs INFO + WARNING if > 50%
```

`active_score` automatically runs this check and records the value in diagnostics.

---

## Core Positioning

scATrans helps users extract **condition-wise nascent RNA relative excess** signals (a lightweight proxy for differential active transcription) from single-cell velocity-style data.

- **Basic pipeline (on by default):** DE + unspliced excess after reference gamma correction + optional light bias correction for length/intron number + composite scoring + gene filtering + enrichment + plotting.
- **Advanced options are opt-in:** They are powerful but add complexity and information overload. New users should start with defaults.
- **Honest by design:** The default `significant` list is deliberately strict (often empty or very small on real data). The primary deliverable is the full ranked table (`all_results`). Diagnostics are always provided so you can judge whether the signals are trustworthy in your data.

---

## Quick Start (Minimal Default Flow) — Recommended

```python
import scanpy as sc
import scatrans as scat

# 1. Load data that contains spliced/unspliced or mature/nascent layers
adata = sc.read_h5ad("your_data.h5ad")

# 2. (Optional but recommended) Attach gene features for bias correction
#    Uses the bundled mouse table by default.
adata = scat.add_gene_features(adata)

# 3. Run the analysis with default parameters — no need to worry about
#    bias_correction, effective_gamma, mixed models, etc.
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
)

# 4. The most important output for almost everyone is all_results (full ranked table)
print(all_results.head())
```

**Key point:** With default parameters you do **not** need to think about `bias_correction`, `effective_gamma`, `use_mixed_model`, or `use_permutation`. The basic pipeline is designed to "just work" and produce clean, easy-to-understand results.

Pseudobulk analysis and choice of differential expression test (e.g. Wilcoxon) **are** basic configuration options — you can freely switch them on when your experimental design calls for it (see below).

`significant` (the second return value) is intentionally conservative and frequently contains 0 or very few genes. This is normal. Use `all_results` for downstream filtering.

---

## Basic Analysis Workflow (Recommended Path)

### 3.1 Run active_score (default parameters)

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    show_plot=True,   # shows a comet plot for quick visual check
)
```

This performs:
- Differential expression between the two groups
- Velocity delta (nascent excess) using a reference-group gamma
- Light Huber bias correction on gene length + intron number (default)
- Composite active_score (0–100)
- Rich diagnostics written to `adata_res.uns["scatrans"]["diagnostics"]`

### 3.1.1 Common basic switches: pseudobulk and DE test method

These two are **standard basic options**, not advanced exploration features. You can turn them on freely depending on your data and analysis preferences:

**Pseudobulk mode** (recommended when you have multiple biological replicates per condition):

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    use_pseudobulk=True,
    sample_col="sample",                    # column identifying biological samples/individuals
    pseudobulk_de_backend="pydeseq2",       # or "scanpy"
    min_cells=5,
    min_counts=100,
    show_plot=True,
)
```

- Requires `sample_col`.
- `pseudobulk_de_backend="pydeseq2"` uses the count-based DESeq2 model (install with `pip install "scatrans[pseudobulk]"`).
- `pseudobulk_de_backend="scanpy"` + `de_method="wilcoxon"` (or `"t-test_overestim_var"`) uses scanpy's rank_genes_groups on the aggregated data.

**Switching the DE statistical test** (works for both single-cell and pseudobulk):

```python
# Use Wilcoxon rank-sum test instead of the default t-test
adata_res, significant, all_results = scat.active_score(
    ...,
    de_method="wilcoxon",                 # any method supported by scanpy.tl.rank_genes_groups
)
```

When using `use_pseudobulk=True` + `pseudobulk_de_backend="scanpy"`, the `de_method` you choose (including `"wilcoxon"`) will be used for the pseudobulk DE step.

These choices are recorded in `adata_res.uns["scatrans"]` (`de_method`, `pseudobulk_de_backend`, `use_pseudobulk`).

The `filter_active_genes` helper has a `preset="pseudobulk"` that applies more lenient default thresholds suitable after aggregation.

### 3.2 Gene filtering with filter_active_genes (core output tool)

Because the built-in `significant` list is strict, most users derive their final list from `all_results` using `filter_active_genes`.

```python
# Start permissive, then tighten based on your data
candidates = scat.filter_active_genes(
    all_results,
    active_score_cutoff=30,
    velocity_residual_cutoff=0.5,
    logfc_cutoff=0.3,
    pval_cutoff=0.05,
)

# Or use presets that choose reasonable defaults for common analysis styles
candidates = scat.filter_active_genes(all_results, preset="heuristic")
# or preset="pseudobulk" after aggregation, or preset="permissive"
```

The helper safely ignores filters for columns that do not exist (e.g. `active_score_fdr` when you did not use `use_permutation`).

### 3.3 Functional enrichment

```python
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process_2023",
    organism="mouse",
    background=adata.var_names.tolist(),
    pval_cutoff=0.05,
)
```

Convenience wrappers `run_kegg` and `simplify_enrichment` are also available.

### 3.4 Visualization

```python
import scatrans as scat

scat.pl.comet_plot(all_results, top_n=12, title="Active Drivers")
scat.pl.volcano_plot(all_results, top_n=10)
scat.pl.bias_diagnostic_plot(all_results)   # before/after bias correction view
```

All plotting functions support `ax=` / `axes=` for embedding in multi-panel publication figures and `save_path=` for high-quality output (300 dpi, tight bbox, vector-friendly fonts).

---

## Result Interpretation and Notes

### Default `significant` is often empty or very small — this is normal

The internal significance mask is a strict conjunction:
- `p_adj < pval_cutoff`
- `logFC > logfc_cutoff`
- `velocity_residual > 0`
- sufficient expression
- `active_score > 0`
- (if `use_permutation`) `active_score_fdr < active_fdr_cutoff`
- (if `use_delta_variance_pval`) `delta_var_pval < cutoff`

On real data this frequently returns 0–few genes. **Use `all_results`** and apply your own biologically motivated filters.

### Always start from `all_results`

It is already sorted by `active_score` descending and contains every gene that passed basic expression filters together with all computed values.

### Diagnostics (always inspect these)

After every run look at:

```python
meta = adata_res.uns["scatrans"]
print(meta["diagnostics"]["unspliced_global_fraction"])
print(meta["diagnostics"]["bias_correction"])
print(meta.get("permutation_approximation_note"))
```

- **unspliced_global_fraction**: > ~50% often indicates technical problems (nuclear enrichment, gDNA contamination).
- **bias_correction**: number of genes used for the fit, coefficients, whether median fallback was used.
- **permutation_approximation_note**: only present when `use_permutation=True`. Records that velocity layers/gamma were fixed for speed.

---

## Optional Advanced Features (Opt-in)

If you are interested in the following, you can turn them on manually:

- `use_permutation=True`: obtain permutation-based FDR for the composite `active_score` (the `permutation_approximation_note` will be recorded for honesty).
- `bias_correction="huber_length_intron"`: (default) explore the influence of gene length and intron number on `velocity_residual`. You can turn it off with `bias_correction="none"`.
- `show_effective_gamma=True`: expose the per-gene reference-group gamma (the U/S ratio used to compute the velocity delta). Useful when you want to filter or interpret the baseline transcription level in the reference group.
- `use_mixed_model=True`: use a mixed linear model (LMM) for replicate-aware differential expression and obtain `delta_variance` (fraction of modeled variance attributable to the condition of interest, variancePartition-style) together with a likelihood-ratio `delta_var_pval`.

These features are off by default so that the basic analysis stays simple and clear. When you enable them, please read the diagnostics and the relevant sections below carefully.

### use_permutation=True

Adds `active_score_pval` and `active_score_fdr` columns. The permutation shuffles only group labels; velocity layers and the reference gamma are computed once on the original data for speed. This approximation is documented in `permutation_approximation_note`.

### bias_correction="none"

Disables the Huber regression (and the median fallback). `velocity_residual` will equal the raw velocity delta. Useful for exploration or when you believe length/intron effects are not relevant or are already handled by your upstream processing. The `bias_diagnostic_plot` will still work and will show that before and after are essentially identical.

### show_effective_gamma=True

Adds the column `effective_gamma` (reference-group shrunk U/S ratio) to `adata.var` and to the results tables. Many genes will have similar values in pure heuristic mode; advanced (moments) mode usually shows more per-gene variation.

Recommended light guard in `filter_active_genes` (when the column is present):

```python
final = scat.filter_active_genes(
    all_results,
    effective_gamma_min=0.05,   # removes genes whose gamma is dominated by the prior
    effective_gamma_max=1.0,    # optional
)
```

### use_mixed_model=True + delta_variance

Requires `sample_col` (the column identifying biological replicates/individuals).

- Replaces the simple DE statistics with LMM estimates (cell-level with sample as random intercept).
- Adds `delta_variance` (fraction of total modeled variance explained by condition) and `delta_var_pval` (LRT).
- `delta_variance` is always available in `all_results` when the flag is on; you can use it post-hoc as an additional filter.
- Use `use_delta_variance_pval=True` only if you want the LRT p-value to participate in the built-in `significant` mask.

**Practical note on small numbers of samples:** With very few biological replicates, pseudobulk aggregation can drive most `velocity_residual` values close to zero. In such regimes the cell-level mixed-model path (`use_mixed_model=True`, `use_pseudobulk=False`) often preserves more of the velocity signal while still respecting sample structure.

The mixed-model settings and median `delta_variance` are recorded in diagnostics.

### mode="advanced"

Uses scVelo moments for local smoothing before computing the group-wise gamma delta. It is still a simple reference-gamma excess calculation on the smoothed moments, not a full stochastic or dynamical model.

- Recommended only when you have a reasonable number of cells and want noise reduction.
- Falls back to heuristic when it fails (`advanced_fallback=True` by default).
- Experimental on pseudobulk data.

---

## Limitations & Honest Interpretation

- The "velocity" signal here is a lightweight, group-contrast-focused proxy (U_target – gamma_ref × S_target). It is **not** scVelo's full dynamical model.
- Results are most interpretable for clear binary biological contrasts. Mixed cell states inside the target group dilute the signal.
- Permutation testing (when enabled) fixes the velocity layers/gamma from the original labeling. Only group labels are permuted. This is documented so you can judge its effect on FDR calibration.
- A very high global unspliced fraction (> ~50%) often indicates technical issues with the velocity layers. The package logs a warning and records the value.
- Bias correction quality depends on having enough genes with reliable length/intron annotations. Check the bias diagnostics.
- With very small numbers of biological replicates, statistical power (especially for velocity delta and permutation FDR) is limited. Inspect the distributions in `all_results` rather than trusting default cutoffs.
- `delta_variance` and mixed-model p-values can be very conservative when there is substantial sample-to-sample variation. They are tools for validity, not automatic power boosters.
- We do not claim that genes called by scATrans are "truly actively transcribed" in an absolute sense. We provide one transparent, comparable view centered on nascent RNA excess relative to a reference baseline.

Always look at the diagnostics, the actual distributions in `all_results`, and (ideally) the raw U/S counts or phase portraits for your top candidates before drawing biological conclusions.

---

## API Reference (Simplified)

### Core functions

- `active_score(...)` — main analysis. Returns `(adata_res, significant, all_results)`.
- `filter_active_genes(results_df, ...)` — post-filter the full ranked table. Supports `preset="heuristic" | "pseudobulk" | "permissive"`.

### Basic parameters (most users only need these)

These are the common "free switches" for the basic pipeline (including pseudobulk and DE method choice):

| Parameter                  | Default                  | Notes |
|----------------------------|--------------------------|-------|
| `adata_input`              | (required)               | AnnData with spliced/unspliced (or mature/nascent) layers |
| `groupby`                  | `"condition"`            | obs column defining the groups |
| `target_group` / `reference_group` | `"GA"` / `"Ctrl"` | The two conditions to compare |
| `use_pseudobulk`           | `False`                  | Set to `True` + provide `sample_col` for pseudobulk analysis |
| `sample_col`               | `None`                   | Required when `use_pseudobulk=True` (biological replicate identifier) |
| `pseudobulk_de_backend`    | `"pydeseq2"`             | `"pydeseq2"` or `"scanpy"` (when `use_pseudobulk=True`) |
| `de_method`                | `"t-test_overestim_var"` | DE method for scanpy path (e.g. `"wilcoxon"`, `"t-test"`, ...) |
| `show_plot`                | `True`                   | Show a comet plot at the end |
| `min_total_counts`         | `50`                     | Minimum total (S+U) counts to consider a gene expressed |

### Opt-in advanced / exploration parameters (see "Optional Advanced Features")

- `use_permutation`, `n_perm`, `active_fdr_cutoff`
- `bias_correction` ("huber_length_intron" or "none")
- `show_effective_gamma`
- `use_mixed_model`, `use_delta_variance_pval`, `mixed_model_pval`
- `mode` ("heuristic" or "advanced")

Full signatures and all parameters are documented in the function docstrings and the source.

### Other commonly used functions

- `add_gene_features(adata, organism="mouse", ...)` — attach length/intron info
- `list_available_gene_features()`
- `run_enrichment(...)`, `run_kegg(...)`, `simplify_enrichment(...)`
- `scat.pl.*` plotting functions (comet_plot, volcano_plot, bias_diagnostic_plot, ...)
- `scat.qc.unspliced_global(adata)`

### Layer names

The package auto-detects `mature`/`nascent` (kb_python) and remaps them internally. You can also pass `spliced_layer=...` and `unspliced_layer=...` explicitly.

---

## Gene Feature Attachment & CLI

```python
adata = scat.add_gene_features(adata, organism="mouse")
# or provide your own table
adata = scat.add_gene_features(adata, gene_features_path="my_features.parquet")
```

After installing the `gene_features` extra, the `generate-gene-features` CLI is available for creating custom tables from GTF files.

---

## Plotting Style

```python
import scatrans as scat
scat.pl.set_style()                 # once, for good defaults
# or temporary:
with scat.pl.style_context(linewidth=0.8):
    scat.pl.comet_plot(...)
```

All functions support `ax=`/`axes=` and `save_path=`.

---

## Command-Line Interface

Only the gene-feature generator is exposed as a CLI (`generate-gene-features`).

---

## License

MIT License.

---

*This README emphasizes the basic, honest, low-ceremony workflow. Advanced capabilities remain available for users who need them and are willing to read the diagnostics.*
