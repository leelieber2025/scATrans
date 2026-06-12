# scATrans

scATrans computes a composite score that integrates differential expression with a simple reference-based measure of excess unspliced (nascent) RNA between two groups. It is designed for users working with single-cell spliced/unspliced or mature/nascent data who want to rank genes according to this combined signal.

The package supplies a basic analysis path together with several optional extensions. All methods have limitations; results should be interpreted in light of the diagnostics and the experimental design. The tool does not claim to be a gold standard or to recover "truly active" genes in an absolute sense.

## Installation

```bash
# Basic installation
pip install scatrans

# With support for scVelo-based advanced mode and the gene feature generation CLI
pip install "scatrans[advanced,gene_features]" gseapy

# With support for pseudobulk differential expression using PyDESeq2
pip install "scatrans[pseudobulk]"

# Optional: Memento (Cell 2024) as an additional cell-level DE backend
pip install "scatrans[memento]"
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

### Preserving raw counts + original spliced/unspliced layers (strongly recommended)

scATrans (especially the Memento backend and velocity/active-transcription calculations) works best when you still have access to the original raw counts and the original spliced/unspliced (or mature/nascent) matrices on as many genes as possible.

Call this **early** (right after loading + basic QC, before any HVG, normalize or log1p):

```python
import scatrans as scat

# Save raw counts + the original velocity layers for later use
scat.store_raw_counts(adata, layer="counts", save_raw=False)

# Now you can safely do the usual Scanpy preprocessing for visualization
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
# ... normalize_total, log1p, neighbors, umap, leiden ...
```

What `store_raw_counts` does:
- Saves the current `.X` (your raw counts at that moment) into `layers["counts"]`.
- If your adata contains `"spliced"` / `"unspliced"` (or `"mature"` / `"nascent"`) layers, it also saves them under `raw_spliced`, `raw_unspliced` etc. These preserved layers survive later HVG subsetting of the main object.
- `save_raw=False` is now the default (we do **not** automatically set `adata.raw` unless you explicitly ask for it with `save_raw=True`).

This way:
- Your visualization pipeline can use a small HVG + log1p `.X`.
- Later you can still run `differential_expression(..., use_memento_de=True)` or `active_score` using the full-gene raw counts and the original spliced/unspliced data from the saved layers.
- When doing enrichment, pass the gene list from the preserved full set as `universe` (see the enrichment section below for details and warnings).

See also the "Additional Capability: Standalone Differential Expression" section and the HVG-vs-velocity-layers note below.

**Impact of HVG filtering on spliced/unspliced layers (important)**

In standard Scanpy operations:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata = adata[:, adata.var.highly_variable].copy()
```

**This will also affect the spliced/unspliced layers**:

- AnnData's `.layers` (including the "spliced" and "unspliced" you stored) are automatically subset together with the genes.
- This is standard AnnData behavior and is usually **desired**, because velocity calculations (gamma estimation, unspliced excess, active_score) require the same gene set as the main expression matrix.
- If you want to use HVGs only for **visualization/clustering**, but use more genes (the full post-QC gene set or a large collection) for **differential analysis (especially Memento)**, the recommended workflow is:
  1. Immediately after loading + basic QC, call `scat.store_raw_counts(adata)` (preserves the full/large gene raw counts into the layer + .raw at that time).
  2. Make a copy for HVG + visualization: `adata_viz = adata.copy(); ... HVG on adata_viz ...`
  3. For DE, use the **original adata** (or the restored version), at which point it can still retrieve the corresponding raw counts from the layer (the number of genes depends on how many genes the adata had when you called store).
  4. If you have already performed HVG subset on the main adata, the layer will also only contain raw counts for those HVGs. In that case DE can only be performed on these genes (consistent with the principle of "user performs filtering before store").

In short: HVG subset will reduce the genes retained in spliced/unspliced, keeping it consistent with .X. If you want to use more genes for DE, you should call the DE function before HVG subset (or on a copy that has not been subset).

Optionally, if you have done HVG + log1p for visualization but later want the raw counts back in `.X` (for the genes currently selected), you can use:

```python
# Restore raw counts into .X (non-destructive by default)
adata_raw = scat.restore_raw_counts(adata, layer="counts", inplace=False)
# or inplace=True to modify the current adata
```

See also the "Additional Capability: Standalone Differential Expression" section below for the pure-DE (no velocity) use case.

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

**Key point:** The default settings run a basic analysis without requiring decisions about `bias_correction`, `effective_gamma`, `use_mixed_model`, or `use_permutation`.

Pseudobulk analysis (`use_pseudobulk`) and choice of differential expression test (`de_method`, e.g. "wilcoxon") are standard configuration options that can be selected according to the experimental design (see the section on common basic switches).

The built-in `significant` list uses a strict conjunction of thresholds and is frequently small or empty. This behavior is expected. The primary output for most users is the full ranked table returned as `all_results`.

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

Over-representation analysis is available via `run_enrichment`:

```python
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",   # or "GO_BP" — automatically resolved to the
                                         # correct organism-specific built-in (Hs/Mm_GO_..._2026)
    organism="mouse",
    adata=adata,   # if you called store_raw_counts(adata) earlier, this will
                   # automatically use the preserved full measured gene list as universe.
                   # Explicit `universe=` still takes precedence.
    pval_cutoff=0.05,
    min_size=5,
    max_size=500,
)
# New in output (for convenience + clusterProfiler compatibility):
#   - "neg_log10_padj" column
#   - res.attrs["universe_info"] with effective_universe_size, dropped_by_annotation_filter, etc.
```

**run_kegg** (convenience wrapper for KEGG pathways):

```python
kegg_res = scat.run_kegg(
    gene_list=candidates.index.tolist(),
    organism="mouse",           # or "human"
    # Defaults to the organism-specific built-in library (Hs_KEGG_2026 or Mm_KEGG_2026)
    adata=adata,   # if store_raw_counts was called earlier, this automatically uses
                   # the preserved full measured gene set as background (best practice).
    pval_cutoff=0.05,
)
```

### Default: use the package's bundled gene sets (clearest logic)

The package now **defaults to the new organism-specific built-in libraries** (4 files added to data/):

- `Hs_GO_Biological_Process_2026.txt` + `Hs_KEGG_2026.txt` for human
- `Mm_GO_Biological_Process_2026.txt` + `Mm_KEGG_2026.txt` for mouse

You only need to specify `organism=` (for KEGG especially). Base names like "GO_Biological_Process", "KEGG", "GO_BP" are automatically resolved to the correct organism + 2026 built-in file.

If you want a specific historical Enrichr version (e.g. GO_Biological_Process_2023), just write the full name — it will be treated as an Enrichr request.

```python
# KEGG — just specify organism, gets the correct built-in (Hs/Mm_2026) automatically
kegg = scat.run_kegg(gene_list=genes, organism="mouse")

# GO — base name is enough (automatically resolved to Hs/Mm_GO_..._2026)
go = scat.run_enrichment(
    gene_list=genes,
    gene_sets="GO_Biological_Process",   # or "GO_BP"
    organism="mouse",
    # universe should be your full measured gene set (the package will use
    # the list saved by store_raw_counts if you pass adata= or the old universe=)
    universe=background,
)
```

### Using original Enrichr versions

To use a specific historical Enrichr/gseapy version, **just write the exact gene set name** (the one that includes the year/version). The system will detect that it is an Enrichr-style versioned library and load it directly via gseapy.

```python
# Specific Enrichr version for KEGG — just write the name
kegg_2021 = scat.run_kegg(
    genes, organism="mouse",
    kegg_library="KEGG_2021"     # or KEGG_2019, KEGG_2016, etc.
)

# Specific version for GO
go_2021 = scat.run_enrichment(
    genes,
    gene_sets="GO_Biological_Process_2021",  # 2023, 2021, 2019, 2018, 2017...
    # `background` or `universe` should be the full set of genes considered
    # in your experiment (not limited to HVGs).
    universe=background,
)

# Even without the year in some cases, but the year-containing names are the clearest signal
```

`gene_set_source` remains as an **explicit override** if you ever need to force one side:

- `gene_set_source="scatrans"` → force the bundled version
- `gene_set_source="enrichr"` → force the gseapy/Enrichr path

Discovery (what bundled sets are available):

```python
print(scat.list_bundled_gene_sets())
# ['Hs_GO_Biological_Process_2026.txt', 'Hs_KEGG_2026.txt', 'Mm_GO_Biological_Process_2026.txt', 'Mm_KEGG_2026.txt', ...]
```

**Motivation**: Default should be the package's own sets with almost no extra parameters (only organism for KEGG). Choosing an Enrichr version should be as simple as writing the gene set name you want.

**Adding your own sets**: Drop `.gmt` files into `src/scatrans/data/`. See `src/scatrans/data/README.md`.

**simplify_enrichment** (reduce redundant terms using Jaccard similarity):

```python
simplified = scat.simplify_enrichment(
    enrich_res,
    similarity_cutoff=0.5,
    min_count=3,
    method="jaccard",           # currently the only supported method
)
```

`run_kegg` and `simplify_enrichment` are convenience wrappers around the core `run_enrichment` function.

### 3.4 Visualization

```python
import scatrans as scat

scat.pl.comet_plot(all_results, top_n=12, title="Active Drivers")
scat.pl.volcano_plot(all_results, top_n=10, label_genes=["YourGene1", "YourGene2"])  # or just top_n
scat.pl.bias_diagnostic_plot(all_results)   # before/after bias correction view
```

All plotting functions support `ax=` / `axes=` for embedding in multi-panel publication figures and `save_path=` for high-quality output (300 dpi, tight bbox, vector-friendly fonts).

---

## Helper Functions

### diagnose_design

`diagnose_design` is a utility to analyze your experimental design (cell counts, number of biological replicates, global unspliced fraction, etc.) and give practical recommendations and warnings. It is especially useful when you have a small number of samples or are deciding whether to use pseudobulk or mixed-model analysis.

**Basic usage:**

```python
import scanpy as sc
import scatrans as scat

adata = sc.read_h5ad("your_velocity_data.h5ad")

diag = scat.diagnose_design(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample"          # highly recommended if you have biological replicates
)

# Inspect the results
print("Warnings:")
for w in diag["warnings"]:
    print("  -", w)

print("\nRecommendations:")
for r in diag["recommendations"]:
    print("  -", r)

print("\nSuggested preset for filter_active_genes:", diag.get("suggested_preset"))
```

**What it returns:**

A dictionary containing:
- `n_cells_target`, `n_cells_reference`
- `n_samples_target`, `n_samples_reference` (when `sample_col` is provided)
- `unspliced_global_fraction`
- `warnings`: list of strings (e.g. low power warnings)
- `recommendations`: list of strings
- `suggested_preset`: "heuristic", "pseudobulk", or None

**Automatic usage:**

`diagnose_design` is automatically called inside `active_score(...)` whenever you pass `sample_col` or set `use_pseudobulk=True`. You will see its output in the log.

### run_kegg and simplify_enrichment

These are convenience functions built on top of `run_enrichment`.

**run_kegg** – Run KEGG pathway enrichment directly:

```python
kegg_res = scat.run_kegg(
    gene_list=significant.index.tolist(),   # or from all_results
    organism="mouse",                       # "mouse" or "human"
    adata=adata,   # preferred: if you called store_raw_counts earlier, this will
                   # automatically use the preserved full measured gene set.
    pval_cutoff=0.05,
    min_size=5,
    max_size=500,
    return_all=False,                       # False = only significant terms
    # Defaults to the organism-specific built-in (Hs/Mm_KEGG_2026).
    # To use the original Enrichr version instead: kegg_library="KEGG_2026"
)

print(kegg_res[["Term", "p.adjust", "Count"]].head())
```

The `gene_set_source` parameter (default `"scatrans"`) controls which KEGG set is used.
See the section "Choosing gene sets explicitly with `gene_set_source`" above for full details
and examples for both GO and KEGG.

**simplify_enrichment** – Remove redundant terms from enrichment results (Jaccard-based):

```python
# After obtaining an enrichment result
simplified = scat.simplify_enrichment(
    kegg_res,                    # or enrich_res from run_enrichment
    similarity_cutoff=0.5,       # Jaccard similarity threshold
    min_count=3,                 # minimum number of genes in a term
    by="p.adjust",               # column to sort by
    ascending=True,
    method="jaccard",            # currently only "jaccard" is supported
)

print(f"Reduced from {len(kegg_res)} to {len(simplified)} terms")
print(simplified[["Term", "p.adjust", "Count"]].head())
```

This function looks for common gene list columns (`Genes`, `Lead_genes`, etc.) automatically.

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

The following options can be enabled when relevant to the analysis goals:

- `use_permutation=True`: compute a permutation-based FDR for the composite score. When enabled, a note describing the approximation (velocity layers and reference gamma are fixed from the original labeling) is stored in the results.
- `bias_correction="none"`: disable the length/intron correction on the velocity delta. The raw delta is then used directly as `velocity_residual`.
- `show_effective_gamma=True`: include the per-gene reference-group U/S ratio (used internally for the delta calculation) in the output tables.
- `use_mixed_model=True`: fit a mixed linear model with sample as random intercept and obtain `delta_variance` (fraction of modeled variance attributed to condition) along with a likelihood-ratio p-value.
- `prioritize_velocity=True`: convenience flag that increases the relative weight given to the velocity_residual (nascent excess) term while decreasing the weights on the differential expression terms. This option is provided for analyses whose primary goal is to highlight differences in unspliced abundance after reference correction. It is documented under advanced features because it changes the balance of the composite score.

A helper function `diagnose_design` is available to summarize cell and sample counts, global unspliced fraction, and to surface warnings and suggestions before or between runs of `active_score`.

These options are not enabled by default. When used, the corresponding diagnostics should be examined.

### use_permutation=True

Adds `active_score_pval` and `active_score_fdr` columns. The permutation shuffles only group labels; velocity layers and the reference gamma are computed once on the original data for speed. This approximation is documented in `permutation_approximation_note`.

### bias_correction

By default the package applies a Huber regression of the raw velocity delta on log(gene length) and log(intron number) and uses the residuals as `velocity_residual`. This step can be disabled by setting `bias_correction="none"`, in which case the raw (reference-gamma corrected) delta is used directly.

The correction is intended to reduce technical contributions from gene length and intron number to the unspliced excess term. Whether length or intron number carry biological signal of interest in a given dataset is a scientific judgment that the user must make; the correction is therefore optional. The `bias_diagnostic_plot` function can be used to inspect the relationship before and after correction.

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

Recommended only when you have a reasonable number of cells and want noise reduction. Falls back to heuristic when it fails (`advanced_fallback=True` by default). Experimental on pseudobulk data.

---

## Limitations

The method implements a composite score based on a simplified, reference-group gamma excess calculation for the unspliced layer together with standard differential expression statistics.

- The unspliced excess term is a group-contrast proxy and is not equivalent to scVelo's full stochastic or dynamical models.
- The approach is most straightforward to interpret for clear binary group contrasts. Heterogeneity within the target group can reduce the observed signal.
- When `use_permutation=True`, only the group labels are permuted; the velocity layers and reference gamma are computed once on the original data for computational efficiency. This approximation is recorded in the results metadata.
- Global unspliced fractions above ~50% are flagged by the package, as they may indicate technical issues affecting the velocity layers.
- Bias correction performance depends on the number and quality of genes with length and intron annotations.
- With small numbers of biological replicates, power for the velocity component and for permutation-based FDR is limited. Users should examine the full distributions in `all_results`.
- `delta_variance` and the associated mixed-model p-values tend to be conservative in the presence of substantial between-sample variation.

Users should examine the diagnostics stored under `adata.uns["scatrans"]["diagnostics"]`, the distributions of scores in the returned tables, and (where possible) the raw spliced/unspliced counts for candidate genes before biological interpretation.

---

## API Reference (Simplified)

### Core functions

- `active_score(...)` — main analysis for active transcription from velocity data. Returns `(adata_res, significant, all_results)`.
- `differential_expression(...)` — standalone DE (no velocity data required). Supports the same backends as `active_score` (including optional Memento). Returns `(adata, results_df)`.
- `filter_active_genes(results_df, ...)` — post-filter the full ranked table. Supports `preset="heuristic" | "pseudobulk" | "permissive"`. Works for both `active_score` and `differential_expression` results.

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
- `diagnose_design(adata, groupby, target_group, reference_group, sample_col=None)` — analyzes cell/sample counts and global unspliced fraction; returns warnings, recommendations, and a suggested `filter_active_genes` preset. Automatically called internally when `sample_col` or `use_pseudobulk=True` is used.
- `run_enrichment(...)`, `run_kegg(...)`, `simplify_enrichment(...)`, `list_bundled_gene_sets()`
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

All `scat.pl.*` functions support `ax=` / `axes=` (for embedding in multi-panel figures) and `save_path=` (high-quality 300 dpi output).

### Main Plotting Functions

- `scat.pl.comet_plot(results_df, top_n=12, ...)`  
  Recommended: log fold change vs. bias-corrected unspliced residual (velocity_residual), sized and colored by active_score.

- `scat.pl.volcano_plot(results_df, top_n=10, label_genes=None, ...)`  
  2D volcano (logFC vs. -log10(p_adj)). Supports `label_genes=[...]` for manual gene labels
  (combined with top_n) — ggVolcano style flexibility. Classic up/down/ns coloring when
  not using active_score. See https://github.com/BioSenior/ggVolcano for style inspiration.

- `scat.pl.bias_diagnostic_plot(results_df, ...)`  
  Before/after view of the effect of length+intron bias correction on the velocity delta.

- `scat.pl.enrich_dotplot(enrich_df, top_n=15, show_terms=None, ...)`  
  Enrichment dot plot (clusterProfiler style). `show_terms` accepts int (top N) or
  list of term strings/Descriptions (exact or partial match, order preserved) —
  directly analogous to `dotplot(..., showCategory=...)` / `showCategory=c("...")`.
  Also available as `enrich_barplot`.

- `scat.pl.volcano_3d(results_df, ...)`  
  3D volcano (logFC × -log10(p) × velocity_residual).

- `scat.pl.active_score_rankplot(results_df, top_n=20, ...)`  
  Simple horizontal barplot of top active scores.

- `scat.pl.active_genes_heatmap(adata, genes, groupby=..., ...)`  
  Convenience wrapper around `scanpy.pl.heatmap` for selected genes.

- `scat.pl.velocity_phase_portraits(adata, genes, groupby=..., ...)`  
  Quick unspliced vs. spliced phase portraits for selected genes (useful for inspecting nascent excess).

- `scat.pl.set_style()` and `scat.pl.style_context()`  
  Control global publication-style settings (vector fonts, minimal ink, etc.).

- `scat.pl.set_nature_style()` (legacy alias for `set_style`).

---

## Command-Line Interface

Only the gene-feature generator is exposed as a CLI (`generate-gene-features`).

---

## Additional Capability: Standalone Differential Expression

While the primary focus of scATrans is composite active transcription scoring from spliced/unspliced (velocity) data via `active_score`, the package also provides a general-purpose differential expression entry point that does **not** require velocity layers.

```python
import scatrans as scat

# Early (right after load + basic QC, before HVG/normalize/log):
scat.store_raw_counts(adata, layer="counts", save_raw=True)

# Works on regular count AnnData (no spliced/unspliced needed)
adata, de_results = scat.differential_expression(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    # de_method="t-test_overestim_var",   # or "wilcoxon", etc. (default)
    # use_memento_de=True,                # optional: use the integrated Memento (Cell 2024) backend
    # memento_capture_rate=0.07,
)

# Then use the same downstream tools as with active_score results
candidates = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3)

# IMPORTANT: when calling run_enrichment / run_kegg, make sure the implicit or
# explicit `universe` / `background` is the full set of genes that were measured
# in the experiment (not just the HVGs in the current adata).
enrich = scat.run_enrichment(candidates.index.tolist(), gene_sets="GO_Biological_Process")  # defaults to built-in Hs/Mm 2026 for the organism
scat.pl.volcano_plot(de_results)
scat.pl.enrich_dotplot(enrich)
```

`differential_expression` supports the same flexible backends as `active_score` (scanpy methods, PyDESeq2 pseudobulk, mixed models, and optionally Memento as a method-of-moments estimator). The returned table is directly compatible with `filter_active_genes`, enrichment functions, and all `scat.pl.*` plotting helpers.

This makes the package useful even if you only need modern DE + enrichment + visualization, while the core `active_score` workflow remains the recommended path when you have velocity information.

See `examples/memento_de_example.py` for a complete demonstration of both the velocity-focused and pure-DE usage patterns.

**Important: raw counts requirement**

Count-based backends (Memento, PyDESeq2) work best with raw integer counts. The very common pattern of

```python
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

leaves `adata.X` as log-transformed HVGs only, which is unsuitable.

**Recommended practice** (do early):

```python
import scatrans as scat

# Before HVG + normalize + log1p
scat.ensure_raw_counts(adata)          # saves raw counts to adata.layers["counts"]

# Then normal Scanpy preprocessing
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Now safe
adata, de_res = scat.differential_expression(adata, use_memento_de=True, ...)
```

`ensure_raw_counts()` will also try to recover from `adata.raw`. The functions emit clear warnings when they detect this situation.

---

## License

MIT License.

---

*This README emphasizes the basic, honest, low-ceremony workflow centered on active transcription analysis from velocity data. Advanced capabilities (including standalone DE with Memento support) remain available for users who need them.*
