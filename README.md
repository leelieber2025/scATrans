# scATrans

scATrans is a Python toolkit for single-cell differential analysis. It is primarily designed for datasets that contain spliced/unspliced (or mature/nascent) RNA layers. In this setting it computes a composite active transcription score that integrates differential expression with reference-based excess unspliced RNA to rank genes.

It also supports conventional differential expression workflows (no velocity data required) using scanpy, PyDESeq2 pseudobulk, linear mixed models, or optional Memento. Functional enrichment (ORA, GSEA, GO, KEGG) uses bundled gene sets with consistent universe handling, and a set of visualization functions is provided.

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

The package ships precomputed gene feature tables (gene length + intron number) for both mouse and human. These are used for optional bias correction in `active_score`. You can also supply custom tables (e.g. from your own GTF).

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

**Quick data quality check.** Before analysis, inspect the global unspliced fraction:

```python
import scatrans as scat
ufrac = scat.qc.unspliced_global(adata)   # logs INFO + WARNING if > 50%
```

`active_score` automatically runs this check and records the value in diagnostics.

## Quick Start (Complete End-to-End Example)

This is a complete, copy-paste friendly workflow for first-time users. It takes you from loaded data to differential results, filtering, enrichment, and visualization of enrichment results.

```python
import scanpy as sc
import scatrans as scat

# 1. Load your data (must contain spliced/unspliced layers or use differential_expression instead)
adata = sc.read_h5ad("your_data.h5ad")

# 2. Store raw counts + original layers early (before HVG/normalization)
scat.store_raw_counts(adata, layer="counts", save_raw=False)

# 3. Standard preprocessing (adjust as needed for your analysis)
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 4. Attach gene features for bias correction (optional)
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# 5. Run differential analysis (active transcription score)
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    show_plot=False,
)

print("Differential analysis results (top rows):")
print(all_results.head())

# 6. Gene filtering (use the full table; the built-in 'significant' is often empty)
candidates = scat.filter_active_genes(
    all_results,
    preset="heuristic",           # or "pseudobulk" / "permissive"
    # active_score_cutoff=30,
    # logfc_cutoff=0.3,
    # pval_cutoff=0.05,
)

print(f"\nFiltered candidate genes: {len(candidates)}")

# 7. Functional enrichment (GO)
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",   # or "GO_BP"
    organism="mouse",                    # or "human"
    adata=adata,                         # uses stored raw genes as background
    pval_cutoff=0.05,
)

print("\nTop GO enrichment terms:")
print(enrich_res.head())

# KEGG enrichment (alternative)
kegg_res = scat.run_kegg(
    gene_list=candidates.index.tolist(),
    organism="mouse",   # or "human"
    adata=adata,
)

# 8. Visualize enrichment results
scat.pl.enrich_dotplot(enrich_res, top_n=15, title="GO Enrichment")
scat.pl.enrich_dotplot(kegg_res, top_n=10, title="KEGG Pathways")

# Optional: save figures
# scat.pl.enrich_dotplot(enrich_res, top_n=12, save_path="enrich_go.pdf")

# Optional: main result plots
# scat.pl.comet_plot(all_results, top_n=12)
# scat.pl.volcano_plot(all_results, top_n=10)
```

You can now explore `all_results`, adjust filters in step 6, try different `gene_sets`, or run `run_go` / `run_gsea`.

For pure differential expression without spliced/unspliced layers, replace step 5 with `scat.differential_expression(...)`.

### Preserving raw counts and layers

Call `store_raw_counts` early (after loading and QC, before HVG or normalization). It writes the current `.X` to `layers["counts"]` and copies the original spliced/unspliced layers. These survive later subsetting and provide the correct background for enrichment and count-based DE.

The default `save_raw=False` avoids populating `adata.raw`.

After HVG-based visualization on a copy, restore or use the preserved layers for full-gene DE, active scoring, or enrichment (pass `adata=` to `run_enrichment` or `run_kegg` to use the stored gene list as background).

HVG subsetting also subsets the saved layers. This keeps velocity calculations consistent with `.X`. To analyze more genes than the HVG set, store before subsetting or operate on the unfiltered object for DE and enrichment steps.

To restore raw counts into `.X` for the current gene set:

```python
adata_raw = scat.restore_raw_counts(adata, layer="counts", inplace=False)
```

See the standalone differential expression section for the no-velocity use case.

---

## Core Workflow

`active_score` performs differential expression, reference-gamma unspliced excess calculation, optional bias correction, composite scoring, and stores results plus diagnostics. Downstream steps commonly include gene filtering with `filter_active_genes`, functional enrichment, and plotting.

The internal `significant` list uses strict thresholds. The complete results table is returned as `all_results`; use `filter_active_genes` for custom criteria. Diagnostics are available under `adata_res.uns["scatrans"]["diagnostics"]`.

---

## Basic Analysis Workflow

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

This computes differential expression, reference-group gamma excess for the unspliced layer, optional Huber bias correction on gene length and intron number, a composite active score, and stores diagnostics in `adata_res.uns["scatrans"]["diagnostics"]`.

### 3.1.1 Common basic switches: pseudobulk and DE test method

These are standard options available for most analyses.

**Pseudobulk mode** (use when you have multiple biological replicates per condition):

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

The internal `significant` list is strict. Users typically filter the full table returned in `all_results` with `filter_active_genes`.

```python
# Start permissive, then tighten based on your data
candidates = scat.filter_active_genes(
    all_results,
    active_score_cutoff=30,
    unspliced_excess_residual_cutoff=0.5,
    unspliced_excess_fdr_cutoff=0.05,
    logfc_cutoff=0.3,
    pval_cutoff=0.05,
)

# Or use presets that choose reasonable defaults for common analysis styles
candidates = scat.filter_active_genes(all_results, preset="heuristic")

# Advanced usage
mask = scat.filter_active_genes(all_results, return_mask=True)  # boolean Series
filtered_inplace = scat.filter_active_genes(all_results, preset="heuristic", inplace=True)
# or preset="pseudobulk" after aggregation, or preset="permissive"
```

The helper safely ignores filters for columns that do not exist (e.g. `unspliced_excess_fdr` when you did not use `use_permutation`). Legacy column names `velocity_residual` / `velocity_delta_raw` remain in `adata.var` as aliases.

### 3.3 Functional enrichment

Over-representation analysis is available via `run_enrichment`:

```python
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",   # or "GO_BP" — automatically resolved to the
                                         # correct organism-specific built-in (Hs/Mm_GO_..._2026)
    organism="mouse",  # or "human"
    adata=adata,   # if you called store_raw_counts(adata) earlier, this will
                   # automatically use the preserved full measured gene list as universe.
                   # Explicit `universe=` still takes precedence.
    pval_cutoff=0.05,
    min_size=5,
    max_size=500,
)
# Additional columns and attrs (clusterProfiler compatibility):
#   - "neg_log10_padj" column
#   - res.attrs["universe_info"] with effective_universe_size, dropped_by_annotation_filter, etc.
```

**run_gsea** (pre-ranked GSEA for ranked gene lists):

```python
# ranked list from active_score / differential_expression results
ranked = all_results["logFC"].sort_values(ascending=False)

gsea_res = scat.run_gsea(
    ranked_genes=ranked,
    gene_sets="GO_Biological_Process",
    organism="mouse",  # or "human"
    nperm=1000,
)
print(gsea_res.head())

scat.pl.enrich_dotplot(gsea_res, x="NES", color_by="NES")

# gseaplot (uses the exact curve stored by run_gsea)
term = gsea_res.iloc[0]["Term"]
scat.pl.gseaplot(ranked, gsea_res, term=term)
```

Requires `pip install "scatrans[gsea]"`.

**run_kegg** (convenience wrapper for KEGG pathways):

```python
kegg_res = scat.run_kegg(
    gene_list=candidates.index.tolist(),
    organism="mouse",  # or "human"
    # Defaults to the organism-specific built-in library (Hs_KEGG_2026 or Mm_KEGG_2026)
    adata=adata,   # if store_raw_counts was called earlier, this automatically uses
                   # the preserved full measured gene set as background.
    pval_cutoff=0.05,
)
```

### Using bundled gene sets

The package defaults to organism-specific bundled sets. Use `organism=` together with base names such as `"GO_Biological_Process"` or `"KEGG"`. Supply a full historical name (e.g. `"GO_Biological_Process_2023"`) to select an Enrichr version.

```python
# KEGG example
kegg = scat.run_kegg(gene_list=genes, organism="mouse")  # or "human"

# GO — base name is enough (automatically resolved to Hs/Mm_GO_..._2026)
go = scat.run_enrichment(
    gene_list=genes,
    gene_sets="GO_Biological_Process",   # or "GO_BP"
    organism="mouse",  # or "human"
    # pass adata= (after you did store_raw_counts early) to use the preserved
    # full measured genes as universe/background automatically.
    # Explicit universe= or background= always takes precedence.
    adata=adata,
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
    # For background: still prefer adata= (from store_raw_counts) over manual universe=.
    adata=adata,
    # universe=background,   # explicit still accepted and takes precedence
)

# Supply the full name containing the year to select a historical version.
```

`gene_set_source` can be used as an explicit override when needed:

- `gene_set_source="scatrans"` → use bundled sets
- `gene_set_source="enrichr"` → use gseapy/Enrichr libraries

List available bundled sets:

```python
print(scat.list_bundled_gene_sets())
```

**Adding your own sets**: Drop `.gmt` files into `src/scatrans/data/`. See `src/scatrans/data/README.md`.

**simplify_enrichment** (reduce redundant enrichment terms):

Two methods are supported:

- **`jaccard`** (default): greedy filtering by Jaccard overlap of enriched gene lists.
- **`pathway_denester`**: combinatorial nested-pathway test adapted from [PathwayDenester](https://github.com/Helmy-Lab/PathwayDenester). Better at removing terms that are significant only because they are nested inside a more significant parent pathway. Requires full pathway gene memberships (auto-loaded from `enrich_res.attrs` when enrichment used bundled/Enrichr libraries; pass `gene_sets=` again if you used a custom dict).

```python
# Jaccard (fast, overlap-based)
simplified = scat.simplify_enrichment(
    enrich_res,
    similarity_cutoff=0.5,
    min_count=3,
    method="jaccard",
)

# PathwayDenester (nested-pathway test)
simplified = scat.simplify_enrichment(
    enrich_res,
    method="pathway_denester",
    min_count=3,
    pval_threshold=0.05,       # independence cutoff
    to_test_threshold=0.0,     # min shared-DEG fraction before testing
    term_size_limit=0,         # e.g. 500 to drop very broad terms
    show_excluded=False,       # True keeps excluded terms + Denester_* diagnostics
)
```

`run_kegg` and `simplify_enrichment` are wrappers around `run_enrichment`.

### run_go (GO enrichment, clusterProfiler-style)

```python
# Biological Process (defaults to the bundled Mm/Hs_GO_Biological_Process_2026)
go_bp = scat.run_go(
    gene_list=markers,
    ontology="BP",          # "BP", "CC", "MF", or "ALL"
    organism="mouse",  # or "human"
    adata=adata,
    return_all=True,
)

# ALL three ontologies + unified multiple-testing correction across them
go_all = scat.run_go(
    markers, ontology="ALL", organism="mouse",  # or "human"
    return_all=True,
    adjust_across_all=True,
)
# go_all.attrs["per_ontology_attrs"] contains full diagnostics for BP/CC/MF separately
```

`run_go` automatically resolves to the organism-specific bundled sets when possible (BP is bundled; CC/MF fall back to gseapy/Enrichr if the library is installed).

### Exporting results

The following helpers export results:

```python
res = scat.run_kegg(genes, organism="mouse", return_all=True, include_gene_list=True)  # or "human"

saved = scat.save_enrichment_report(
    res,
    prefix="cluster1_kegg",   # or "results/suppl/my_enrich" (directories created automatically)
    save_excel=True,
    save_csv=True,
    save_tsv=True,            # often preferred for gene symbols + Excel locale safety
    save_metadata=True,
    save_term_gene_table=True,
)

# saved -> {'results_csv': ..., 'results_tsv': ..., 'term_gene_table_csv': ..., 'metadata_json': ..., 'results_xlsx': ...}

# Long-format term–gene table (one row per gene; perfect for networks, follow-up stats, etc.)
long_table = scat.expand_enrichment_genes(res)
# If the input was from run_go(ontology="ALL"), long_table will have an "Ontology" column first.
```

`save_enrichment_report` also writes a rich `metadata.json` (and a "metadata" sheet in the xlsx) containing:
- `analysis_info` (package, version, timestamp)
- `gene_set_info` (requested/resolved + `requested_source` vs `actual_source`: "bundled", "gseapy", "gmt", "dict")
- `universe_info` (effective N, dropped genes, restrict behavior, etc.)
- Full `.attrs` from the enrichment call (including per-ontology details for GO ALL)

All empty results still carry diagnostic `.attrs` (`reason`, `gene_set_info`, `universe_info`, etc.) so you never lose information when a call returns no terms.

### Additional enrichment plot options

For basic usage see the Quick Start example above. Additional controls:

```python
import scatrans as scat

# Dot plot for ORA results from run_enrichment / run_kegg / run_go
# x-axis defaults to "GeneRatio"; other options: "Count", "FoldEnrichment", "-log10(p.adj)"
scat.pl.enrich_dotplot(
    enrich_res,
    top_n=15,
    title="GO Biological Process enrichment",
)

# For KEGG
scat.pl.enrich_dotplot(
    kegg_res,
    top_n=10,
    title="KEGG pathways",
)

# Save figure (vector-friendly, 300 dpi)
scat.pl.enrich_dotplot(
    enrich_res,
    top_n=12,
    save_path="go_enrichment.pdf",
)

# GSEA results (auto-switches to NES on x and color when NES column present)
scat.pl.enrich_dotplot(
    gsea_res,
    top_n=15,
    x="NES",
    color_by="NES",
)

# GSEA running-sum plot (uses curves stored by run_gsea)
if len(gsea_res) > 0:
    term = gsea_res.iloc[0]["Term"]
    scat.pl.gseaplot(
        ranked_genes=ranked,
        gsea_result=gsea_res,
        term=term,
        save_path="gsea_running_score.pdf",
    )

# Embed in multi-panel figure with ax=
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 5))
scat.pl.enrich_dotplot(enrich_res, top_n=8, ax=ax, show=False)
fig.savefig("enrich_panel.pdf", dpi=300, bbox_inches="tight")
```

Additional options:
- `show_terms=15` or `show_terms="auto"` or `show_terms=["term A", "term B"]`
- `use_style=True` to apply publication style for that call only

All `scat.pl.*` functions accept `save_path`, `ax`, `figsize`, `show`, and `dpi`.

### 3.4 Visualization

```python
import scatrans as scat

scat.pl.comet_plot(all_results, top_n=12, title="Active Drivers")
scat.pl.volcano_plot(all_results, top_n=10, label_genes=["YourGene1", "YourGene2"])
scat.pl.bias_diagnostic_plot(all_results)
```

All plotting functions support `ax=` / `axes=` for multi-panel figures and `save_path=` (300 dpi output).

---

## Helper Functions

### diagnose_design

`diagnose_design` analyzes experimental design (cell counts, replicate numbers, global unspliced fraction) and returns warnings and recommendations. It is called automatically inside `active_score` when `sample_col` or `use_pseudobulk=True`.

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
    sample_col="sample"          # required for pseudobulk and mixed-model paths when replicates exist
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

Enrichment functions are covered in section 3.3. Run KEGG pathway enrichment (see 3.3):

```python
kegg_res = scat.run_kegg(
    gene_list=significant.index.tolist(),   # or from all_results
    organism="mouse",  # or "human"
    adata=adata,   # preferred: if you called store_raw_counts earlier, this will
                   # automatically use the preserved full measured gene set.
    pval_cutoff=0.05,
    min_size=5,
    max_size=500,
    return_all=False,
    # Defaults to the organism-specific built-in (Hs/Mm_KEGG_2026).
    # To use a historical Enrichr version: kegg_library="KEGG_2021"
)

### run_gsea (pre-ranked GSEA)

For ranked-list enrichment (the classic GSEA / prerank approach):

```python
# Obtain a ranked gene list (e.g. from active_score or differential_expression results).
# Convention: higher value = more associated with the target group.
ranked = all_results["logFC"].sort_values(ascending=False)

gsea_res = scat.run_gsea(
    ranked_genes=ranked,
    gene_sets="GO_Biological_Process",
    organism="mouse",  # or "human"
    nperm=1000,
    min_size=15,
)
print(gsea_res[["Term", "NES", "p.adjust", "leading_edge"]].head())

# Works with existing plotting helpers (auto-detects GSEA columns)
scat.pl.enrich_dotplot(gsea_res, x="NES", color_by="NES")

# Dedicated running-sum plot (uses curves stored by run_gsea)
scat.pl.gseaplot(ranked, gsea_res, term=gsea_res.iloc[0]["Term"])
```

`run_gsea` stores the full enrichment score curves in `.attrs["gsea_details"]` so that `gseaplot` renders exactly the same RES that produced the reported NES/p-values.

Requires the optional extra:
```bash
pip install "scatrans[gsea]"   # pulls in gseapy
```

**simplify_enrichment** (see full documentation and examples in 3.3):

```python
# Jaccard: drop terms whose enriched gene sets overlap strongly with a kept term
simplified = scat.simplify_enrichment(
    kegg_res,
    similarity_cutoff=0.5,
    min_count=3,
    by="p.adjust",
    method="jaccard",
)

# PathwayDenester: drop nested pathways explained by a more significant parent
simplified = scat.simplify_enrichment(
    kegg_res,
    method="pathway_denester",
    min_count=3,
    by="p.adjust",
    gene_sets="KEGG",            # optional if kegg_res.attrs records the library
    pval_threshold=0.05,
    to_test_threshold=0.0,
)

See 3.3 for usage examples and parameters.

---

## Result Interpretation

### Column naming (v0.9+)

Primary result columns use **unspliced / nascent excess** terminology (not RNA velocity):

| Primary column | Legacy alias (deprecated) | Meaning |
|----------------|---------------------------|---------|
| `unspliced_excess_delta` | `velocity_delta_raw` | Raw U − γ_ref·S in target group |
| `unspliced_excess_residual` | `velocity_residual` | Bias-corrected excess residual |
| `unspliced_excess_pval` | — | One-sided permutation p-value on residual |
| `unspliced_excess_fdr` | — | BH-FDR on `unspliced_excess_pval` |

`active_score` (0–100) is a **heuristic ranking score** (weighted soft-scaled composite of logFC + unspliced excess residual + -log p_adj). It is intended **for ranking and visualization only** and should **not** be interpreted or reported as a p-value or statistical significance measure. Use the permutation-derived `unspliced_excess_fdr` (when enabled) or your own post-hoc statistics for claims.

### Built-in `significant` gene list

When `use_permutation=True`, the internal mask requires **all** of:

- `logFC > logfc_cutoff` (default 0.5)
- `p_adj < pval_cutoff` (default 0.05)
- `unspliced_excess_residual > 0`
- `unspliced_excess_fdr < unspliced_excess_fdr_cutoff` (default 0.05)

Without `use_permutation=True`, the built-in `significant` list is **empty** (FDR on unspliced excess cannot be computed). Use `all_results` + `filter_active_genes` for custom thresholds.

On real data the built-in list often returns zero or few genes. Use the full table in `all_results`, sorted by `active_score` descending.

After each run inspect the diagnostics:

```python
meta = adata_res.uns["scatrans"]
print(meta["diagnostics"]["unspliced_global_fraction"])
print(meta["diagnostics"]["bias_correction"])
print(meta.get("permutation_approximation_note"))
```

Global unspliced fractions above ~50% frequently indicate technical issues. Bias-correction diagnostics report the number of genes used and any fallback behavior. The permutation note records that unspliced/spliced layers and the reference gamma were fixed for speed while labels were shuffled.

---

## Optional Advanced Features

The following flags are disabled by default and should be enabled only when required by the experimental design:

- `use_permutation=True`
- `bias_correction="none"`
- `show_effective_gamma=True`
- `gamma_method="robust_median"` (or "raw")
- `use_mixed_model=True`
- `prioritize_velocity=True`

`diagnose_design` summarizes cell and sample counts plus global unspliced fraction and returns warnings and a suggested `filter_active_genes` preset. It runs automatically when `sample_col` or `use_pseudobulk=True` is supplied.

Inspect the corresponding diagnostics after enabling any advanced option.

### use_permutation=True

**Required for the built-in `significant` list** (via `unspliced_excess_fdr`).

Adds:

- `unspliced_excess_pval` / `unspliced_excess_fdr` — permutation significance on the bias-corrected unspliced excess residual (one-sided, positive direction). **Use these for active-gene calls.**
- `active_score_pval` / `active_score_fdr` — permutation on the composite heuristic score (ranking aid only).

The permutation shuffles only group labels; unspliced/spliced layers and the reference gamma are fixed from the original labeling for speed. **This is a conditional permutation** (conditioned on the observed velocity structure and gamma). It is a speed/tractability tradeoff and **not an unconditional permutation of the full data**. In small reference groups or strong batch effects, interpret the resulting FDR with extra caution; always inspect diagnostics and consider biological replicates.

See diagnostics["velocity"] for the actual gamma_method and prior_weight used.

```python
adata_res, significant, all_results = scat.active_score(
    adata,
    use_permutation=True,
    n_perm=500,
    unspliced_excess_fdr_cutoff=0.05,
)
```

### bias_correction

By default the package applies a Huber regression of the raw unspliced excess delta on log(gene length) and log(intron number) and uses the residuals as `unspliced_excess_residual`. This step can be disabled by setting `bias_correction="none"`, in which case the raw (reference-gamma corrected) delta is used directly.

The correction is intended to reduce technical contributions from gene length and intron number to the unspliced excess term. Whether length or intron number carry biological signal of interest in a given dataset is a scientific judgment that the user must make; the correction is therefore optional. The `bias_diagnostic_plot` function can be used to inspect the relationship before and after correction.

### gamma_method and reference gamma robustness

The core unspliced excess uses a per-gene reference gamma = U_ref / S_ref (shrunk).

- Default: `gamma_method="heuristic_shrink"` + `prior_weight=5.0` (additive pseudo-count shrinkage toward a global ratio).
- For small reference groups, try `gamma_method="robust_median"`: uses the **median** ratio across reference genes as the anchor. This reduces sensitivity to a few outlier genes in the reference and can yield more stable residuals.
- `gamma_method="raw"` disables most shrinkage (exploratory only).

The chosen method, prior_weight, and summary stats of the realized effective_gamma are **always** written to diagnostics:

```python
v = adata_res.uns["scatrans"]["diagnostics"]["velocity"]
print(v["gamma_method"], v["prior_weight"], v["effective_gamma_stats"])
```

Shrinkage strength and stability are now visible without `show_effective_gamma`.

### show_effective_gamma=True

Adds the column `effective_gamma` (reference-group shrunk U/S ratio) to `adata.var` and to the results tables. Many genes will have similar values in pure heuristic mode; advanced (moments) mode usually shows more per-gene variation.

Example filter using the column (when present):

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

**Practical note on small numbers of samples:** With very few biological replicates, pseudobulk aggregation can drive most `unspliced_excess_residual` values close to zero. In such regimes the cell-level mixed-model path (`use_mixed_model=True`, `use_pseudobulk=False`) often preserves more of the nascent-excess signal while still respecting sample structure.

The mixed-model settings and median `delta_variance` are recorded in diagnostics.

### mode="advanced"

Uses scVelo moments for local smoothing before computing the group-wise gamma delta. It is still a simple reference-gamma excess calculation on the smoothed moments, not a full stochastic or dynamical model.

Use when you have sufficient cells and want local smoothing. The function falls back to heuristic mode on failure (`advanced_fallback=True` by default).

---

## API Reference (Simplified)

### Core functions

- `active_score(...)` — main analysis for active transcription from velocity data. Returns `(adata_res, significant, all_results)`.
- `differential_expression(...)` — standalone DE (no velocity data required). Supports the same backends as `active_score` (including optional Memento). Returns `(adata, results_df)`.
- `filter_active_genes(results_df, ...)` — post-filter the full ranked table. Supports `preset="heuristic" | "pseudobulk" | "permissive"`. Works for both `active_score` and `differential_expression` results.
- `store_raw_counts` / `ensure_raw_counts` / `restore_raw_counts` — preserve the original full count matrix and spliced/unspliced layers (call early, before HVG/normalize). Essential for correct backgrounds in enrichment and for count-based DE backends.

### Common parameters

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

- `use_permutation`, `n_perm`, `unspliced_excess_fdr_cutoff` (and deprecated `active_fdr_cutoff`)
- `bias_correction` ("huber_length_intron" or "none")
- `show_effective_gamma`, `gamma_method`, `prior_weight`
- `use_mixed_model`, `use_delta_variance_pval`, `mixed_model_pval`
- `mode` ("heuristic" or "advanced")

Full signatures and all parameters are documented in the function docstrings and the source.

### Other commonly used functions

- `add_gene_features(adata, organism=..., ...)` — attach length/intron info (`"mouse"` or `"human"`)
- `generate_gene_features_from_gtf(gtf_path, output_name, ...)` — build a custom table from a GTF (requires `[gene_features]`)
- `list_available_gene_features()` — list bundled tables
- `generate-gene-features` (CLI) — same as above, for the shell
- `store_raw_counts(adata, layer="counts", save_raw=False)`, `ensure_raw_counts(adata)`, `restore_raw_counts(adata, ...)` — preserve full raw counts + original spliced/unspliced layers before HVG/normalization (critical for correct DE, enrichment background, and Memento/PyDESeq2)
- `diagnose_design(adata, groupby, target_group, reference_group, sample_col=None)` — analyzes cell/sample counts and global unspliced fraction; returns warnings, recommendations, and a suggested `filter_active_genes` preset. Automatically called internally when `sample_col` or `use_pseudobulk=True` is used.
- `run_enrichment(...)`, `run_kegg(...)`, `run_go(...)`, `run_gsea(...)`, `simplify_enrichment(...)`, `save_enrichment_report(...)`, `expand_enrichment_genes(...)`, `list_bundled_gene_sets()`
- `scat.pl.*` plotting functions (comet_plot, volcano_plot, bias_diagnostic_plot, enrich_dotplot, gseaplot, active_score_rankplot, active_genes_heatmap, velocity_phase_portraits, ...)
- `scat.qc.unspliced_global(adata)`
- `scat.pl.set_style()` / `scat.pl.style_context()` — publication-friendly matplotlib style (opt-in, off by default per-plot via use_style=)
- Submodules `scat.pl` and `scat.qc` (scanpy-style)

### Layer names

The package auto-detects `mature`/`nascent` (kb_python) and remaps them internally. You can also pass `spliced_layer=...` and `unspliced_layer=...` explicitly.

---

## Gene Feature Attachment & CLI

Gene length and intron count are used for optional bias correction inside `active_score`.

```python
# Use bundled tables
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# or provide your own table
adata = scat.add_gene_features(adata, gene_features_path="my_features.parquet")
```

The package includes tables for mouse and human. Use `organism="mouse"` (default) or `organism="human"` when calling `add_gene_features`. For other species or custom annotations use the gene feature generator CLI.

### Generating a custom table from GTF

Install the generator:

```bash
pip install "scatrans[gene_features]"
```

Use the CLI (works with 10x `genes.gtf` or GENCODE GTFs):

```bash
# Mouse
generate-gene-features --gtf /path/to/genes.gtf \
                       --output my_mouse_features.parquet \
                       --organism mouse

# Human (GENCODE or 10x)
generate-gene-features --gtf gencode.v49.primary_assembly.annotation.gtf \
                       --output human_GRCh38_2024A_gene_features.parquet \
                       --organism human
```

Then use it:

```python
import scatrans as scat

adata = scat.add_gene_features(
    adata,
    gene_features_path="human_GRCh38_2024A_gene_features.parquet"
)

# bias correction will now be able to use length + intron_number
adata_res, significant, all_results = scat.active_score(adata, ...)
```

You can also call the generator programmatically:

```python
from scatrans import generate_gene_features_from_gtf

df = generate_gene_features_from_gtf(
    "path/to/genes.gtf",
    output_name="my_custom_features.parquet",
    organism="human"
)
```

See also `scat.list_available_gene_features()` (for bundled tables) and the full signature of `add_gene_features`.

**Tip**: The generated parquet must contain a `gene_name` column (plus `gene_length` and `intron_number`). `add_gene_features` does a `reindex` on your `adata.var_names`.

---

## Plotting Style

```python
import scatrans as scat
scat.pl.set_style()                 # once early (opt-in)
# or (to limit scope):
with scat.pl.style_context(linewidth=0.8):
    scat.pl.comet_plot(...)         # inside block or pass use_style=True
# Default for pl.* functions is use_style=False (prevents surprising rcParams changes in notebooks).
```

All `scat.pl.*` functions support `ax=` / `axes=` (for embedding in multi-panel figures), `save_path=`, `show=`, `use_style=`, `figsize=` for consistency.
Most return `(fig, ax)` (or `(fig, axes_list)` for grids like phase portraits).

### Main Plotting Functions

- `scat.pl.comet_plot(results_df, top_n=12, point_scale=1.0, min_size=2, max_size=180, s=None, ...)`  
  Plots log fold change vs. bias-corrected unspliced excess residual (`unspliced_excess_residual`), sized and colored by `active_score`.
  - `s=3` (or 1-5): force **fixed** small point size for everything (direct, simple control).
  - `point_scale=0.2` + `min_size=1`: for variable sizing, make tiniest background points truly small.

- `scat.pl.volcano_plot(results_df, top_n=10, label_genes=None, point_scale=1.0, min_size=2, s=None, ...)`  
  2D volcano (logFC vs. -log10(p_adj)). Supports `label_genes=[...]` for manual gene labels
  (combined with top_n). Classic up/down/ns coloring when not using active_score.
  Use `s=2` for uniformly small points, or min_size + point_scale for score/p-value sized tiny backgrounds.
  Especially helpful for pure DE results (no active_score).

- `scat.pl.bias_diagnostic_plot(results_df, point_size=10, ...)`  
  Before/after view of the effect of length+intron bias correction on the velocity delta.
  `point_size` controls the gene cloud density. 

- `scat.pl.volcano_3d(results_df, point_scale=..., min_size=2, s=None, ...)`  
  3D version of the volcano. Same size controls (`s` for fixed size).

- `scat.pl.enrich_dotplot(enrich_df, ...)` now also works well with GSEA results (auto defaults to `x="NES"`, diverging cmap for `color_by="NES"`).
- `scat.pl.gseaplot(ranked_genes, gsea_result, term=...)` — classic GSEA running-sum plot (uses precomputed curves from `run_gsea` when available).
- `scat.pl.enrich_dotplot(enrich_df, top_n=15, show_terms=None, x="GeneRatio", size_by="Count", color_by="Adjusted P-value", ...)`  
  Enrichment dot plot (clusterProfiler style). 
  - `x`: x-axis variable — "GeneRatio" (default for ORA), "FoldEnrichment", **"Count"**, "-log10(p.adj)", or "NES" (for GSEA).
  - `size_by` (dot size, default "Count"), `color_by` (default adjusted p-value; "NES" for GSEA uses diverging colormap).
  - `show_terms` accepts int (top N), "auto" (p.adjust <0.05 + Count>=2 smart selection), or list of term strings/Descriptions (exact or partial match, order preserved) —
    directly analogous to `dotplot(..., showCategory=...)`.
  Also available as `enrich_barplot`.

- `scat.pl.volcano_3d(results_df, ...)`  
  3D volcano (logFC × -log10(p) × unspliced_excess_residual).

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

The only console script is the gene-feature table generator:

```bash
pip install "scatrans[gene_features]"
generate-gene-features --gtf /path/to/genes.gtf --output my_features.parquet --organism human
```

See the "Gene Feature Attachment & CLI" section for full examples (mouse/human + how to use the output with `add_gene_features`).

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

# After scat.store_raw_counts(adata) early in the workflow,
# just pass adata= here. It auto-supplies the full measured gene list as background/universe.
enrich = scat.run_enrichment(
    candidates.index.tolist(),
    gene_sets="GO_Biological_Process",  # auto → correct Hs/Mm 2026 bundled
    adata=adata,
)
scat.pl.volcano_plot(de_results)
scat.pl.enrich_dotplot(enrich)
```

`differential_expression` supports the same flexible backends as `active_score` (scanpy methods, PyDESeq2 pseudobulk, mixed models, and optionally Memento as a method-of-moments estimator). The returned table is directly compatible with `filter_active_genes`, enrichment functions, and all `scat.pl.*` plotting helpers.

The package therefore supports both velocity-based active transcription analysis and conventional DE + enrichment workflows. See `examples/memento_de_example.py` for a complete demonstration of the pure-DE path.

**Important: raw counts requirement**

Count-based backends (Memento, PyDESeq2) expect raw integer counts. A common pattern that leaves unsuitable data is:

```python
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

leaves `adata.X` as log-transformed HVGs only, which is unsuitable.

**Early in the workflow:**

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

## Limitations

The unspliced excess term (used by the primary `active_score` workflow) is a group-contrast proxy derived from a reference-group gamma calculation. It is not a full stochastic or dynamical model.

The unspliced excess term is most directly applicable to binary group contrasts. Within-group heterogeneity can reduce observed signal. When `use_permutation=True`, labels are shuffled while unspliced/spliced layers and the reference gamma remain fixed; this is noted in the results. Global unspliced fractions above ~50% are reported in diagnostics. Bias correction effectiveness depends on annotation coverage. Small replicate numbers limit power for the unspliced excess term and FDR estimates. Mixed-model results tend to be conservative with large between-sample variation.

When used purely as a differential expression + enrichment toolkit (via `differential_expression`, `run_enrichment`, etc.), scATrans relies on established backends (scanpy, PyDESeq2, etc.) whose standard statistical caveats apply.

Always examine diagnostics, score distributions, and (when available) the original spliced/unspliced counts before biological interpretation.

---

## License

Apache License 2.0.

---


