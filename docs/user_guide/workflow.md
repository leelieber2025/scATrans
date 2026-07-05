# Core Workflow

## Run `active_score` (default parameters)

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    show_plot=True,   # shows a comet plot for quick visual check
)
```

This computes differential expression, reference-group gamma excess for the
unspliced layer, optional Huber bias correction on gene length and intron
number, a composite active score, and stores diagnostics in
`adata_res.uns["scatrans"]["diagnostics"]`.

### Common basic switches: pseudobulk and DE test method

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

When using `use_pseudobulk=True` + `pseudobulk_de_backend="scanpy"`, the
`de_method` you choose (including `"wilcoxon"`) will be used for the
pseudobulk DE step.

These choices are recorded in `adata_res.uns["scatrans"]` (`de_method`,
`pseudobulk_de_backend`, `use_pseudobulk`).

The `filter_active_genes` helper has a `preset="pseudobulk"` that applies
more lenient default thresholds suitable after aggregation.

### Choosing a DE backend (decision guide)

| Your design | Recommended backend | Caveats |
|-------------|---------------------|---------|
| Exploratory / default | scanpy `wilcoxon` or `t-test` on normalized data | Fast; standard pseudoreplication limits |
| ≥2 biological replicates per group, aggregated counts | `use_pseudobulk=True` + `pydeseq2` | Requires raw counts (`store_raw_counts`); DESeq2 assumptions |
| Few pseudobulk samples, no DESeq2 | `use_pseudobulk=True` + `pseudobulk_de_backend="scanpy"` | Non-parametric on aggregated profiles |
| Cell-level data + true sample replicates | `use_mixed_model=True` + `sample_col` | Lightweight LMM (log1p); check `diagnostics["mixed_model"]["failed_fit_rate"]` — not NB-GLMM/voom |
| Method-of-moments cell-level DE | `use_memento_de=True` | Raw integer counts required; compare `memento_p_adj_native` vs package `p_adj` if auditing |

Always run `recommend_workflow(...)` first; inspect
`adata.uns["scatrans"]["diagnostics"]` (bias, gamma, permutation
`disabled_reason`) before publication claims.

## Gene filtering with `filter_active_genes` (core output tool)

The internal `significant` list is strict. Users typically filter the full
table returned in `all_results` with `filter_active_genes`.

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

# Reproduce the built-in `significant` list exactly (requires use_permutation=True upstream)
builtin_again = scat.filter_active_genes(all_results, preset="significant")
assert builtin_again.index.tolist() == significant.index.tolist()

# Advanced usage
mask = scat.filter_active_genes(all_results, return_mask=True)  # boolean Series
filtered_inplace = scat.filter_active_genes(all_results, preset="heuristic", inplace=True)
# or preset="pseudobulk" after aggregation, or preset="permissive"
```

**`preset="significant"`** (aliases: `"builtin"`, `"active_score_significant"`)
replays the built-in `significant` mask from `active_score` using metadata in
`all_results.attrs["scatrans_filter_context"]`. It requires
`use_permutation=True` on the upstream run. When permutation FDR was
disabled (e.g. too few pseudobulk shuffles), `preset="heuristic"` is often a
better exploratory fallback than `preset="significant"`.

For pure `differential_expression()` results you can also select
downregulated genes:

```python
down_cands = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="down")
both = scat.filter_active_genes(de_results, pval_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="both")
```

The helper safely ignores filters for columns that do not exist (e.g.
`unspliced_excess_fdr` when you did not use `use_permutation`). Legacy column
names `velocity_residual` / `velocity_delta_raw` remain in `adata.var` as
aliases.

## `diagnose_design`

`diagnose_design` analyzes experimental design (cell counts, replicate
numbers, global unspliced fraction) and returns warnings and
recommendations. It is called automatically inside `active_score` when
`sample_col` or `use_pseudobulk=True`.

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

print("Warnings:")
for w in diag["warnings"]:
    print("  -", w)

print("\nRecommendations:")
for r in diag["recommendations"]:
    print("  -", r)

print("\nSuggested preset for filter_active_genes:", diag.get("suggested_preset"))
```

**What it returns**, a dictionary containing:

- `n_cells_target`, `n_cells_reference`
- `n_samples_target`, `n_samples_reference` (when `sample_col` is provided)
- `unspliced_global_fraction`
- `warnings`: list of strings (e.g. low power warnings)
- `recommendations`: list of strings
- `suggested_preset`: `"heuristic"`, `"pseudobulk"`, or `None`

`diagnose_design` is automatically called inside `active_score(...)`
whenever you pass `sample_col` or set `use_pseudobulk=True`. You will see its
output in the log.

## Layer names

The package auto-detects `mature`/`nascent` (kb_python) and remaps them
internally. You can also pass `spliced_layer=...` and `unspliced_layer=...`
explicitly.
