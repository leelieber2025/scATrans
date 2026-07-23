# Core Workflow

## Primary entry point

{func}`~scatrans.partition_de_by_mechanism` is the recommended entry point:

1. DE selects changed genes  
2. residual-based annotation assigns soft transcription-driven versus stabilization-driven labels  
3. reliability pre-flight scales confidence (and may suppress hard gene labels at extremes)  
4. optional program-level table when `gene_sets=` is supplied  
5. optional induction-matched program tests when `induction_matched=True`  

It composes lower-level pieces documented below (`active_score`,
`filter_active_genes(select_by="de")`, `annotate_mechanism_class` /
`program_mechanism` / `program_mechanism_induction_matched`).

```python
res = scat.partition_de_by_mechanism(
    adata, groupby="condition", target_group="Disease", reference_group="Control",
    organism="mouse",
    de="builtin",            # or a DE method name / precomputed DE table / callable
    sample_col="sample",     # preferred when biological replicates exist
    # add_nascent_score=True,  # optional detection columns (not used for mechanism)
    gene_sets=my_pathways,   # optional -> program-level mechanism table
    # induction_matched=True,  # also run induction-controlled program tests
)
res.regime          # reliability pre-flight
res.selected        # DE-selected genes + soft per-gene mechanism annotation
res.programs        # program-level table (if gene_sets=)
res.programs_induction_matched  # if induction_matched=True
res.summary()       # program-first overview
res.meta.get("nascent_score")  # if add_nascent_score=True: enabled / status / …
```

### Detection columns: `add_nascent_score`

Mechanism labels **always** use the induction-normalized residual
(`unspliced_excess_residual` / abnorm residual when present). Separately, you can
opt in to **active-transcription detection** columns:

| Parameter | Effect |
|-----------|--------|
| `add_nascent_score=False` (default) | No extra columns |
| `add_nascent_score=True` | Append output of {func}`~scatrans.nascent_activity_score` to the gene table |

| Column | Meaning |
|--------|---------|
| `nascent_poisson_z` | Pseudobulk variance-stabilized nascent increase (length-robust **detection** score) |
| `dlog_unspliced` / `dlog_spliced` | CPM log fold-changes (diagnostic) |
| `de_reproducible` / `de_repro_frac` | Spliced-side DE-reproducibility flag (annotation only — never gates membership) |

The Poisson-z is **not** fed into `annotate_mechanism_class` / `program_mechanism`
(it is induction-coupled and would collapse the stabilization signal). Fail-soft:
errors land in `meta["nascent_score"]` without breaking the residual mechanism path.

Standalone scoring:

```python
nz = scat.nascent_activity_score(
    adata, groupby="condition", target_group="Disease", reference_group="Control",
    sample_col="sample",  # optional; missing name raises (no silent random folds)
)
# columns: nascent_poisson_z, dlog_unspliced, dlog_spliced, de_reproducible, de_repro_frac
```

Layers auto-resolve (`spliced`/`unspliced` or kb_python `mature`/`nascent`).

The rest of this page documents lower-level building blocks. Deprecated composite
ranking and pure-DE alternatives: {doc}`../faq`.

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
number (the bias-corrected `unspliced_excess_residual`), and stores diagnostics
in `adata_res.uns["scatrans"]["diagnostics"]`.

### Pseudobulk and DE method

**Pseudobulk mode** (multiple biological replicates per condition):

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    use_pseudobulk=True,
    sample_col="sample",                    # column identifying biological samples/individuals
    pseudobulk_de_backend="pydeseq2",       # or "scanpy"
    min_cells=5,                            # explicit override (code default is 10)
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
| Cell-level data + true sample replicates | `use_mixed_model=True` + `sample_col` | Lightweight LMM (log1p); `logFC` = sample-mean-of-means (not LMM coef); check `diagnostics["mixed_model"]` (`failed_fit_rate`, sign-discordant count) — not NB-GLMM/voom. **Do not combine with** `use_memento_de` |
| Method-of-moments cell-level DE | `use_memento_de=True` | Raw integer counts required; mutually exclusive with MixedLM; compare `memento_p_adj_native` vs package `p_adj` if auditing |

Always run `recommend_workflow(...)` first; inspect
`adata.uns["scatrans"]["diagnostics"]` (bias, gamma, permutation
`disabled_reason`) before publication claims.

## Gene filtering with `filter_active_genes`

The built-in `significant` list is strict and often empty on modestly powered
designs. Filter the full `all_results` table with `filter_active_genes`:

```python
# Explicit cutoffs (tighten for your design)
candidates = scat.filter_active_genes(
    all_results,
    unspliced_excess_residual_cutoff=0.5,
    unspliced_excess_fdr_cutoff=0.05,
    logfc_cutoff=0.3,
    padj_cutoff=0.05,  # preferred over legacy pval_cutoff=
)

# Presets for common analysis styles
candidates = scat.filter_active_genes(all_results, preset="heuristic")

# DE-only membership (padj/logFC defaults padj<0.05 & |log2FC|>1 when cutoffs
# are omitted). Nascent columns remain annotations and do not gate the list.
de_list = scat.filter_active_genes(all_results, select_by="de")

# Replay the built-in significant mask (requires use_permutation=True upstream)
builtin_again = scat.filter_active_genes(all_results, preset="significant")
assert builtin_again.index.tolist() == significant.index.tolist()

mask = scat.filter_active_genes(all_results, return_mask=True)
filtered_inplace = scat.filter_active_genes(
    all_results, preset="heuristic", inplace=True
)
# Also: preset="pseudobulk" after aggregation, or preset="permissive"
```

**`select_by="composite"` (default)** vs **`select_by="de"`**

| Mode | Who decides membership | Proxy gates | Sort |
|------|------------------------|-------------|------|
| `"composite"` | DE **and** nascent-proxy cutoffs | Applied | `p_adj` then `logFC` |
| `"de"` | DE only (`p_adj`, `logFC`, optional MixedLM coef direction) | **Skipped** (columns remain) | `p_adj` then `logFC` |

`select_by="de"` is incompatible with `preset="significant"`. The same flag is
accepted by `run_default_pipeline(..., select_by="de")` and recorded in
`meta["select_by"]`.

**`preset="significant"`** (alias: `"builtin"`)
replays the built-in `significant` mask from `active_score` using metadata in
`all_results.attrs["scatrans_filter_context"]`. It requires
`use_permutation=True` on the upstream run. When permutation FDR was
disabled (e.g. too few pseudobulk shuffles), `preset="heuristic"` is often a
better exploratory fallback than `preset="significant"`.

For pure `differential_expression()` results you can also select
downregulated genes:

```python
down_cands = scat.filter_active_genes(de_results, padj_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="down")
both = scat.filter_active_genes(de_results, padj_cutoff=0.05, logfc_cutoff=0.3, logfc_direction="both")
# or the DE-only defaults: select_by="de"
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

Returns a dictionary with:

- `n_cells_target`, `n_cells_reference`
- `n_samples_target`, `n_samples_reference` (when `sample_col` is provided)
- `unspliced_global_fraction`
- `warnings`: list of strings (e.g., low-power warnings)
- `recommendations`: list of strings
- `suggested_preset`: `"heuristic"`, `"pseudobulk"`, or `None`

`active_score` calls `diagnose_design` automatically when `sample_col` or
`use_pseudobulk=True` is set; warnings appear in the log and under
`adata.uns["scatrans"]["diagnostics"]["design"]`.

## Input data and layers

What you need depends on the path you take:

| Goal | Required in AnnData | Entry point |
|------|---------------------|-------------|
| Mechanism partition (transcription vs stabilization) | Spliced and unspliced counts as layers, **or** kb-python `mature` / `nascent` | {func}`~scatrans.partition_de_by_mechanism` |
| DE + enrichment only | Count matrix in `.X` or a counts layer (no velocity layers) | {func}`~scatrans.differential_expression` — {doc}`standalone_de` |
| Residual / diagnostics only | Same layers as mechanism | {func}`~scatrans.active_score` / {func}`~scatrans.active_score_simple` |

### Layer names

By default scATrans looks for `spliced` and `unspliced`. Layers named
`mature` / `nascent` (kb-python) are remapped internally. Override with
`spliced_layer=` and `unspliced_layer=` when your names differ.

### Raw counts and preprocessing

Pseudobulk PyDESeq2 and several DE backends need **integer raw counts**. Call
this early, before HVG filtering or log-normalization overwrites the matrix:

```python
scat.store_raw_counts(adata, layer="counts")  # or mode="auto"
```

Gene length / intron tables for optional bias correction:

```python
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"
```

Group labels must live in `adata.obs` (`groupby`, `target_group`,
`reference_group`). For pseudobulk, also set `sample_col` to a biological
replicate / individual column.

### Regime pre-flight (velocity layers)

When spliced/unspliced layers are present, run a cheap capture-quality check
before interpreting residual or mechanism annotations:

```python
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
```

`partition_de_by_mechanism` always runs this check (`result.regime`).
`run_default_pipeline` writes the same block to `result.meta["regime"]`
(fail-soft). See {doc}`advanced` for how reliability scales
`mechanism_confidence`.

Runnable end-to-end demo (`select_by="de"`, `annotate_mechanism`,
`threshold_sensitivity`, `program_mechanism`, `regime_diagnosis`) on synthetic
and optional real data: `examples/select_annotate_workflow_example.py`.

