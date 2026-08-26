# Core Workflow

If you have nascent layers and want the default DE-to-mechanism path, stay
here. No nascent layers? {doc}`standalone_de`. Brand new to the package?
{doc}`../quickstart`.

## Start here: `partition_de_by_mechanism`

For most mechanism analyses, call this once:

```python
import scatrans as scat

res = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    de="builtin",          # or method name / kwargs / DE DataFrame / callable
    sample_col="sample",   # set when you have biological replicates
    # gene_sets=my_pathways,
    # induction_matched=True,
    # add_nascent_score=True,
)
```

What it does, in order:

1. Checks nascent-layer quality (`res.regime`)
2. Runs DE and keeps changed genes (`padj < 0.05` and `logFC > 1.0` unless you
   change the cutoffs)
3. Scores the unspliced residual and labels mechanism
4. Optionally builds a program table (`gene_sets=`)
5. Optionally runs induction-matched program tests (`induction_matched=True`)

`res.summary()` prints the cutoffs that were actually used.

```python
res.regime                         # capture quality
res.selected                       # DE genes + mechanism labels
res.programs                       # if gene_sets=  (relative vs background)
res.programs_induction_matched     # if induction_matched=True
res.summary()
```

Internally this composes `active_score` → DE filter → mechanism helpers.
You do not need those pieces unless you want a custom path.

### Common options

| Option | When to use it |
|--------|----------------|
| `sample_col=...` | Biological replicates or libraries. Auto-pseudobulk + PyDESeq2 only if ≥3 samples per group; otherwise Wilcoxon stays cell-level |
| `de="builtin"` | Default DE inside the call |
| `de="wilcoxon"` / kwargs / DataFrame | A specific DE backend or your own DE table |
| `gene_sets={...}` | Program-level relative tests → `res.programs` |
| `induction_matched=True` | Induction strength varies a lot across genes |
| `add_nascent_score=True` | Extra detection columns only |
| Absolute placement | Separate call below (not a partition argument) |

### Absolute program placement (permutation-calibrated)

Competitive / induction-matched tables are **relative**. For absolute placement
on the mechanism axis (manuscript recommendation), subtract the same program’s
expectation under shuffled condition labels:

```python
de = res.gene_table[["logFC", "p_adj", "p_val"]]  # freeze membership
cal = scat.program_mechanism_permutation_calibrated(
    adata,
    gene_sets=my_programs,
    de=de,                         # required — freezes membership
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    restrict_to_selected=True,     # same membership as res.programs figures
    n_perm=200,                    # package default; floor p ≈ 1/201; raise for a finer grid
    random_state=0,
    # block_col="donor",           # optional: shuffle within batch/donor
)
print(cal[["observed_mean", "null_mean", "calibrated", "p_perm"]])
```

| Column | Meaning |
|--------|---------|
| `observed_mean` | Mean `transcription_support` for the program |
| `null_mean` | Same genes under shuffled labels (often ≠ 0) |
| `calibrated` | `observed − null` — absolute displacement (usual report) |
| `p_perm` | Permutation p (never zero; floor ≈ 1/(n+1)) |
| `null_sd` / `z` | Optional: offset stability / standardized offset |

Negative `calibrated` → more stabilization-weighted; positive → more
transcription-weighted. Pass the **same** DE table as the observed run so
membership stays fixed. Details: {doc}`../method`.

### Optional detection columns

Mechanism labels always use the induction-normalized residual. Separately,
you can add **detection** columns (absolute nascent increase):

| Parameter | Effect |
|-----------|--------|
| `add_nascent_score=False` (default) | No extra columns |
| `add_nascent_score=True` | Appends {func}`~scatrans.nascent_activity_score` columns |

| Column | Meaning |
|--------|---------|
| `nascent_poisson_z` | Detection score (not used for mechanism) |
| `dlog_unspliced` / `dlog_spliced` | CPM log fold-changes |
| `de_reproducible` / `de_repro_frac` | Spliced-side DE agreement flag (annotation only) |

Failures are recorded in `meta["nascent_score"]` without breaking mechanism
labels. Standalone:

```python
nz = scat.nascent_activity_score(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
```

Layers resolve automatically (`spliced`/`unspliced` or kb-python
`mature`/`nascent`).

---

The rest of this page is for custom or lower-level pipelines. Pure DE without
nascent layers: {doc}`standalone_de`. Scope questions: {doc}`../faq`.

## Lower-level: `active_score`

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    show_plot=True,
)
```

This computes DE, the reference-corrected unspliced residual (optional Huber
length/intron correction), and stores diagnostics in
`adata_res.uns["scatrans"]["diagnostics"]`.

The second return value (`significant`) is a strict DE and residual
conjunction. It is often empty. Use `all_results` with
`filter_active_genes(..., select_by="de")` for production gene lists, or
prefer `partition_de_by_mechanism`.

### Pseudobulk and DE method

With biological replicates:

```python
adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    use_pseudobulk=True,
    sample_col="sample",
    pseudobulk_de_backend="pydeseq2",  # or "scanpy"
    min_cells=5,
    min_counts=100,
    show_plot=True,
)
```

- Requires `sample_col`.
- `pydeseq2` needs raw counts (`store_raw_counts`) and
  `pip install "scatrans[pseudobulk]"`.
- `pseudobulk_de_backend="scanpy"` uses scanpy methods on aggregated profiles.
- Aggregation is **internal**. The returned `adata_res` stays cell-level
  (same `obs` columns, embeddings, and layers). Sample-level summary is in
  `adata_res.uns["scatrans"]["pseudobulk_obs"]`. The same contract applies
  to `differential_expression` and the `*_simple` / pipeline wrappers.

Change the scanpy test (cell-level or scanpy-on-pseudobulk):

```python
adata_res, significant, all_results = scat.active_score(
    ...,
    de_method="wilcoxon",
)
```

Settings are stored under `adata_res.uns["scatrans"]`.

### Choosing a DE backend

| Design | Suggested backend | Notes |
|--------|-------------------|--------|
| Quick look | scanpy Wilcoxon or t-test | Fast; cell-level pseudoreplication |
| ≥2 replicates / group | `use_pseudobulk=True` + PyDESeq2 | Needs raw counts |
| Pseudobulk without DESeq2 | `use_pseudobulk=True` + `scanpy` | Non-parametric on aggregates |
| Cell-level + true sample IDs | `use_mixed_model=True` + `sample_col` | Do not combine with Memento; see {doc}`advanced` |
| Method-of-moments cell DE | `use_memento_de=True` | Raw integer counts; exclusive with MixedLM |

```python
rec = scat.recommend_workflow(adata, groupby="condition", sample_col="sample")
# inspect rec and adata.uns["scatrans"]["diagnostics"] before publication claims
```

## Gene filtering with `filter_active_genes`

On the primary path, membership is already DE-only (`result.selected`). Use
this helper when you work from `all_results` or a pure DE table.

**Recommended:**

```python
de_list = scat.filter_active_genes(all_results, select_by="de")

candidates = scat.filter_active_genes(
    all_results,
    select_by="de",
    logfc_cutoff=0.3,
    padj_cutoff=0.05,
)
```

**Legacy / exploratory** (also applies residual cutoffs — not for production
discovery lists):

```python
candidates = scat.filter_active_genes(all_results, preset="heuristic")
# preset="significant" replays active_score's strict mask (needs permutation)
# preset="pseudobulk" / "permissive" also available
```

| Mode | Membership | Proxy gates |
|------|------------|-------------|
| `"de"` (**use this**) | DE only | Skipped |
| `"composite"` (API default for compatibility) | DE + residual gates | Applied |

Prefer `select_by="de"` even though the function default is still
`"composite"`. Partition and `run_default_pipeline(..., select_by="de")`
already use DE-only membership.

Direction for pure DE tables:

```python
down = scat.filter_active_genes(
    de_results,
    select_by="de",
    padj_cutoff=0.05,
    logfc_cutoff=0.3,
    logfc_direction="down",
)
both = scat.filter_active_genes(
    de_results,
    select_by="de",
    padj_cutoff=0.05,
    logfc_cutoff=0.3,
    logfc_direction="both",
)
```

Missing columns (e.g. residual FDR when permutation was off) are skipped
safely. Legacy names `velocity_residual` / `velocity_delta_raw` still exist as
aliases in `adata.var`.

## Design diagnostics

`diagnose_design` checks cell counts, replicates, and global unspliced
fraction. It runs automatically inside `active_score` when `sample_col` or
`use_pseudobulk=True`.

```python
diag = scat.diagnose_design(adata, groupby="condition", sample_col="sample")
print(diag["warnings"])
print(diag["recommendations"])
```

## Comparing DE across groups (e.g. cell types)

`compare_de_across_groups` runs the same contrast independently within every
level of a grouping column (cell type, cluster, ...) and summarizes
up/down-regulated gene counts for comparison:

```python
cmp = scat.compare_de_across_groups(
    adata,
    split_by="cell_type",          # loop over this column
    groupby="sample",              # the DE contrast
    target_group="RTX",
    reference_group="VEH",
    de_kwargs={                    # forwarded to differential_expression
        "use_pseudobulk": True,
        "sample_col": "individual",
        "pseudobulk_de_backend": "pydeseq2",
    },
    padj_cutoff=0.05,
    logfc_cutoff=0.25,
)
cmp.summary                  # per-cell-type up / down / total counts
cmp.up["T"]                  # upregulated genes table for cell type "T"
cmp.failed                   # {level: error} for levels DE could not run on
scat.pl.de_summary_barplot(cmp, mode="stacked", sort_by="total")
scat.pl.composition_barplot(cmp.results["T"], groupby=..., hue="mechanism_class")
```

Levels with fewer than `min_cells` cells, or where DE raises (e.g. too few
replicates for pseudobulk dispersion estimation), are skipped and recorded in
`cmp.failed` rather than aborting the whole comparison
(`raise_on_error=True` to fail fast instead). Significant up/down genes reuse
`filter_active_genes(..., select_by="de")` internally, so the same
`padj_cutoff` / `p_type` / `logfc_cutoff` semantics apply.

### Synthetic pseudo-replicates

Pseudobulk DE (PyDESeq2 dispersion estimation, `use_mixed_model`, ...) needs
more than one "sample" per condition. If you have no real donor/individual
column, `assign_pseudo_replicates` randomly splits the cells of each
condition into `n_replicates` synthetic subgroups:

```python
scat.assign_pseudo_replicates(
    adata, groupby="sample", n_replicates=3, key_added="individual", random_state=0,
)
```

**These are pseudo-replicates, not biological replicates.** Random splits of
the same cells only resample cell-level (technical/sampling) noise — they do
not capture true donor-to-donor variance, so pseudobulk p-values on this
column remain anti-conservative relative to a real multi-donor design. Use
real donor/sample identifiers whenever they exist; reach for this only when
none are available and you need an approximate multi-sample pseudobulk run.

## Layer requirements

| Goal | Layers | Entry point |
|------|--------|-------------|
| Mechanism partition | `spliced`/`unspliced` or `mature`/`nascent` | {func}`~scatrans.partition_de_by_mechanism` |
| Residual / diagnostics only | Same | {func}`~scatrans.active_score` / {func}`~scatrans.active_score_simple` |
| DE + enrichment only | Counts in `.X` (raw snapshot recommended) | {func}`~scatrans.differential_expression` |

## Regime check

```python
r = scat.qc.regime_diagnosis(adata)
print(r["regime"], r["reliability"], r["message"])
```

Partition always stores this as `result.regime` and scales mechanism
confidence. Low reliability is a caution, not a hard stop.
