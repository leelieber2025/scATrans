# Quickstart

:::{important}
**Scope.** Gene-list membership always comes from **DE**; the nascent
(spliced/unspliced) signal only **annotates** mechanism and is a low-confidence
hint per gene — conclude at the **program** level and heed the reliability
pre-flight (weak on low-capture / 3′-biased data). Differential expression,
functional enrichment, and plotting that do **not** rely on velocity layers are
stable; see {doc}`user_guide/standalone_de` for a pure-DE path.
:::

## Minimal API (recommended default path)

The recommended entry point is {func}`~scatrans.partition_de_by_mechanism`: a
standard DE test **selects** the changed genes, then scATrans **partitions** them
into transcription-driven vs stabilization-driven — a soft per-gene annotation and
a decisive **program-level** call. It runs the mandatory reliability pre-flight for
you.

```python
import scatrans as scat

# DE selects membership; scATrans annotates mechanism (never filters on the proxy).
result = scat.partition_de_by_mechanism(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    organism="mouse",
    de="builtin",            # or "wilcoxon"/"t-test"/pseudobulk kwargs, a precomputed
                             # DE DataFrame (de_logfc_col/de_padj_col), or a callable
    gene_sets=my_pathways,   # optional -> program-level transcription-vs-stabilization table
)
print(result.regime)         # reliability pre-flight (regime_diagnosis)
print(result.selected.head())    # DE-selected genes + transcription_support / mechanism_class
print(result.programs)       # decisive program-level calls (when gene_sets given)

# Just the core scoring step (advanced):
adata_res, significant, all_results = scat.active_score_simple(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
```

`partition_de_by_mechanism` / `active_score_simple` auto-attach gene features,
pick Wilcoxon (single-cell) or pseudobulk+PyDESeq2 (when replicates allow), and
keep permutation off by default. The reliability pre-flight
(`scat.qc.regime_diagnosis`) scales every per-gene `mechanism_confidence`.

:::{note}
The legacy composite pipeline `scat.run_default_pipeline(..., select_by="composite")`
(the default of `run_default_pipeline`) still works but now emits a
`DeprecationWarning` — its composite `active_score` mixes the DE and proxy legs
and does not out-discover DE. For a **pure DE gene list** without the mechanism
layer, use `run_default_pipeline(..., select_by="de")` or
{doc}`user_guide/standalone_de`. Use `active_score(...)` directly for advanced
options (permutation, mixed models, Memento) — see the {doc}`user_guide/index`.
:::

## Complete end-to-end example (lower-level `active_score` path)

This walks the **lower-level** `active_score` → `filter_active_genes` → enrichment
path in detail (QC, raw-count preservation, gene features, scoring, filtering,
plotting). For the recommended primary workflow use
{func}`~scatrans.partition_de_by_mechanism` (above) — it composes these steps and
adds the mechanism partition. The example below is useful for understanding the
building blocks and for the residual/scoring diagnostics `active_score` computes.

```python
import scanpy as sc
import scatrans as scat

# 1. Load your data (must contain spliced/unspliced layers or use differential_expression instead)
adata = sc.read_h5ad("your_data.h5ad")

# 2. Store raw counts + original layers early (before HVG/normalization).
#    sidecar=True (default) also snapshots the full-gene counts + velocity layers
#    into .uns so they survive later HVG/cell subsetting.
scat.store_raw_counts(adata, layer="counts")

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
    # padj_cutoff=0.05,           # preferred; legacy pval_cutoff= still works
)

print(f"\nFiltered candidate genes: {len(candidates)}")

# 7. Functional enrichment (GO)
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",   # or "GO_BP"
    organism="mouse",                    # or "human"
    adata=adata,                         # uses stored raw genes as background
    padj_cutoff=0.05,
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

You can now explore `all_results`, adjust filters in step 6, try different
`gene_sets`, or run `run_go` / `run_gsea` (see {doc}`user_guide/enrichment`).

For pure differential expression without spliced/unspliced layers, replace
step 5 with `scat.differential_expression(...)` — see
{doc}`user_guide/standalone_de`.

## Preserving raw counts and layers

Call `store_raw_counts` early (after loading and QC, before HVG or
normalization). It does two things:

- writes the current `.X` to `layers["counts"]` (axis-aligned, like any layer);
- with `sidecar=True` (the default), also writes a **label-indexed snapshot** of
  the full obs × var counts — plus any velocity layers
  (`spliced`/`unspliced` or `mature`/`nascent`) — to
  `adata.uns['scatrans']['raw_snapshot']`.

The counts layer and velocity layers follow normal AnnData behavior: HVG or cell
subsetting trims them. The **snapshot** lives in `.uns`, so it is *not* tied to
the obs/var axes — it survives HVG subsetting, cell subsetting, `copy()`, and
`write_h5ad()`, and is aligned back by cell/gene name.

To restore raw counts into `.X` for the **current** gene set (also absorbs cell
subsetting and gene reordering):

```python
adata_raw = scat.restore_raw_counts(adata, inplace=False)
```

To recover the **full pre-HVG gene universe** (counts + velocity layers) as a new
AnnData — for full-gene DE, active scoring, or enrichment — even after HVG or cell
subsetting:

```python
adata_full = scat.restore_raw_counts(adata, full_genes=True)
```

For enrichment you can also just pass `adata=` to `run_enrichment` / `run_kegg`
to use the stored full gene list as the background universe.

For large datasets, use `store_raw_counts(adata, sidecar="ondisk",
snapshot_path="raw_snapshot.h5ad")` to keep the full matrix on disk and only a
lightweight pointer in `.uns`.

```{note}
`save_raw=True` is deprecated: `adata.raw` is commonly reserved for
log-normalized data, and the sidecar snapshot already preserves full-gene raw
counts. Use `restore_raw_counts(..., full_genes=True)` instead.
```

See {doc}`user_guide/standalone_de` for the no-velocity use case.
