# Quickstart

## Minimal API (recommended default path)

If you want the recommended default path without dozens of parameters, use
the simple wrappers:

```python
import scatrans as scat

# One-liner pipeline: score → filter → GO enrichment
result = scat.run_default_pipeline(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",   # optional; auto-selects pseudobulk when >=3 replicates/group
    organism="mouse",
)
print(result["candidates"].head())
print(result["enrichment"].head())

# Or just the core scoring step:
adata_res, significant, all_results = scat.active_score_simple(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",
)
```

`active_score_simple` / `run_default_pipeline` auto-attach gene features,
pick Wilcoxon (single-cell) or pseudobulk+PyDESeq2 (when replicates allow),
and keep permutation off by default. Use `active_score(...)` directly for
advanced options (permutation, mixed models, Memento, etc.) — see the
{doc}`user_guide/index`.

## Complete end-to-end example

This is a complete, copy-paste friendly workflow for first-time users. It
takes you from loaded data to differential results, filtering, enrichment,
and visualization of enrichment results.

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

You can now explore `all_results`, adjust filters in step 6, try different
`gene_sets`, or run `run_go` / `run_gsea` (see {doc}`user_guide/enrichment`).

For pure differential expression without spliced/unspliced layers, replace
step 5 with `scat.differential_expression(...)` — see
{doc}`user_guide/standalone_de`.

## Preserving raw counts and layers

Call `store_raw_counts` early (after loading and QC, before HVG or
normalization). It writes the current `.X` to `layers["counts"]` and copies
the original spliced/unspliced layers. These survive later subsetting and
provide the correct background for enrichment and count-based DE.

The default `save_raw=False` avoids populating `adata.raw`.

After HVG-based visualization on a copy, restore or use the preserved layers
for full-gene DE, active scoring, or enrichment (pass `adata=` to
`run_enrichment` or `run_kegg` to use the stored gene list as background).

HVG subsetting also subsets the saved layers. This keeps velocity
calculations consistent with `.X`. To analyze more genes than the HVG set,
store before subsetting or operate on the unfiltered object for DE and
enrichment steps.

To restore raw counts into `.X` for the current gene set:

```python
adata_raw = scat.restore_raw_counts(adata, layer="counts", inplace=False)
```

See {doc}`user_guide/standalone_de` for the no-velocity use case.
