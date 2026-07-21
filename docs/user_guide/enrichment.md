# Functional Enrichment

Over-representation analysis with `run_enrichment`:

```python
enrich_res = scat.run_enrichment(
    gene_list=candidates.index.tolist(),
    gene_sets="GO_Biological_Process",   # or "GO_BP" — automatically resolved to the
                                         # correct organism-specific built-in (Hs/Mm_GO_..._2026)
    organism="mouse",  # or "human"
    adata=adata,   # if you called store_raw_counts(adata) earlier, this will
                   # automatically use the preserved full measured gene list as universe.
                   # Explicit `universe=` still takes precedence.
    padj_cutoff=0.05,  # preferred; legacy pval_cutoff= still accepted (filters adjusted p)
    min_size=5,
    max_size=500,
)
# Additional columns and attrs (clusterProfiler compatibility):
#   - "neg_log10_padj" column
#   - res.attrs["universe_info"] with effective_universe_size, dropped_by_annotation_filter, etc.
```

## `run_gsea` (pre-ranked GSEA)

For ranked-list enrichment (classic preranked GSEA; Subramanian et al. 2005 —
see {doc}`../references`), via [GSEApy](https://github.com/zqfang/GSEApy):

```python
# Prefer a *signed* ranking metric (logFC). Non-negative scores such as
# active_score cannot produce negative NES (depletion) and are not auto-selected.
ranked = all_results["logFC"].sort_values(ascending=False)

gsea_res = scat.run_gsea(
    ranked_genes=ranked,  # or all_results (auto-picks logFC when present)
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

GSEA needs **signed** ranks. Passing a full `all_results` table auto-prefers
`logFC` (and similar t-stat columns); `active_score` is never auto-selected.
If you force `score_column="active_score"` (or pass a one-sided Series), a
warning is emitted.

**Low mapping rate and ID cleanup (same gate as ORA):** after loading gene
sets, ranked genes are intersected with gene-set members via
`_check_gene_set_mapping_rate`. Mapping rate **&lt; 20%** emits a
`UserWarning` with input-symbol and gene-set examples; **0% overlap**
returns an empty frame with `reason="no_ranked_genes_mapped"` (avoids
opaque gseapy failures from case/ID mismatch). Duplicate gene IDs (common
after case-folding or multi-mapped symbols) are collapsed by keeping the
score with **largest absolute value** (`max |score|` per gene). Mismatched
case is a common failure mode with Enrichr libraries (UPPERCASE) vs mouse
symbols (`Tp53`) — pass `gene_case="upper"` when needed.

`run_gsea` stores the full enrichment score curves in
`.attrs["gsea_details"]` so that `gseaplot` renders exactly the same RES that
produced the reported NES/p-values.

Requires the optional extra:

```bash
pip install "scatrans[gsea]"   # pulls in gseapy
```

## `run_kegg` (convenience wrapper for KEGG pathways)

```python
kegg_res = scat.run_kegg(
    gene_list=candidates.index.tolist(),
    organism="mouse",  # or "human"
    # Defaults to the organism-specific built-in library (Hs_KEGG_2026 or Mm_KEGG_2026)
    adata=adata,   # if store_raw_counts was called earlier, this automatically uses
                   # the preserved full measured gene set as background.
    padj_cutoff=0.05,
)
```

## Using bundled gene sets

The package defaults to organism-specific bundled sets. Use `organism=`
together with base names such as `"GO_Biological_Process"` or `"KEGG"`.
Supply a full historical name (e.g. `"GO_Biological_Process_2023"`) to
select an Enrichr version.

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

## Using original Enrichr versions

To use a specific historical Enrichr/gseapy library, pass the full versioned
name (including the year). Versioned Enrichr-style names load through gseapy.

```python
# Specific Enrichr version for KEGG
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

**Adding your own sets**: Drop `.gmt` files into `src/scatrans/data/`. See
`src/scatrans/data/README.md`.

## `simplify_enrichment` (reduce redundant enrichment terms)

Two methods are supported:

- **`jaccard`** (default): greedy filtering by Jaccard overlap of enriched gene lists.
- **`pathway_denester`**: combinatorial nested-pathway test adapted from
  [PathwayDenester](https://github.com/Helmy-Lab/PathwayDenester). Better at
  removing terms that are significant only because they are nested inside a
  more significant parent pathway. Requires full pathway gene memberships
  (auto-loaded from `enrich_res.attrs` when enrichment used bundled/Enrichr
  libraries; pass `gene_sets=` again if you used a custom dict).

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

## `run_go` (GO enrichment, clusterProfiler-style)

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

`run_go` automatically resolves to the organism-specific bundled sets when
possible (BP is bundled; CC/MF fall back to gseapy/Enrichr if the library is
installed).

## Exporting results

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

# Long-format term–gene table (one row per gene)
long_table = scat.expand_enrichment_genes(res)
# If the input was from run_go(ontology="ALL"), long_table will have an "Ontology" column first.
```

`save_enrichment_report` also writes a rich `metadata.json` (and a
"metadata" sheet in the xlsx) containing:

- `analysis_info` (package, version, timestamp)
- `gene_set_info` (requested/resolved + `requested_source` vs `actual_source`: "bundled", "gseapy", "gmt", "dict")
- `universe_info` (effective N, dropped genes, restrict behavior, etc.)
- Full `.attrs` from the enrichment call (including per-ontology details for GO ALL)

Empty results still carry diagnostic `.attrs` (`reason`, `gene_set_info`,
`universe_info`, and related fields).

## Additional enrichment plot options

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

# Multi-group compareCluster-style grid (groups as columns, terms as rows).
# Input: long table from compare_enrichment / concat_compare_results.
scat.pl.compare_dotplot(
    compare_df,          # must include Cluster + Term/Description + p.adjust
    top_n=5,             # top terms kept *per group*, then unioned
    size_by="GeneRatio",
    color_by="p.adjust",
    save_path="compare_dotplot.pdf",
)

# Overlap views across contrasts
scat.pl.enrich_upsetplot(compare_df, pval_cutoff=0.05)
scat.pl.enrich_vennplot(compare_df, pval_cutoff=0.05)  # best for 2–3 groups

# Embed in multi-panel figure with ax=
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 5))
scat.pl.enrich_dotplot(enrich_res, top_n=8, ax=ax, show=False)
fig.savefig("enrich_panel.pdf", dpi=300, bbox_inches="tight")
```

Additional options:

- `show_terms=15` or `show_terms="auto"` or `show_terms=["term A", "term B"]`
- `use_style=True` to apply publication style for that call only
- `context="paper"` on major plotters for journal-sized defaults

All `scat.pl.*` functions accept `save_path`, `ax`/`axes`, `figsize`, `show`,
and `dpi`. Files written via `save_path` are always exported at ≥300 dpi.
See {doc}`plotting` for display defaults and style helpers.