# Visualization

Worked gallery (synthetic data, nearly all `scat.pl.*` functions):
{doc}`../tutorials/t_synthetic_visualization`.

```python
import scatrans as scat

scat.pl.comet_plot(all_results, top_n=12, title="Active Drivers")
scat.pl.volcano_plot(all_results, top_n=10, label_genes=["YourGene1", "YourGene2"])
scat.pl.bias_diagnostic_plot(all_results)
```

**ggVolcano-style volcano plots**
([BioSenior/ggVolcano](https://github.com/BioSenior/ggVolcano)) are available
via `style=`:

```python
# Classic three-color volcano (Down=teal, Normal=gray, Up=orange; theme_bw; labels by FDR)
scat.pl.volcano_plot(
    all_results,
    style="ggvolcano",
    top_n=12,
    logfc_cutoff=0.35,
    pval_cutoff=0.05,
    legend_position="UL",  # UL / UR / DL / DR
    save_path="volcano_ggvolcano.png",
)

# Gradient fill + point size by -log10(FDR) (gradual_volcano)
scat.pl.volcano_plot(all_results, style="gradual", top_n=10)

# Legacy scATrans look (active_score colormap when present) — default
scat.pl.volcano_plot(all_results, style="auto")
```

`style="ggvolcano"` labels the top `top_n` genes by smallest `p_adj` (FDR),
and accepts `label_genes=[...]` for manual labels. Custom palettes: `fills=`
/ `colors=` (Down, Normal, Up hex tuples).

All plotting functions support `ax=` / `axes=` for multi-panel figures and
`save_path=` (files are written at **≥300 dpi** even when on-screen `dpi`
is lower).

## Display defaults (notebook-first)

Defaults are tuned for notebooks and docs, not full-page journal panels:

| Knob | Default | Paper / export |
|------|---------|----------------|
| `dpi` | **150** (figure) | `context="paper"` → 300, or set `dpi=` |
| `figsize` | modest (~6×4.5) | `context="paper"` → ~8×6 |
| `fontsize` / gene labels | **10** / **8**, normal weight | larger under `context="paper"` |
| `top_n` (comet / volcano) | **8** / **6** | raise if you want more labels |
| `save_path` | — | always ≥**300** dpi on disk |

```python
scat.pl.volcano_plot(all_results)                 # notebook defaults
scat.pl.volcano_plot(all_results, context="paper")  # larger fonts / dpi=300
scat.pl.volcano_plot(all_results, save_path="fig.pdf")  # sharp export
```

Default volcano **`style` remains `"auto"`** (color by `active_score` when
present). For a classic two-sided teal/orange look use `style="ggvolcano"`.
Dense tables without `s=` log a warning suggesting fixed small points.

## Plotting style

```python
import scatrans as scat
scat.pl.set_style()                 # once early (opt-in)
# or (to limit scope):
with scat.pl.style_context(linewidth=0.8):
    scat.pl.comet_plot(...)         # inside block or pass use_style=True
# Default for pl.* functions is use_style=False (prevents surprising rcParams changes in notebooks).
```

All `scat.pl.*` functions support `ax=` / `axes=` (for embedding in
multi-panel figures), `save_path=`, `show=`, `use_style=`, `figsize=`, and
(on major plotters) `context=` for consistency. Most return `(fig, ax)`
(or `(fig, axes_list)` for grids like phase portraits). When embedding
`enrich_dotplot(..., ax=...)`, the size legend is drawn **inside** the axes
so multipanel layouts do not clip an exterior legend.

## Main plotting functions

- `scat.pl.comet_plot(results_df, top_n=8, point_scale=1.0, min_size=2, max_size=180, s=None, ...)`
  Plots log fold change vs. bias-corrected unspliced excess residual
  (`unspliced_excess_residual`), sized and colored by `active_score`.
  - `s=3` (or 1-5): force **fixed** small point size for everything (direct, simple control).
  - `point_scale=0.2` + `min_size=1`: for variable sizing, make tiniest background points truly small.

- `scat.pl.volcano_plot(results_df, top_n=6, label_genes=None, style="auto", ...)`
  2D volcano (logFC vs. -log10(p_adj)).
  - **`style="auto"`** (default): scATrans legacy — `active_score` continuous colormap when present; otherwise up/down/ns.
  - **`style="ggvolcano"`**: [ggVolcano](https://github.com/BioSenior/ggVolcano) classic — teal Down / gray Normal / orange Up, `theme_bw` grid, dashed cutoffs, FDR-ranked labels, in-axes legend (`legend_position="UL"`).
  - **`style="gradual"`**: ggVolcano `gradual_volcano` — gradient color and point size by `-log10(FDR)`.
  - `label_genes=[...]` merges with `top_n` auto-labels; `label_by="p_adj"` (default for ggvolcano) or `"active_score"`.
  - `s=` / `point_scale` / `min_size` / `max_size` for size; `color_by=` + `cmap=` (auto) or `fills=` / `colors=` (ggvolcano).

- `scat.pl.bias_diagnostic_plot(results_df, point_size=10, ...)`
  Before/after view of the effect of length+intron bias correction on the
  velocity delta. `point_size` controls the gene cloud density.

- `scat.pl.volcano_3d(results_df, point_scale=..., min_size=2, s=None, ...)`
  3D version of the volcano. Same size controls (`s` for fixed size).

- `scat.pl.enrich_dotplot(enrich_df, ...)` works well with GSEA results too
  (auto defaults to `x="NES"`, diverging cmap for `color_by="NES"`).

- `scat.pl.gseaplot(ranked_genes, gsea_result, term=...)` — classic GSEA
  running-sum plot (uses precomputed curves from `run_gsea` when available).

- `scat.pl.enrich_dotplot(enrich_df, top_n=15, show_terms=None, x="GeneRatio", size_by="Count", color_by="Adjusted P-value", ...)`
  Enrichment dot plot (clusterProfiler style).
  - `x`: x-axis variable — "GeneRatio" (default for ORA), "FoldEnrichment", **"Count"**, "-log10(p.adj)", or "NES" (for GSEA).
  - `size_by` (dot size, default "Count"), `color_by` (default adjusted p-value; "NES" for GSEA uses diverging colormap).
  - `show_terms` accepts int (top N), "auto" (p.adjust <0.05 + Count>=2 smart selection), or list of term strings/Descriptions (exact or partial match, order preserved) — directly analogous to `dotplot(..., showCategory=...)`.
  - Also available as `enrich_barplot`.

- `scat.pl.compare_dotplot(enrich_df, top_n=5, cluster_col=None, size_by="GeneRatio", color_by="p.adjust", ...)`
  Multi-group comparison grid in the style of clusterProfiler
  `dotplot(compareCluster(...))`: **groups as columns** (x-axis), **terms as
  rows** (y-axis), with dot size/color encoding enrichment strength.
  - Designed for long tables from `scat.compare_enrichment` /
    `scat.concat_compare_results` (any table with a `Cluster` column plus
    `Term`/`Description` and `p.adjust`).
  - `top_n` keeps the top terms **per group** before taking the union across
    groups; `show_terms=` can pin an explicit term list.
  - Use `enrich_dotplot(..., facet_by_cluster=True)` when you prefer separate
    panels rather than a single compare grid.

- `scat.pl.enrich_upsetplot(enrich_df, ...)` / `scat.pl.enrich_vennplot(enrich_df, ...)`
  Set-overlap views across clusters/contrasts (significant **terms** by
  `pval_cutoff` / `min_count`). Prefer UpSet for more than three groups.

- `scat.pl.gene_upsetplot(...)` — the **gene-level** UpSet: how genes overlap
  across several DE results or gene lists (companion to the term-level
  `enrich_upsetplot`). Pure matplotlib, no external `upsetplot` dependency. The
  workflow is three small pieces:

  ```python
  # 1) tidy multiple DE results into a gene x set membership matrix
  mem = scat.pl.build_gene_membership(
      {"wilcoxon": de_wilcox, "ttest": de_ttest, "pseudobulk": de_pb},
      direction="separate",          # each DE -> name::up and name::down
      pval_cutoff=0.05, logfc_cutoff=0.5,
  )
  # 2) draw it (common-up and common-down show up as their own columns)
  scat.pl.gene_upsetplot(membership=mem, save_path="gene_upset.png")
  # 3) pull the intersection genes back out for enrichment
  up   = scat.pl.common_genes(mem, direction="up")      # up in *every* method
  down = scat.pl.common_genes(mem, direction="down")
  enr  = scat.run_enrichment(up, adata=adata)           # straight into ORA
  ```

  `build_gene_membership` also accepts ready-made `{name: [gene, ...]}` lists
  (no thresholds applied). `common_genes(..., min_sets=2)` relaxes the strict
  intersection to a "recovered by at least *k* sets" signature, and `sets=[...]`
  intersects an explicit subset of columns. Every visual element is recolorable
  (`set_color`, `intersection_color`, `dot_color`, `inactive_color`,
  `line_color`); pass a per-column list to `intersection_color` / `dot_color` to
  highlight specific intersections.

- `scat.pl.active_score_rankplot(results_df, top_n=15, context=None, ...)` —
  horizontal bar plot of top active scores (gradient fill by magnitude).
  Default `top_n` is 15 under notebook display defaults.

- `scat.pl.active_genes_heatmap(adata, genes, groupby=..., ...)` —
  convenience wrapper around `scanpy.pl.heatmap` for selected genes.

- `scat.pl.velocity_phase_portraits(adata, genes, groupby=..., ...)` — quick
  unspliced vs. spliced phase portraits for selected genes (useful for
  inspecting nascent excess).

- `scat.pl.gamma_shrinkage_plot(results_df, ...)` — empirical-Bayes gamma
  shrinkage weight vs expression depth. Requires a
  `gamma_shrinkage_weight` column (`gamma_method="empirical_bayes"`). Recolor
  with `cmap=` (when `effective_gamma` is present) or `color=` (single-color
  fallback).

- `scat.pl.set_style()` and `scat.pl.style_context()` — control global
  publication-style settings (vector fonts, minimal ink, etc.).

- `scat.pl.set_nature_style()` (legacy alias for `set_style`).

- `scat.pl.figure_export_context(directory=...)` /
  `scat.pl.save_all_figures(figures, directory=...)` — batch figure export
  helpers for multipanel manuscript figure packs.

Label repel on comet/volcano uses optional `adjustText`
(`label_repel=True` by default). Pass `label_repel=False` to skip if the
dependency is unavailable or labels should stay at data coordinates.
`comet_plot(..., positive_logfc_only=True)` (default) restricts the classic
active-driver view to upregulated genes; set `False` for the full logFC
range.