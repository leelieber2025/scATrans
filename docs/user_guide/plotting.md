# Visualization

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
`save_path=` (300 dpi output).

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
multi-panel figures), `save_path=`, `show=`, `use_style=`, `figsize=` for
consistency. Most return `(fig, ax)` (or `(fig, axes_list)` for grids like
phase portraits).

## Main plotting functions

- `scat.pl.comet_plot(results_df, top_n=12, point_scale=1.0, min_size=2, max_size=180, s=None, ...)`
  Plots log fold change vs. bias-corrected unspliced excess residual
  (`unspliced_excess_residual`), sized and colored by `active_score`.
  - `s=3` (or 1-5): force **fixed** small point size for everything (direct, simple control).
  - `point_scale=0.2` + `min_size=1`: for variable sizing, make tiniest background points truly small.

- `scat.pl.volcano_plot(results_df, top_n=10, label_genes=None, style="auto", ...)`
  2D volcano (logFC vs. -log10(p_adj)).
  - **`style="auto"`** (default): scATrans legacy — `active_score` continuous colormap when present; otherwise up/down/ns.
  - **`style="ggvolcano"`**: [ggVolcano](https://github.com/BioSenior/ggVolcano) classic — teal Down / gray Normal / orange Up, `theme_bw` grid, dashed cutoffs, FDR-ranked labels, in-axes legend (`legend_position="UL"`).
  - **`style="gradual"`**: ggVolcano `gradual_volcano` — gradient color and point size by `-log10(FDR)`.
  - `label_genes=[...]` merges with `top_n` auto-labels; `label_by="p_adj"` (default for ggvolcano) or `"active_score"`.
  - `s=2` for fixed small points; `fills=` / `colors=` override the ggVolcano palette.

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

- `scat.pl.active_score_rankplot(results_df, top_n=20, ...)` — simple
  horizontal barplot of top active scores.

- `scat.pl.active_genes_heatmap(adata, genes, groupby=..., ...)` —
  convenience wrapper around `scanpy.pl.heatmap` for selected genes.

- `scat.pl.velocity_phase_portraits(adata, genes, groupby=..., ...)` — quick
  unspliced vs. spliced phase portraits for selected genes (useful for
  inspecting nascent excess).

- `scat.pl.set_style()` and `scat.pl.style_context()` — control global
  publication-style settings (vector fonts, minimal ink, etc.).

- `scat.pl.set_nature_style()` (legacy alias for `set_style`).
