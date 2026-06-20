# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.2] - 2026-06-20

### Added / Changed
- `filter_active_genes` now accepts `logfc_direction="up"|"down"|"both"` (default remains `"up"` for backward compatibility and "active" semantics).
  - `logfc_cutoff` is interpreted as a positive magnitude in all modes.
  - `"down"`: selects logFC < -cutoff (downregulated genes from differential_expression results).
  - `"both"`: selects |logFC| > cutoff.
  - Sorting for pure-DE tables is now direction-aware (most-negative-first for down, largest |logFC| for both).
  - This directly supports the common request for downregulated candidates in standalone DE workflows.
- Updated docstrings, tests, and examples.

## [0.9.0] - 2026-06-19

### Added
- `run_gsea(ranked_genes, gene_sets=..., nperm=..., ...)` — pre-ranked GSEA (via gseapy.prerank wrapper).
  Reuses the same gene-set loading, `gene_case`, diagnostics, and `.attrs` system as ORA.
  Returns DataFrame with `NES`, `ES`, `pvalue`, `p.adjust`, `leading_edge`, etc.
  Optional dependency: `pip install "scatrans[gsea]"`.
- `scat.pl.gseaplot(ranked_genes, gsea_result, term=...)` — classic GSEA running enrichment score plot.
  Automatically uses pre-computed RES curves + hits stored in `run_gsea` results (`.attrs["gsea_details"]`).
- `enrich_dotplot` now auto-detects GSEA results (defaults `x="NES"`, uses diverging colormap when `color_by="NES"`).
- Added `gsea` extra in `pyproject.toml`.

### Changed
- Minor internal cleanups and test coverage for the new GSEA path.
- All new functions follow the existing consistent signatures (`ax=`, `use_style=`, `save_path=`, etc.).

## [0.8.0] - 2026-06-14

### Added (enrichment module — major paper-readiness upgrade)
- `run_go(ontology="BP"|"CC"|"MF"|"ALL", ...)` — direct wrapper analogous to clusterProfiler `enrichGO`. Supports `adjust_across_all=True` for a single BH correction across all GO terms when using "ALL".
- `save_enrichment_report(res, prefix=..., save_excel=True, save_csv=True, save_tsv=True, save_metadata=True, save_term_gene_table=True)` — one-call export of main table, term-gene long table (via `expand_enrichment_genes`), and rich `metadata.json` + xlsx sheet. Auto-creates parent directories. List columns (e.g. `Genes_list`) are sanitized to `;` strings for clean export.
- `expand_enrichment_genes(res)` — expands the `Genes` (semicolon) column into a long-format Term–Gene table (one row per gene). Preserves `Ontology` column when input came from `run_go(..., "ALL")`.
- Rich provenance in every result `.attrs` (success and empty):
  - `analysis_info`: package, version, timestamp, module
  - `gene_set_info`: `requested`/`resolved`, `requested_source` vs `actual_source` ("bundled", "gseapy", "gmt", "dict"), `library_name`, `n_terms`, `n_unique_genes`
  - `universe_info`: full details of background handling (provided size, restricted, dropped_by_annotation, force_universe, mapping counts)
  - Empty results now carry `reason` ("gene_list_empty", "universe_empty", "no_term_overlap_after_filters", ...) + the above fields so users can diagnose why nothing came back.
- New `run_enrichment` / `run_kegg` / `run_go` parameters: `padj_cutoff` (preferred modern name), `include_gene_list` (adds `Genes_list` python-list column), `adjust_across_all`.
- `list_bundled_gene_sets()` now clearly documents the 2026 organism-specific defaults.
- Improved low-mapping-rate warning (includes input examples + gene-set examples).
- `background` is now a documented deprecated alias of `universe`; passing both raises immediately.
- All empty-result DataFrames preserve consistent columns (including optional `Genes_list` when requested) and full diagnostic attrs.

### Changed / Improved
- `_load_gene_sets` now returns `(term_to_genes, term_to_desc, load_info)` so `actual_source` is always recorded accurately (even on gseapy fallback after bundled attempt).
- `run_kegg` fully synchronized with new parameters (`padj_cutoff`, `include_gene_list`, etc.).
- `enrich_dotplot` (pl.py) and various tl.py flows updated for new columns/attrs.
- Version unified to 0.8.0 for this release.
- README and docstrings extensively updated with manuscript-export examples, `run_go`, provenance details, and `adjust_across_all` guidance.
- Full test coverage for new paths (per-ontology attrs, within_ontology p.adjust, save+tsv+dir creation, expand with Ontology, dual-cutoff warning, etc.). All tests pass.

## [0.9.0] - 2026-06-18

### Added
- **Independent permutation statistics for unspliced excess**: `unspliced_excess_pval` and `unspliced_excess_fdr` (one-sided test on bias-corrected `unspliced_excess_residual`). Computed alongside existing `active_score_pval` / `active_score_fdr` when `use_permutation=True`.
- New parameter `unspliced_excess_fdr_cutoff` (default 0.05) for the built-in `significant` gene list and `filter_active_genes`.
- `filter_active_genes` parameters `unspliced_excess_residual_cutoff` and `unspliced_excess_fdr_cutoff`; heuristic/pseudobulk presets updated accordingly.
- `adata.uns["scatrans"]["significant_criteria"]` metadata documenting the built-in significance conjunction.

### Changed
- **Terminology**: primary result columns renamed from velocity to unspliced/nascent excess:
  - `unspliced_excess_delta` (was `velocity_delta_raw`)
  - `unspliced_excess_residual` (was `velocity_residual`)
  - Legacy `velocity_*` columns remain in `adata.var` as deprecated aliases.
- **Built-in `significant` gene list** now requires:
  - `logFC > logfc_cutoff`, `p_adj < pval_cutoff`, `unspliced_excess_residual > 0`, `unspliced_excess_fdr < unspliced_excess_fdr_cutoff`
  - `active_score` is no longer used for significance (ranking/visualization only).
  - Without `use_permutation=True`, the built-in `significant` list is empty (logged warning).
- Plotting functions accept primary or legacy column names; axis labels updated.
- README rewritten for the new significance model and column names.

### Deprecated
- `active_fdr_cutoff` (no longer used for built-in significance; use `unspliced_excess_fdr_cutoff`).
- `velocity_residual_cutoff` in `filter_active_genes` (use `unspliced_excess_residual_cutoff`).

## [Unreleased]

### Added
- **Rich runtime diagnostics** (high priority improvement): `active_score` now automatically computes global unspliced fraction (via integrated `qc.unspliced_global`), captures detailed bias-correction fit results (coefficients, n_genes_used, fallback status), stores per-gene `effective_gamma`, and records everything under `adata.uns["scatrans"]["diagnostics"]`. A concise run summary is logged at completion.
- New public helper `scat.pl.velocity_phase_portraits(adata, genes, groupby=...)` for quick visual inspection of U vs S relationships on top hits (lower-priority diagnostic aid).
- `examples/real_data_template.py` — heavily commented, non-runnable but copy-adaptable template demonstrating the full recommended real-data workflow, QC, diagnostics inspection, and publication plotting.
- Explicit documentation of the permutation approximation (velocity layers fixed from original data; only labels shuffled) in code metadata, logs, and the new "Choosing mode" section of the README.
- `effective_gamma` column is now added to `adata.var` (and included in result tables when present) for transparency of the gamma used in the delta calculation.

### Changed / Improved
- Major usability & paper-readiness upgrade to diagnostics, metadata, and user guidance.
- `qc.unspliced_global` is now called automatically inside `active_score` (result stored); the function remains directly usable as a pre-flight check (`scat.qc.unspliced_global(adata)`).
- Added prominent "Choosing `mode`: heuristic vs advanced (and common pitfalls)" section + decision guide to README.

### Added (0609 refactor)
- Mixed model support (`use_mixed_model=True` + `sample_col`): statsmodels LMM (~ condition + (1 | sample)) as cell-level replicate-aware DE backend (lightweight Python analogue to dreamlet/variancePartition/dreampy LMM + Libra mixed options). Replaces scanpy/pydeseq2 p/logFC when enabled.
- Delta Variance: `delta_variance` (condition-attributable variance fraction, 0-1, variancePartition-style) and `delta_var_pval` (LRT) computed during mixed fit; always surfaced in `all_results` / adata.var / diagnostics when mixed used.
- New `active_score` options: `use_delta_variance_pval` (bool) + `delta_var_pval_cutoff` to optionally include the delta LRT pval as supplementary filter in significant gene selection.
- Full backward compat; mixed path documented vs. pseudobulk; guidance + references to Libra, dreampy, NEBULA, dreamlet in README.
- Existing package backed up to `backup/backup0609/`.
- Tests + direct verification cover mixed, delta col, filter option, incompatibility with pseudobulk, and no regression on legacy paths. All tests pass.
- Enhanced `README.md` with real-data workflow pointers and stronger emphasis on inspecting diagnostics.
- `scat.pl` now documents the new phase-portrait helper; all existing `ax=` support preserved.
- Internal: `_fit_huber_bias_correction` and velocity helpers now return extra fit/quality information (used for the new diagnostics) while preserving full backward compatibility of public results.

### Changed
- **Major internal refactor** (2025): Core logic in `tl.py` (`active_score`) was extracted into private supporting modules (`_utils.py`, `_de.py`, `_velocity.py`, `_bias.py`, `_permutation.py`). The public `active_score` function is now a thin, readable orchestrator while preserving 100% identical behavior, return values, and side effects on `adata`.
- Logging discipline: Library code in `pp_bias.py` and the gene-features generator now uses the `scatrans` logger instead of direct `print()` calls with emojis. CLI entrypoint configures basic logging for user-friendly output.
- Package data access (`add_gene_features`, etc.) now uses `importlib.resources` (with backport for Python < 3.9) for robustness across wheel, sdist, and editable installs.
- Plotting module docstrings and a few source comments cleaned of exaggerated language.
- `set_style` and plot functions now produce more consistent, vector-friendly output by default (pdf/ps fonttype 42, clean spines, etc.).
- README completely rewritten in calm, scientific English with detailed parameter tables, plotting settings, and usage guidance. All hype language removed.
- Test suite significantly expanded (heuristic + permutation, layer remapping, enrichment, plotting headless, error paths, etc.). All tests pass when using local source.

### Fixed
- `pydeseq2` is now properly declared as an optional dependency under the `pseudobulk` extra (was previously only imported dynamically).
- Various minor style and unused-variable issues addressed via ruff.
- Namespace pollution in `scat.*` reduced (internal modules are still importable for advanced use but top-level `dir(scat)` is cleaner).

### Deprecated / Notes
- The previous "review fixes applied" comment in `tl.py` was removed as part of the refactor.
- Old code is preserved in the `backup/` directory for reference.

## [0.7.0] - Previous release (pre-refactor baseline)

See git history and the original README for prior changes. The 0.7 series introduced the dual-track (heuristic/advanced) engine, kb_python layer support, permutation testing, and the current public API surface.

[Unreleased]: https://github.com/scATrans/scatrans/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/scATrans/scatrans/releases/tag/v0.7.0
