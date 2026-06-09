# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
