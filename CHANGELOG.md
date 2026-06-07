# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- New optional dependency group `[pseudobulk]` for PyDESeq2-based pseudobulk DE.
- `list_available_gene_features()` is now exported at the top level (`scat.list_available_gene_features`).
- Comprehensive parameter documentation in README (active_score, plotting functions, etc.).
- CI workflow (`.github/workflows/ci.yml`) with matrix testing across Python versions and extras.
- Internal module reorganization for better maintainability (see "Changed").

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
