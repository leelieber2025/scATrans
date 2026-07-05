# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## 2026-07-05

### Fixed
- **`run_gsea` made downstream `DataFrame` operations pathologically slow**: the
  per-term running-enrichment-score curves stored in `res_df.attrs["gsea_details"]`
  (tens of millions of floats for a genome-wide ranked list against thousands of
  gene sets) were deep-copied by pandas' `NDFrame.__finalize__` on essentially
  every subsequent `DataFrame` operation (`.head()`, `.copy()`, slicing, column
  assignment), including inside `scat.pl.enrich_dotplot` and `scat.pl.gseaplot`.
  This made those two plotting functions appear to hang indefinitely (in practice:
  several minutes) on real GSEA results, and made even a plain `gsea_res.head()`
  slow. `gsea_details` and `ranking` are now wrapped in a small dict subclass
  whose `__deepcopy__` is an O(1) identity return (safe: these payloads are
  write-once/read-only after `run_gsea` returns). `enrich_dotplot`/`gseaplot`
  on a real ~3700-term result went from indefinite/minutes to <0.1s / ~6s.
- **`run_gsea` returned several numeric columns (`ES`, `NES`, `pvalue`, `p.adjust`,
  `FWER_pval`, `Tag_percent`, `Gene_percent`) as `object` dtype** instead of
  `float`, inherited as-is from `gseapy`'s `res2d`. Now coerced explicitly via
  `pd.to_numeric`. Contributed to (but was not the root cause of) the hang above.

## 2026-07-04

### Added
- **`scat.pl.volcano_plot(style=...)`**: ggVolcano-inspired styles from
  [BioSenior/ggVolcano](https://github.com/BioSenior/ggVolcano) — ``style="ggvolcano"``
  (teal/grey/orange Up-Down-Normal, theme_bw, FDR labels) and ``style="gradual"``
  (gradient by ``-log10 FDR``). Default ``style="auto"`` keeps the previous look.
  Documented in README §3.4 and API reference.

### Fixed
- **PyDESeq2 pseudobulk tests**: use ``pb_x_layer="counts"`` + ``pb_use_total_for_x=False``
  (``spliced+unspliced`` sums are non-integer and correctly fail ``strict_pydeseq2_counts``).
- **MixedLM edge-case test**: design meets ≥4 samples/group after stricter mixed-model gates.
- **Small-sample edge tests** (`tests/test_small_sample_edges.py`): 1–2 cells/group, all-zero genes.
- **Regression tests**: ``max_avoid_points`` volcano subsampling, numpy ``raw_gene_list`` enrichment
  universe, and ``recommend_workflow`` auto-disabling ``use_permutation`` on small pseudobulk designs.
- **README**: mixed-model small-sample guidance (≥4 samples/group); documented ``paired_replicates``
  and ``filter_active_genes(preset='significant')``.
- **`run_enrichment` / `run_kegg` DataFrame `gene_list`**: gene symbols read from index
  (or ``gene`` / ``names`` columns), not column names.
- **`run_gsea` DataFrame `ranked_genes`**: index-based ``all_results`` support.
- **MixedLM**: NaN neutral-fill on degenerate fits; composite ``condition::sample``
  random-effect groups when replicate labels are reused across conditions (``paired_replicates=True``
  for paired designs).
- **CI pytest failure on matrix legs without `pydeseq2`**: two new pseudobulk regression tests
  (`test_strict_pydeseq2_counts_rejects_log_normalized_*`) now use the same
  `@pytest.mark.skipif(importlib.util.find_spec("pydeseq2") is None, ...)` guard as
  `tests/test_de_backends.py`, instead of relying on production-code check ordering.
  Restored PyDESeq2 `ImportError` (missing dependency) ahead of data-validation `ValueError`
  in `_run_de_wrapper`, so users without `pydeseq2` see the install hint first.
- **Pseudobulk `strict_pydeseq2_counts` check ran on already-rounded data**: `_pseudobulk_with_layers`
  always rounds aggregated sums to integers, so the PyDESeq2 count-likeness check in
  `_run_de_wrapper` (which inspected the rounded pseudobulk `.X`) could never detect that the
  underlying per-cell source (e.g. a `pb_x_layer` pointing at log-normalized data) was not raw
  counts. The check now runs on the pre-aggregation, pre-rounding source matrix inside
  `_pseudobulk_with_layers` and the verdict is carried via `adata.uns["pb_x_is_count_like"]` to
  `_run_de_wrapper`, so `strict_pydeseq2_counts=True` correctly rejects non-count pseudobulk input.
- **`de_preprocess="auto"` false "already log-normalized" detection regressed by real cell
  heterogeneity**: the per-cell library-size CV heuristic added to detect raw counts when
  `uns['log1p']` is missing (e.g. after `anndata.concat()`, which drops `.uns`) could misfire on
  correctly log-normalized data from heterogeneous cell populations (different cell types/states
  varying in expression breadth inflate per-cell CV even post-log1p), triggering a spurious
  second `normalize_total` + `log1p` pass that systematically compresses effect sizes. Added a
  cross-gene mean-variance dispersion check (`_x_gene_dispersion_looks_raw`), which is not
  confounded by cell-level heterogeneity, and now require it to corroborate the library-size CV
  signal before concluding data still needs normalization.

### Documentation
- README: noted that `anndata.concat()` drops `uns['log1p']` and recommends re-setting the marker
  or passing `de_preprocess="none"` explicitly after combining pre-normalized samples.

## 2026-07-03

### Changed (user-visible default behavior)
- **`active_score()` built-in `significant` list** now uses the same default thresholds as
  `filter_active_genes(preset="heuristic")`: `logFC > 0.35`, `unspliced_excess_residual > 1.0`,
  `active_score >= 55`, `active_score_fdr < 0.25`, plus `p_adj` / `unspliced_excess_fdr` cutoffs.
  Re-running analyses with defaults may return a **different** significant gene set than the
  2026-06-20 release (typically more genes when signal is present; aligns with the documented
  post-hoc workflow).
- **`active_score(logfc_cutoff=...)` default** lowered from `0.5` to `0.35` to match the heuristic
  preset (only affects the built-in significant mask and metadata; ranking in `all_results` unchanged).
- **`gamma_method="empirical_bayes"` + `prior_weight`**: fixed a numeric floor so `prior_weight` in
  the usual tuning range (0.5–5.0, including the default 5.0) now changes `count_pseudocount` and
  shrinkage weights instead of being pinned at 1.0. EB gamma values may differ from the 2026-06-20
  release at the same `prior_weight`.

### Fixed
- **`scat.pl.gseaplot`**: RES curve / hit ticks and the bottom ranked-metric bar now share the same
  gene order by preferring `gsea_result.attrs["ranking"]` over caller-supplied row order.
- **`scat.pl` robustness**: placeholder figures for empty or incomplete volcano/comet inputs;
  `bias_diagnostic_plot` validates external `axes`; `active_genes_heatmap` default `show=True` no
  longer hits a Python `UnboundLocalError` on `plt`.
- **PyDESeq2 pseudobulk DE**: log warning + `de_df.attrs` diagnostics when genes are skipped by
  `min_counts_per_gene` or marked NaN by DESeq2 independent filtering (neutral fill is not
  "tested and non-significant").
- **Memento DE**: log warning + `n_genes_not_returned_by_memento` when memento drops genes
  internally; reindexed neutral rows are documented in diagnostics.
- **`restore_raw_counts`**: reject `layers["counts"]` restore when `raw_gene_list` order differs
  from current `var_names` (same guard as `adata.raw`).
- **`active_score` pseudobulk**: `pb_x_layer="X"` now means `adata.X`, matching
  `differential_expression()`.
- **`_velocity._estimate_eb_prior_from_reference`**: removed unreachable trim fallback branch;
  trim-skipped path labeled `empirical_bayes_median_mad`.
- **`_utils._get_group_mean`**: removed redundant sparse/dense branches.

### Added
- Shared `HEURISTIC_FILTER_DEFAULTS` in `tl.py` (single source for heuristic preset + built-in
  significant mask).
- Tests: gseaplot ranking alignment, `prior_weight` EB sensitivity, significant vs heuristic
  filter parity, PyDESeq2/Memento diagnostic paths, `restore_raw_counts` order guard,
  `active_genes_heatmap` default show, `pb_x_layer="X"` sentinel.

### Documentation
- README built-in `significant` section updated to match `HEURISTIC_FILTER_DEFAULTS`.

## 2026-07-02

### Fixed
- **compare_enrichment / concat_compare_results metadata**: `attrs["clusters"]` now lists only clusters that contributed rows to the combined table (skipped, failed, and empty per-cluster runs remain in `per_cluster` diagnostics only). Fixes ghost cluster names in both the main compare API and the concat wrapper.
- **DE schema validation**: `_run_de_wrapper` validates that all backends return `logFC`, `p_val`, `p_adj` with at least one finite value per column.
- **MixedLM robustness**: near-constant genes use variance threshold; non-converged fits and missing condition coefficients are counted as failed fits instead of guessing the second coefficient; `failed_fit_rate` added to result attrs with percentage in warning log.
- **DE warnings**: selective `_de_warning_context()` suppresses deprecation/future noise only (scanpy/PyDESeq2/MixedLM/BH), no longer blanket `ignore` on all warning categories.
- **Memento audit column**: native adjusted p-values preserved as `memento_p_adj_native` when returned by memento-de; package `p_adj` still uses BH for cross-backend consistency.
- **DE edge-case tests**: new `tests/test_de_edge_cases.py` (schema validation, constant-gene MixedLM failure counting, scanpy direction agreement, all-zero gene finiteness).
- **Diagnostics / docs**: `failed_fit_rate` exposed in `active_score` mixed-model diagnostics and `differential_expression` metadata; README adds a short DE backend decision guide.
- **copy_input performance**: `active_score` and `differential_expression` now honor `copy_input=False` — combined obs filters (subset + target/reference) perform at most one `AnnData.copy()` when `copy_input=True`, and zero when `copy_input=False` and no obs filtering is required (previously an unconditional `[keep_mask].copy()` always ran).
- **filter_active_genes permissive mode**: default/permissive thresholds for `pval_cutoff`, `active_score_fdr_cutoff`, and `unspliced_excess_fdr_cutoff` now use `float("inf")` instead of `1.0`, so genes with adjusted p-value or permutation FDR exactly equal to 1.0 are no longer silently dropped (strict `<` vs `1.0` bug).
- **CI lint**: removed unused `import scanpy as sc` from `tl.py` and `_permutation.py` (leftover from permutation refactor); fixed `test_enrich_api.py` formatting.
- **CI tests**: `test_pp_bias_cli` GTF generator tests now `pytest.importorskip("gtfparse")` so base installs without `scatrans[gene_features]` skip instead of failing.

### Added
- Regression test `test_filter_active_genes_permissive_keeps_padj_one`.

## 2026-06-27

### Added / Improved
- Clarified and documented the `gamma_method="empirical_bayes"` implementation as **hierarchical (分层) gamma estimation** for the reference U/S ratio (README keeps the CN term; source now English-only).
- Stronger emphasis on "always pass explicit target_group/reference_group".

### Changed / Fixed (critical)
- **Bug 1**: Eliminated dead code — `tl.py` now calls the canonical `run_permutation_test(...)` in `_permutation.py` instead of duplicating the Parallel loop. Removed ~duplicated logic and the maintenance trap. `valid_expr` is passed explicitly for consistent behavior.
- **Bug 2**: Fixed double normalization (double log1p) in permutations. When `de_preprocess="normalize_log1p"` (or auto that applies), the value passed to permutation tasks is forced to "none" so that perm copies of the already-transformed adata are not re-normalized. Prevents systematically biased FDR.
- **Scientific Error 1**: Fixed EB gamma: `sigma2` in `_apply_empirical_bayes_gamma` now correctly includes the `n_ref` factor: `1.0 / (n_ref * U_r + c) + ...`. This was causing n_r-fold over-estimate of observation noise and excessive shrinkage (especially bad for the small-ref case EB is meant to help). `n_ref` is computed from r_mask in the caller.
- **Scientific Error 2**: Clarified `robust_median` docs in code/README: it is a heuristic variant of `heuristic_shrink` (different base_gamma estimator) and is **not** Bayesian/EB/hierarchical. Renamed descriptions to prevent user confusion.
- **Design fixes**:
  - Tightened Memento raw counts check from `>=` to exact `== n_vars` + `var_names` equality (prevents misaligned HVG/raw usage).
  - Renamed shadowing local `ad = ...` in `active_score_simple` / `differential_expression_simple` (avoids hiding `import anndata as ad`).
  - Made `min_cells_per_sample` private (`_min_cells_per_sample`) + doc note (was public but did nothing).
  - Removed all Chinese comments ("分层") from source (tl.py, _velocity.py); retained in README only.
  - Replaced risky bare `except TypeError` for PyDESeq2 design_factors/design compat with explicit `_pydeseq2_uses_design_factors()` using `importlib.metadata` version check.

- Reduced anndata category storage log noise during internal DE/perm (narrow logger level bump).
- All tests + targeted EB/perm/explicit-norm cases pass after fixes.

### Additional fixes (round 2)
- **Scientific Error 3**: `generate_gene_features_from_gtf()` now computes **true exon union length per gene** (merge overlapping intervals) instead of naively summing all exons across all transcripts. Prevents ~N_transcript-fold overestimation of gene_length for multi-isoform genes (affects only users who generate their own tables; bundled .parquet files are unaffected).
- **Scientific Error 4**: `logFC` is now normalized toward consistent **log2 scale** across backends inside `_run_de_wrapper`:
  - wilcoxon, mixedlm, memento: divided by ln(2)
  - t-test / PyDESeq2: left as native log2
  - Updated docstring. Makes `logfc_cutoff` semantically comparable.
- **Bug 3**: `_pseudobulk_with_layers` now uses a **per-run UUID-based private separator** (instead of fragile "||") to build internal keys. Completely eliminates split errors when sample names contain "||".
- **Bug 4**: Removed redundant `import warnings as _w` inside `with warnings.catch_warnings()` in `_fit_huber_bias_correction`.
- **Design 6**: `pval_cutoff` deprecation warning now fires **on every use** of the legacy name (previously only when != 0.05).
- **Design 7**: `run_go(ontology="CC"/"MF")` now emits a clear INFO log explaining that only BP is bundled and the others require gseapy + network.
- **Design 8**: Made the "max > 50" heuristic in `_prepare_log_normalized_expression` (mixed model prep) more robust: checks for negatives, uses lower threshold (20), better warning message.
- **Design 9**: `valid_expr` is now **explicitly passed** from `tl.py` into `run_permutation_test` (no more hidden reliance on adata.var after the fact).

### Post-review hardening
- **Fragility fix (high)**: Added explicit schema validation + safe column access + clear warning in `_run_memento_de` for expected 'de_coef'/'de_pval' from memento.binary_test_1d. Prevents silent fallback (or crash) if upstream memento-de changes column names/structure. Updated optional dep pin to "memento-de>=0.1.0,<0.3.0".
- **Cleanup + consistency (high)**: `run_permutation_test` now owns the small-space FDR decision (use_fdr=False + disabled_reason="small_permutation_space" when n_perm < 100). Removed vestigial dead `if not perm_use_fdr` and useless max_perm check/assign in tl.py. `disabled_reason` (and `perm_disabled_reason`) now properly returned and stored in adata.uns["scatrans"] metadata.
- **Diagnostics improvement**: MixedLM per-gene fits now count genes hitting the neutral-except fallback (`n_genes_failed_fit`). Recorded in diagnostics["mixed_model"] (active_score) and DE metadata (differential_expression). Warning emitted when >0 so users are not unaware of silent neutral (delta_var=0, p=1) genes.
- Only confirmed high-impact issues from the review list were addressed; medium/nuance items (e.g. EB small-ref prior, advanced fallback diag completeness) were already mitigated by existing diagnostics/fallbacks or not correctness bugs.

All new issues addressed. Full test suite green.

## 2026-06-20

### Added / Changed
- `filter_active_genes` now accepts `logfc_direction="up"|"down"|"both"` (default remains `"up"` for backward compatibility and "active" semantics).
  - `logfc_cutoff` is interpreted as a positive magnitude in all modes.
  - `"down"`: selects logFC < -cutoff (downregulated genes from differential_expression results).
  - `"both"`: selects |logFC| > cutoff.
  - Sorting for pure-DE tables is now direction-aware (most-negative-first for down, largest |logFC| for both).
  - This directly supports the common request for downregulated candidates in standalone DE workflows.
- Updated docstrings, tests, and examples.

## 2026-06-19

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

## 2026-06-18

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

## 2026-06-14

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
- README and docstrings extensively updated with manuscript-export examples, `run_go`, provenance details, and `adjust_across_all` guidance.
- Full test coverage for new paths (per-ontology attrs, within_ontology p.adjust, save+tsv+dir creation, expand with Ontology, dual-cutoff warning, etc.). All tests pass.

## 2026-06-09

### Added
- **Rich runtime diagnostics**: `active_score` now automatically computes global unspliced fraction (via integrated `qc.unspliced_global`), captures detailed bias-correction fit results (coefficients, n_genes_used, fallback status), stores per-gene `effective_gamma`, and records everything under `adata.uns["scatrans"]["diagnostics"]`. A concise run summary is logged at completion.
- New public helper `scat.pl.velocity_phase_portraits(adata, genes, groupby=...)` for quick visual inspection of U vs S relationships on top hits.
- `examples/real_data_template.py` — heavily commented, non-runnable but copy-adaptable template demonstrating the full recommended real-data workflow, QC, diagnostics inspection, and publication plotting.
- Explicit documentation of the permutation approximation (velocity layers fixed from original data; only labels shuffled) in code metadata, logs, and the new "Choosing mode" section of the README.
- `effective_gamma` column is now added to `adata.var` (and included in result tables when present) for transparency of the gamma used in the delta calculation.
- Mixed model support (`use_mixed_model=True` + `sample_col`): statsmodels LMM (~ condition + (1 | sample)) as cell-level replicate-aware DE backend (lightweight Python analogue to dreamlet/variancePartition/dreampy LMM + Libra mixed options). Replaces scanpy/pydeseq2 p/logFC when enabled.
- Delta Variance: `delta_variance` (condition-attributable variance fraction, 0-1, variancePartition-style) and `delta_var_pval` (LRT) computed during mixed fit; always surfaced in `all_results` / adata.var / diagnostics when mixed used.
- New `active_score` options: `use_delta_variance_pval` (bool) + `delta_var_pval_cutoff` to optionally include the delta LRT pval as supplementary filter in significant gene selection.
- Full backward compat; mixed path documented vs. pseudobulk; guidance + references to Libra, dreampy, NEBULA, dreamlet in README.
- Tests + direct verification cover mixed, delta col, filter option, incompatibility with pseudobulk, and no regression on legacy paths. All tests pass.
- Enhanced `README.md` with real-data workflow pointers and stronger emphasis on inspecting diagnostics.
- `scat.pl` now documents the new phase-portrait helper; all existing `ax=` support preserved.
- Internal: `_fit_huber_bias_correction` and velocity helpers now return extra fit/quality information (used for the new diagnostics) while preserving full backward compatibility of public results.

### Changed / Improved
- Major usability & paper-readiness upgrade to diagnostics, metadata, and user guidance.
- `qc.unspliced_global` is now called automatically inside `active_score` (result stored); the function remains directly usable as a pre-flight check (`scat.qc.unspliced_global(adata)`).
- Added prominent "Choosing `mode`: heuristic vs advanced (and common pitfalls)" section + decision guide to README.
- **Major internal refactor** (2026-06-09): Core logic in `tl.py` (`active_score`) was extracted into private supporting modules (`_utils.py`, `_de.py`, `_velocity.py`, `_permutation.py`, plus bias correction in `pp_bias.py`). The public `active_score` function is now a thin, readable orchestrator while preserving 100% identical behavior, return values, and side effects on `adata`.
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

## Earlier

Pre-refactor baseline introduced the dual-track (heuristic/advanced) engine, kb_python layer support, permutation testing, and the current public API surface. See git history and the original README for prior changes.