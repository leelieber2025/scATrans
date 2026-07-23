# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Unreleased

### Added
- `program_mechanism_induction_matched`: program-level mechanism tests that
  control for induction strength (OLS `support ~ logFC + membership`, optional
  nearest-logFC matching). Exposed as
  `partition_de_by_mechanism(induction_matched=True)` →
  `PartitionResult.programs_induction_matched`.
- `annotate_mechanism_class(..., flag_induction_confound=True)`: marks high-
  induction stabilization calls in `induction_confounded` and down-weights
  `mechanism_confidence` (does not relabel or change program tests). Penalty
  shapes: `graded` (default) or `smooth`.
- `annotate_mechanism_class(preset="high_precision")`: sets
  `class_threshold=1.0` (lower hard-call rate; explicit `class_threshold`
  overrides). Wired as `partition_de_by_mechanism(mechanism_preset=...)`.
- `annotate_mechanism_class(suppress_hard_labels_when_unreliable=True)`: when
  regime `reliability` is below `min_reliability_for_hard_labels` (default
  0.05), hard per-gene classes become `ambiguous` (support unchanged).

### Changed
- `PartitionResult.summary()`: program-level counts first; documents that
  per-gene classes are soft (`per_gene_labels_are_soft`).
- `partition_de_by_mechanism`: warns when `sample_col` is missing
  (`meta["pseudoreplication_warning"]`).
- `run_enrichment(..., allow_mechanism_class_ora=False)`: warns if the query
  table carries a `mechanism_class` column (ORA on mechanism subsets is
  discouraged). Set `allow_mechanism_class_ora=True` to silence.

### Fixed
- `program_mechanism_induction_matched`: drop `zip(..., strict=True)` so the
  function runs on Python 3.9 (CI matrix).
- `pl.enrich_barplot`: prefer non-null `Description`, else `Term` (avoids
  literal `"nan"` labels on bundled GO/KEGG tables).
- `pl.compare_dotplot`: rotate long or many cluster x-labels to reduce overlap.

### Documentation
- `program_mechanism`: docstring notes that arbitrary KEGG/GO screens are
  confounded by gene length; prefer mechanism-coherent sets.
- User guide / FAQ / API reference aligned with induction-matched programs,
  reliability hard-label suppression, and enrichment guards (this release).

## [0.10.7] - 2026-07-20

### Changed
- `active_score` / `active_score_simple` outputs now sort DE-first (`p_adj`, then
  `logFC`) instead of by the composite score. Display-order only — no value or
  metric changes. The composite `active_score` column is legacy; prefer
  `partition_de_by_mechanism` (or `filter_active_genes(select_by="de")`).
- `pl.volcano_plot(style="auto")` now colors and sizes points by the nascent
  excess residual by default (`color_by="unspliced_excess_residual"`, legacy
  `velocity_residual` accepted) instead of the composite `active_score`. Pass
  `color_by="active_score"` to restore the legacy composite coloring. Pure-DE
  tables (no residual column) fall back silently to up/down/ns categories.

### Documentation
- Removed the legacy 0–100 composite `active_score` score from the docs
  (API reference, quickstart, user guide, FAQ, statistical guidance, domain
  assumptions, method page). Its derived columns/cutoffs
  (`active_score_pval` / `active_score_fdr` / `active_score_cutoff` /
  `active_score_fdr_cutoff`) and "rank/color by `active_score`" guidance are no
  longer documented; docs now direct selection/ranking to DE (`p_adj`, `logFC`)
  and the nascent residual (`unspliced_excess_residual` / `unspliced_excess_fdr`).
  The `active_score()` function itself is unchanged and still documented. The
  method page renumbers the permutation equation (old Eq. 6 → Eq. 5) after
  dropping the composite-score equation. No code changes: the composite column
  is still produced by `active_score()`.
- Tutorials index: **Start here** learning path, comparison of the three
  partition notebooks, approximate runtimes/extras, and data-file table.
- `references.md`: GSE226488 / Derbois et al. 2023 citation, tutorial-object
  notes, and how large `.h5ad` files are distributed.
- User guide workflow: **Input data and layers** (what AnnData needs, layer
  names, raw counts, regime pre-flight).
- Tutorial notebooks: shorter intros (less composite-ranking / migration
  prose); standalone DE GSEA note points to the changelog instead of a long
  bug narrative; GSE226488 **Reproduce** section no longer depends on
  machine-local paths.

### Added
- **`scatrans.tl.nascent_activity_score(...)` and
  `partition_de_by_mechanism(add_nascent_score=True)`** — optional active-transcription
  **detection** columns, kept **decoupled** from the mechanism partition.
  - Score: pseudobulk (target vs reference) variance-stabilized Poisson-z of the
    nascent increase (`nascent_poisson_z`), without length/intron residualization
    (unlike the default `unspliced_excess_residual` used for mechanism).
  - Also: `dlog_unspliced` / `dlog_spliced`, and spliced-side DE-reproducibility
    flags `de_reproducible` / `de_repro_frac` (annotation only; do not gate
    membership). The flag cannot distinguish a DE false positive from a genuine
    stabilization-driven gene (both can show near-zero nascent excess).
  - **Decoupling:** the Poisson-z is an absolute, induction-coupled nascent
    increase. It is **never** used for transcription-vs-stabilization labels or
    program-level mechanism calls (those remain residual-based). Using it as
    `residual_col` would mis-label highly induced stabilization targets as
    transcription-driven.
  - `add_nascent_score=True` only **appends** columns; it does not change
    `transcription_support`, `mechanism_class`, or program directions. Fail-soft
    via `meta["nascent_score"]`. Layers auto-resolve (`spliced`/`unspliced` or
    `mature`/`nascent`). Exported top-level and from `tl`. Tests:
    `tests/test_nascent.py`.

## [0.10.6] - 2026-07-20

### Added
- **`scatrans.partition_de_by_mechanism(...)`** — primary DE → mechanism workflow.
  (1) a pluggable **DE step selects** changed genes; (2) scATrans **partitions**
  them into transcription-driven vs stabilization-driven classes using the nascent
  residual. `de=` accepts `"builtin"`, a `de_method` name / kwargs dict routed to
  `differential_expression`, a precomputed DataFrame (`de_logfc_col` /
  `de_padj_col`), or a callable `adata -> DataFrame`. Always runs a reliability
  pre-flight (`regime_diagnosis`, scales `mechanism_confidence`), soft per-gene
  annotation (never gates membership), and optional program-level table when
  `gene_sets=` is provided. Returns `PartitionResult` (`adata`, `regime`,
  `gene_table`, `selected`, `programs`, `enrichment`, `meta`). Down-regulation is
  not yet mechanism-resolved (`unclassified_down`). Exported top-level and from
  `tl`. Tests: `tests/test_partition.py`.
- **`scatrans.qc.regime_diagnosis`** — maps global unspliced fraction to a
  dataset-level `reliability` in [0, 1] (U-shaped: high in a normal band, lower
  at low- or high-unspliced extremes) plus `regime` (`ok` / `low_unspliced` /
  `high_unspliced`) and a message. `run_default_pipeline` records `meta["regime"]`
  (fail-soft) and, when `annotate_mechanism=True`, scales mechanism confidence by
  `reliability`. This is a data-quality / gamma check only; dynamic-vs-steady-state
  detection is not implemented.
- **Mechanism annotation module `scatrans.tl.mechanism`** (annotation-only; never
  gates membership):
  - `annotate_mechanism_class(...)` — `transcription_support`, `mechanism_class`,
    `mechanism_confidence` (scaled by dataset `reliability`). Per-gene labels are
    intentionally low-confidence; prefer program-level pooling.
  - `program_mechanism(results, gene_sets, ...)` — threshold-free program-level
    inference (competitive Mann–Whitney + BH FDR on pooled support).
  - `threshold_sensitivity(results, ...)` — DE list size and Jaccard vs reference
    across a (padj, logFC) grid.
  - `run_default_pipeline(..., annotate_mechanism=True)` — fail-soft add-on;
    diagnostics in `meta["mechanism"]` (default off).
- **`select_by="de"`** on `filter_active_genes` and `run_default_pipeline` —
  membership from DE gates only (`p_adj` / `logFC` + direction; defaults
  `padj<0.05` and `|log2FC|>1` when cutoffs omitted); nascent/composite gates are
  skipped (columns remain as annotations). Sort by `p_adj` then `logFC`. Default
  `select_by="composite"` preserves prior behavior; mode recorded in
  `meta["select_by"]`. Incompatible with `preset="significant"`.
- **First-class pipeline options** `bias_method=` and `adaptive_weighting=` /
  `adaptive_anchor=` on `run_default_pipeline` (additive, fail-soft, default off):
  - `bias_method="abundance"` / `"abundance_length"` →
    `unspliced_excess_residual_abnorm` (`meta["bias"]`).
  - `adaptive_weighting=True` → `adaptive_score` / `adaptive_score_pct`
    (`meta["adaptive"]`).
- **Configurable reliability anchor** for adaptive weighting (`anchor=` on
  `add_adaptive_score` / `adaptive_active_score`): `"de"` (default) or a callable /
  array / Series. Helper `scat.labeling_anchor(column="new_log2fc", threshold=1.0)`
  anchors reliability on a metabolic-labeling truth column when available.
  Diagnostics include `anchor` and `n_anchor_induced` (`n_anchor_de_induced`
  retained as a back-compat alias).

### Deprecated
- **`run_default_pipeline(select_by="composite")`** as a gene-discovery path —
  composite ranking mixes DE and nascent legs and is not recommended for production
  discovery. Still available (default unchanged) but emits `DeprecationWarning`;
  prefer `partition_de_by_mechanism(...)` or
  `run_default_pipeline(..., select_by="de", annotate_mechanism=True)`.

### Fixed
- **External DE path of `partition_de_by_mechanism`**: reported `logFC` / `p_adj`
  / `p_val` and mechanism direction now come from the *selecting* DE (not the
  builtin `active_score` DE). Missing external raw p clears `p_val` to NaN so a
  stale builtin p does not sit next to external stats.
- **`run_default_pipeline`**: invalid `select_by` raises; add-on columns
  (mechanism / bias / adaptive) are re-attached to `candidates` after annotate
  (with an explicit copy); regime pre-flight failure uses cautious
  `reliability=0.5` instead of implying full confidence.
- **`filter_active_genes`**: residual cutoff of `-inf` truly disables the residual
  gate (no longer drops genes with missing residual).
- **`annotate_mechanism_class`**: missing/NaN `logFC` is `unknown`, not
  `unclassified_down`.
- **`program_mechanism`**: empty generic background falls back to competitive
  background with a warning; equal mean support is `ns` (not a directional call).
- **Enrichment**: `expand_enrichment_genes` accepts GSEA `leading_edge` /
  `Lead_genes` and warns when no gene-list column is present; Entrez-style
  numeric `Series` indexes are treated as gene IDs (aligned with DataFrame path).
- **Adaptive DE anchor** uses strict `logFC > 1` to match `filter_active_genes`
  / partition DE gates.

## [0.10.5] - 2026-07-18

### Added
- **Reliability-adaptive weighting of the nascent leg** (`scatrans.tl.adaptive`,
  exported as `scat.adaptive_active_score`, `scat.add_adaptive_score`,
  `scat.adaptive_weight`). An **additive** wrapper (core `active_score` is
  unchanged) that estimates how informative the unspliced-excess leg is on the
  data at hand — the AUC of `unspliced_excess_residual` for recovering the
  obvious DE-induced genes (`logFC >= 1` & `p_adj < 0.05`) — and produces an
  `adaptive_score` whose nascent leg is weighted by that reliability:
  `w = clip(k*(reliability-0.5), 0, w_max)` (defaults `k=4`, `w_max=2`). The
  proxy is shrunk to 0 when it anti-correlates with induction (e.g. steady-state
  / late velocity snapshots) and up-weighted (>1) when it is highly reliable
  (e.g. metabolic-labeling data). Returns diagnostics (`reliability_auc`,
  `w_proxy`, `verdict`). Heuristic ranking, not a calibrated FDR; the DE anchor
  is isolated in `_de_induced_anchor` for future swapping to a labeling anchor.
- **Abundance-/length-normalized unspliced-excess residual**
  (`scatrans.tl.bias`, exported as `scat.add_abundance_normalized_residual`).
  Adds `unspliced_excess_residual_abnorm`. `method="abundance"` (default) uses a
  scale-free excess `delta / (total_us_counts + quantile(total_us_counts,
  floor_quantile))` (default 0.75) that removes the abundance / nuclear-retention
  artifact — nuclear-retained lncRNAs (e.g. *MALAT1*) and very high-abundance
  genes no longer dominate the top of the ranking. `method="abundance_length"`
  additionally applies a gentle robust residualization on `log1p(gene_length)`
  and `log1p(intron_number)` to suppress long-gene artifacts. Improves
  interpretability of the residual ranking; it does not change the residual's
  reliability on steady-state snapshots (a kinetic, not a bias, limitation).

## [0.10.4] - 2026-07-15

### Added
- **Raw-count sidecar snapshot (survives HVG / cell subsetting)**:
  `store_raw_counts(..., sidecar=True)` (the new default) now also writes a
  label-indexed snapshot of the full obs × var count matrix to
  `adata.uns['scatrans']['raw_snapshot']`. Because `uns` is not tied to the
  obs/var axes, the snapshot survives HVG gene subsetting, cell subsetting,
  `copy()`, and `write_h5ad()` — unlike `layers['counts']` (trimmed on both
  axes) and `adata.raw` (trimmed on cells). It is aligned by cell/gene **name**,
  so it also tolerates reordering.
  - `restore_raw_counts(..., full_genes=True)` reconstructs the complete
    pre-HVG gene universe as a new AnnData (for DE / enrichment on all genes),
    even when called on an HVG-subsetted or cell-subsetted object.
  - `restore_raw_counts(..., prefer_snapshot=True)` (default) restores counts
    aligned to the current cells/genes, absorbing cell subsetting and gene
    reordering.
  - `store_raw_counts(sidecar='ondisk', snapshot_path=...)` writes the full
    matrix to a standalone `.h5ad` and keeps only a lightweight pointer in
    `uns`, avoiding a doubled count matrix in memory / in the main file for
    large datasets.
  - **Velocity layers preserved too**: `spliced`/`unspliced` (or
    `mature`/`nascent`) are captured into the snapshot and restored alongside
    the counts, so active-transcription analysis can run on the full gene set
    after HVG subsetting via `restore_raw_counts(full_genes=True)`.

### Changed
- **`store_raw_counts` gained `mode="force"|"auto"`**: `"auto"` is idempotent
  and recovery-aware (reuse an existing integer counts layer, or recover from
  `adata.raw` when `.X` is already normalized/log-transformed). This folds the
  former `ensure_raw_counts` logic into one function.
- **Density preserved through store → restore**: restoring no longer forces a
  dense matrix to sparse (or vice versa); the original matrix format round-trips.

### Fixed
- **`restore_raw_counts` no longer misreads a log-normalized `adata.raw` as
  counts**: the `adata.raw` fallback is used only when it actually looks like
  integer counts. Following the common scanpy convention where `adata.raw`
  holds log-normalized data, restoring it into `.X` as raw counts would have
  been a silent correctness bug; it now raises a clear error pointing to the
  sidecar snapshot instead.
- **Snapshot restore with duplicate cell/gene names**: label-aligned restore now
  fails with an actionable message (call `obs_names_make_unique()` /
  `var_names_make_unique()`) instead of pandas' cryptic `InvalidIndexError` when
  barcodes are duplicated (e.g. multiple batches before de-duplication).
- **Dead velocity `raw_*` layers removed**: `store_raw_counts` no longer writes
  `raw_spliced`/`raw_unspliced`/`raw_mature`/`raw_nascent` layers. They were
  position-aligned (so trimmed by HVG/cell subsetting) and never read back —
  the sidecar snapshot now preserves velocity correctly instead.

### Deprecated
- **`ensure_raw_counts()`** — use `store_raw_counts(..., mode="auto")`. The old
  name remains as a thin alias and emits a `DeprecationWarning`.
- **`store_raw_counts(save_raw=True)`** — `adata.raw` is commonly reserved for
  log-normalized data, so writing raw integer counts there is ambiguous. The
  sidecar snapshot supersedes it; `save_raw=True` now emits a
  `DeprecationWarning`.

## [0.10.3] - 2026-07-14

### Added
- **Gene-level UpSet (`scat.pl`)**: three new functions for comparing gene
  overlap across multiple DE results / gene lists — the gene-level companion to
  the term-level `enrich_upsetplot`.
  - `build_gene_membership(de_results, ...)` — tidies a `{name: de_df}` mapping
    (or ready-made `{name: [gene, ...]}` lists) into a gene × set 0/1 membership
    matrix. `direction="separate"` (default) splits each DE result into
    `name::up` / `name::down` sets so common-up and common-down genes are both
    visible; `"up"`/`"down"`/`"both"` give one set per result. DataFrame inputs
    are filtered via `filter_active_genes` (`pval_cutoff` / `padj_cutoff` /
    `logfc_cutoff`); per-set gene lists are stored in `membership.attrs["gene_sets"]`.
  - `gene_upsetplot(...)` — draws the UpSet (pure matplotlib, no external
    `upsetplot` dependency) from either a `{name: de_df}` mapping or a pre-built
    membership matrix. Fully recolorable (`set_color`, `intersection_color`,
    `dot_color`, `inactive_color`, `line_color`; `intersection_color`/`dot_color`
    also accept a per-column list to highlight specific intersections).
  - `common_genes(membership, direction="up"|"down", ...)` — pulls the
    intersection genes (strict, or relaxed via `min_sets=`, or an explicit
    `sets=` subset) back out as a list ready for `run_enrichment`.

### Changed
- **Plotting color customization**: `enrich_upsetplot` gained the same color
  parameters as `gene_upsetplot`; `bias_diagnostic_plot` gained
  `raw_color`/`corrected_color`/`trend_color`; `gamma_shrinkage_plot` gained
  `cmap`/`color`. All previously hardcoded their colors. Defaults are unchanged,
  so existing figures look identical.

### Fixed
- **Package layout / CI**: GitHub trees must not ship pre-split flat modules
  `src/scatrans/tl.py` or `src/scatrans/enrich.py` beside the `tl/` and
  `enrich/` packages (Python imports the package; the flat files are dead
  code and trip `tests/test_package_layout.py`). Release packaging refuses
  to include them and fails if they reappear in the source tree.

### Changed
- **Versioning (single source of truth)**: package version lives only in
  `src/scatrans/_version.py` (`__version__`). `pyproject.toml` reads it via
  `[tool.setuptools.dynamic]`; runtime import, Sphinx `release`, and
  `scripts/make_release_zips.py` all consume the same value. Removed
  `setuptools_scm` / `SETUPTOOLS_SCM_PRETEND_VERSION` and the multi-place
  `0.10.2` fallbacks. **To bump a release, edit `_version.py` (and
  CHANGELOG), then rebuild** — the release helper also syncs `CITATION.cff`
  and `packaging/ecosystem-packages/meta.yaml`.
- **Documentation (plotting)**: `docs/user_guide/plotting.md`,
  `docs/user_guide/enrichment.md`, and `docs/api/index.md` now document the
  full `scat.pl` surface aligned with current `pl.py`, including
  `compare_dotplot`, `context=` / notebook vs paper defaults, `label_repel`,
  `positive_logfc_only`, UpSet/Venn helpers, `gamma_shrinkage_plot`, and
  batch export helpers. American English spelling in public plot docstrings.
- **CI PyPI publish is manual only**: `.github/workflows/publish.yml` no longer
  runs on GitHub Release `published`. Creating a release/tag will not upload to
  PyPI. Use Actions → "Publish to PyPI" → Run workflow, or publish dist
  artifacts via your separate path. Build version is taken from `_version.py`.
- **Documentation / citations**: expanded `docs/references.md` with method papers
  (GSEA, Enrichr, DESeq2, RNA velocity, BH FDR, Phipson–Smyth permutation *p*,
  Huber, clusterProfiler, GSEApy paper) and verified external links; fixed
  broken GO license URL and README Read-the-Docs markdown; domain-assumptions
  README link points at in-repo `docs/domain_assumptions.md`.

### Fixed
- **`active_score` direction consistency**: the significance leg
  (`-log10(p_adj)`) is gated by upregulation (`logFC > 0`, or
  `mixedlm_coef > 0` under MixedLM so the gate matches what `p_adj` tests).
  Previously a strongly downregulated gene could still score ~50/100 from
  the directionless p-value term alone. Observed and permutation paths
  share `_composite_active_score_terms`.
- **`lambda_pval` scale**: estimated on direction-positive genes only, so
  extreme downregulated p-values no longer inflate λ and shrink s3 for
  true up genes.
- **`filter_active_genes` residual direction**: residual magnitude cutoffs
  now follow `logfc_direction` (`up` → residual > c; `down` → residual < -c;
  `both` → residual sign concordant with logFC). One-sided
  `unspliced_excess_fdr` is skipped for `down`/`both` (with a warning).
- **`padj_cutoff` alias**: `filter_active_genes` and `extract_gene_lists`
  accept `padj_cutoff` (preferred); `pval_cutoff` remains as a legacy name
  for the adjusted-p filter.
- **GSEA ranking metric**: `active_score` and other non-negative columns are
  no longer auto-selected for preranked GSEA (requires signed ranks for
  bidirectional NES). Auto-pick prefers `logFC` / t-stat-like columns; using
  a one-sided metric emits an explicit warning. Guarded by
  `tests/test_package_layout.py`.
- **Package layout guard**: `scatrans/tl.py` and `scatrans/enrich.py` must not
  reappear beside the `tl/` and `enrich/` packages (unreachable dead code after
  the package split). `MANIFEST.in` excludes them; layout tests enforce
  import paths resolve to package `__init__.py`.
- **Huber `valid_feat`**: require `gene_length > 0` (not `>= 0`). Length 0 is a
  missing-annotation sentinel; `log1p(0)=0` is an extreme x-leverage point that
  biased Huber slopes and all residuals. `intron_number >= 0` unchanged (0
  introns is valid). Warn when zero/missing lengths are present; GTF feature
  tables leave missing length as NaN instead of filling 0.
- **MixedLM significant / filter**: require `mixedlm_coef > 0` (or `< 0` for
  `logfc_direction="down"`) so p_adj direction matches the tested coefficient;
  sign-discordant genes no longer enter the built-in significant list. Discordance
  is logged at **warning** level.
- **`extract_gene_lists` p/logFC resolution**: warn when falling back from
  adjusted p to raw p (cutoff would otherwise silently inflate false positives);
  recognize Seurat `avg_log2FC` / `avg_logFC`; warn when no logFC column is found
  instead of silently returning empty lists.
- **`active_score` design diagnosis**: no longer discards `diagnose_design(...)`
  return value. Warnings (e.g. &lt;3 samples/group) are forwarded via
  `logger.warning` and stored under
  `adata.uns["scatrans"]["diagnostics"]["design"]`. Failures are logged at
  warning level instead of silent debug.
- **`run_gsea` mapping-rate check**: same 20% gate as ORA (`_check_gene_set_mapping_rate`).
  Warns with input vs gene-set symbol examples when overlap is low; returns empty
  with `reason="no_ranked_genes_mapped"` at zero overlap (avoids opaque
  gseapy "No gene sets passed" RuntimeWarnings from case/ID mismatches).
  Mapping stats stored under `attrs["gsea_info"]` / `gene_set_info["mapping"]`.
- **GTF gene features (`transcript_id`)**: missing/empty `transcript_id` on exons
  raises a clear `ValueError` (no bare KeyError). Non-positive `gene_length` from
  empty exon unions is stored as **NaN** (not 0). Simple-path
  `_maybe_add_gene_features` fills sparse tables (&lt;50% usable length) from the
  bundle without overwriting existing length&gt;0 — closes the "1% real / 99%
  missing never auto-completes" trap with partial GTF feature tables.
- **Regression tests for #1 / #4**: Huber fit must not let length=0 genes bias
  the slope; downregulated genes must not outrank mild-up on active_score when
  residual weight is 0.
- **Docs / FAQ**: top-N `active_score` is not DE-gated (nascent excess can rank
  alone); GSEA mapping/`gene_case` troubleshooting; design warnings under
  `diagnostics["design"]`; MixedLM coef gate on significant list.
- **CI statistical guards** (`tests/test_statistical_guards.py`, ~3s):
  (1) null-label permutation Type I @0.05 in [0.03, 0.07];
  (2) planted up/down ground-truth ranking (up in top-N, down not);
  (3) Huber `n_genes_used_for_fit` excludes length 0/NaN; active_score
  diagnostics match.
- **Domain assumptions doc** (`docs/domain_assumptions.md`): explicit list of
  product/domain conventions (active = upregulation-oriented; s2 independent
  of DE; `pval_cutoff` → adjusted p; GSEA signed ranks; Huber length &gt; 0;
  etc.) so semantic intent is not left implicit.
- **`active_score` scale is within-run only**: soft-scale λ =
  `median(positive)/ln(2)` from the run's gene vectors. Documented as
  high-risk if used for cross-dataset / cross-subset / HVG comparisons;
  lambdas + note stored in `diagnostics["scoring"]`; info log on each run;
  guard test in `test_statistical_guards`.
- **Domain-assumption verification suite**
  (`tests/test_domain_assumptions_verified.py`): locks documented semantics
  (λ data-dependence, s1/s3 gate, s2 DE-independence, MixedLM coef gate,
  Huber length filter, padj vs raw-p, GSEA signed ranks, empty significant
  without permutation).

## [0.10.2] - 2026-07-11

Patch release over PyPI `0.10.1` (CI / Python 3.9 compatibility).

### Fixed
- **Python 3.9 import**: `pl.py` was missing `from __future__ import annotations`,
  so PEP 604 annotations such as `str | None` raised `TypeError` at import
  time on 3.9 (CI matrix). Fixed; other `src/scatrans` modules already had
  the future import.
- **Python 3.9 `zip(strict=...)`**: `volcano_plot(style="ggvolcano")` used
  `zip(..., strict=True)` (3.10+ only), causing
  `TypeError: zip() takes no keyword arguments` on 3.9. Replaced with plain
  `zip` (fixed-length category tuples). Full `src/scatrans` scan found no other
  3.10+ constructs (`match`/`case`, `bit_count`, etc.).
- **PyDESeq2 `DeseqStats(n_cpus=...)`**: some installed pydeseq2 versions accept
  `n_cpus` only on `DeseqDataSet`, not `DeseqStats`. Init kwargs are now filtered
  by `inspect.signature(cls)` so optional args (`n_cpus`, `quiet`, …) are dropped
  when unsupported instead of raising `TypeError` (and without mypy
  ``[misc]`` on ``cls.__init__``).
- **CI base matrix (no optional extras)**: PyDESeq2 replicate-count gate and
  Memento non-integer `counts=` check now raise **before** importing the
  optional package (so `.[dev]`-only CI jobs get `ValueError` not
  `ImportError`). GSEA regression that needs `gseapy` is
  `@pytest.mark.skipif` when the extra is missing.

## [0.10.1] - 2026-07-11

### Added
- **API stability doc** (`docs/api_stability.md`): public surface vs
  implementation modules (`scatrans.tl.active`, …); linked from README and
  API reference.
- **Roadmap** (`docs/ROADMAP.md`): done vs open work (1.0 prep, coverage,
  packaging) for continuing after the 0.9.x / 0.10.x refactor track.
- **DE core tests** (`tests/test_de_core_coverage.py`): scanpy / PyDESeq2 /
  MixedLM / memento branches for default-suite coverage of `scatrans._de`.
- **High-risk combo regression tests** (`tests/test_combo_regressions.py`):
  `mode="advanced"` + permutation, `allow_advanced_pseudobulk` gate,
  Memento consistent/inconsistent null metadata, `ranking_mode="nascent_excess"`
  residual-only score (including DE-poisoned MixedLM null independence).
- **Sample / RE-cluster label shuffle** for MixedLM permutation nulls
  (`_shuffle_condition_labels` in `_permutation.py`).

### Changed
- **CI**: mypy on `src/scatrans` is a **hard fail** (Python 3.11 core job);
  coverage reports split into `coverage-default.xml` (no slow/plot) and
  `coverage-with-slow.xml`; default suite must keep `scatrans._de` line
  coverage ≥ 70%. Coverage steps install optional science deps
  (`pydeseq2`, `memento-de`, `gseapy`, `gtfparse`) so Codecov is not
  under-counted for DE / GSEA / gene-features when the matrix cell is only
  `.[dev]`.
- **`tl` / `enrich` package layout**: implementation split into submodules
  (`scatrans.tl.active`, `scatrans.enrich.ora`, …). Public imports
  (`import scatrans as scat`, `from scatrans.tl import active_score`) are
  unchanged.
- **`run_default_pipeline` returns `PipelineResult`**: a **read-only**
  `dict` subclass (so `isinstance(result, dict)` is still `True`) with
  attribute access (`result.candidates`). Item/attribute assignment and
  `update`/`pop`/`clear` raise `TypeError`; use `result.to_dict()` for a
  mutable shallow copy. New optional key `meta` always has
  `scatrans_version` / `organism`, and now also surfaces
  `diagnostics` plus selected run flags from `adata.uns["scatrans"]`
  (`use_permutation`, `gamma_method`, `mode`, …) when `active_score` ran.
- **Public `__all__` tightened** on `scatrans.pl`, `scatrans.qc`,
  `scatrans.tl`, and `scatrans.enrich` (private helpers no longer listed;
  import them from implementation submodules if needed).
- **Top-level API**: `generate_gene_features_main` is no longer re-exported
  from `scatrans` (CLI entry point `generate-gene-features` unchanged).

### Fixed
- **`use_mixed_model=True` + `use_permutation=True`**: with
  `perm_de_backend='same'` (default), each permutation now refits the same
  MixedLM (`condition + (1|sample)`) used for the observed DE so
  `active_score_pval` / `active_score_fdr` compare like with like. Previously
  the null always used scanpy `t-test_overestim_var` while logging that it
  used the “same” backend — invalid active_score FDR under sample structure.
  **Also:** MixedLM nulls shuffle labels at the **sample / RE-cluster** level
  (not independently per cell). Cell-level shuffle splits a biological sample
  across conditions and breaks hierarchical exchangeability; unpaired nulls
  reassign conditions to whole observed RE clusters and pin those cluster IDs
  as `sample_col` so `n_groups` cannot collapse. Paired designs
  (`paired_replicates=True`) shuffle within subject. Non-MixedLM backends keep
  cell-level shuffle. `perm_de_backend='fast'` still uses t-test (with warning);
  `unspliced_excess_fdr` remains DE-backend-independent.
  **Scope:** DE-null / shuffle-unit mismatches only affect composite
  `active_score_pval`/`fdr` when `weight_fc`/`weight_pval` are non-zero
  (default `ranking_mode="composite"`). Under
  `ranking_mode="nascent_excess"` those weights are forced to 0, so
  active_score FDR tracks residual only.
- **`extract_gene_lists`**: accepts `padj` / `p.adjust` / `log2FoldChange` etc.;
  warns when no p-value column exists (avoids silent empty lists).
- **`gseaplot` fallback RES**: sorts the ranked list high→low before the
  running ES walk (unsorted `all_results['logFC']` no longer draws a wrong curve).
- **`run_gsea` duplicate gene IDs**: keep max |score| per gene (not first row).
- **`active_score_simple` / `differential_expression_simple`**: always isolate
  a working copy before mutation so `copy_input=False` cannot write gene
  features into the caller's `.var`.
- **`perm_de_backend='fast'` on pseudobulk counts**: scanpy path now
  `normalize_log1p`s integer pb `.X` before `rank_genes_groups` (avoids
  expm1 overflow / inf logFC under `pb_use_total_for_x` + PyDESeq2 main).
- **Scanpy DE non-finite fill**: ±inf / NaN logFC and p-values neutral-filled
  after `rank_genes_groups` (defense if raw counts still reach scanpy).
- **Memento counts safety**: refuses non-integer `counts=` / `layers['counts']`
  (uses `_resolve_aligned_raw_counts` + integer check); no silent log-scale DE.
- **`volcano_plot` style='auto'**: labels respect `label_by` (default `p_adj`)
  instead of always preferring `active_score`.
- **`enrich_vennplot` 4-set**: fixed `NameError` (`groups` → `clusters`) that
  crashed any 4-cluster Venn after the region-legend addition.
- **`active_score` metadata**: `sample_col` is recorded under
  `uns["scatrans"]` for pseudobulk runs (not only MixedLM), matching
  `differential_expression`.
- **MixedLM + Memento**: mutually exclusive — raising instead of silently
  preferring MixedLM and ignoring Memento.
- **`paired_replicates=True` without shared sample IDs**: warning that
  pairing has no effect.
- **`ranking_mode='nascent_excess'`**: always forces residual-only weights
  `(0, 0, 1)`, overriding custom `weight_*` with a warning (no silent
  composite score under a residual-only mode name).
- **`gamma_method='raw'`**: zero-expression reference genes use the global
  U/S sum ratio instead of `eps/eps ≈ 1`.
- **PathwayDenester**: hypergeometric tail via `scipy.stats.hypergeom.sf`
  (exact combinatorial sum for large GO terms no longer hangs on `n!`).
- **`enrich_barplot`**: real horizontal barplot (not a dotplot alias).
- **4-set `enrich_vennplot`**: side legend lists all exclusive region sizes
  (pairwise/triple/all), not only on-circle uniques.
- **`copy_input=False`**: always isolates a working copy before mutation so
  the caller's AnnData is never modified (active_score / differential_expression).
- **`store_raw_counts` after HVG**: preserves prior full universe as
  `raw_gene_list_full` (sticky); ORA prefers it for enrichment background.
- **`robust_median` gamma**: median anchor now excludes zero-expression
  genes in the reference (same mask as EB); sparse data no longer pulls
  `base_gamma` toward 1 via `(eps/eps)`.
- **Size-factor U/S**: zero-total cells keep factor=1 (not `target/1e-8`).
- **`extract_gene_lists` / `_expand_gene_list_input`**: RangeIndex tables
  without a gene column no longer treat `logFC` as gene IDs (return empty +
  warning).
- **`diagnose_design`**: missing/wrong `sample_col` is a warning (not “no
  sample_col provided”); `workflow_preset="pseudobulk_report"` only when
  samples/group ≥ mixed-model minimum (avoids recommending permutation
  while also warning against it).
- **gseapy Enrichr organism**: map `mouse`→`Mouse` / `human`→`Human`;
  organism-less fallback logs a loud warning.
- **GSEA**: Entrez-style numeric index accepted; log when both `logFC`
  and `active_score` exist and default rank is logFC.
- **`expand_enrichment_genes`**: split genes on `;` or `,`.
- **`copy_input=False`**: docstring + runtime warning when the caller's
  object is mutated in place.
- **MixedLM diagnostics export**: `logFC_method` and
  `n_genes_logFC_mixedlm_sign_discordant` from `de_df.attrs` are now copied
  into `adata.uns["scatrans"]["diagnostics"]["mixed_model"]` for both
  `active_score` and `differential_expression` (previously only logged).
- **Permutation DE matrix matches observed**: each permutation re-resolves
  `layers['counts']` and forwards `min_counts_per_gene=pydeseq2_min_counts`
  into `_run_de_wrapper`, so null DE uses the same count source as the main
  PyDESeq2/Memento path (not post-aggregation `.X` = U+S). Fixes invalid
  `active_score_fdr` under recommended `store_raw_counts` + pseudobulk.
- **`subset_col` / `subset_values` label normalize**: numeric /
  CSV-style `"1.0"` labels match the same way as `groupby` contrasts.
- **`de_preprocess='normalize_log1p'` on PyDESeq2 path**: skipped (with
  warning) when `skip_auto=True` so counts are not log-transformed before
  DESeq2.
- **`uns['log1p']` with max>20**: trust the marker for non-integer data
  (high-depth / bulk-like log matrices) instead of clearing it and risking
  double `log1p`; still reject when the matrix looks raw + library-size
  dominated (stale marker).
- **DE label Categorical defense**: MixedLM / Memento / perm inject labels
  via `_as_contrast_categorical` (float `1.0` ↔ `"1"`).
- **`filter_active_genes` unknown preset message**: lists `significant`
  and aliases.
- **`strict_pydeseq2_counts` gap for aggregated `layers['counts']`**:
  `_pseudobulk_with_layers` now records pre-aggregation count-likeness for
  every aggregated layer (`pb_layer_is_count_like`, plus
  `pb_counts_is_count_like` for the `counts` layer). PyDESeq2 and
  `_resolve_aligned_raw_counts` honor that flag so a non-count matrix stored
  under `layers['counts']` cannot pass the safety net after `np.round()`
  (previously only `.X` / `pb_x_is_count_like` was protected).
- **MixedLM `logFC` is sample-aware**: effect size is the scanpy-style log2FC
  of *mean of per-sample means* within each condition (equal weight per
  biological sample), not a cell-weighted group mean that a large sample can
  dominate. `mixedlm_coef` remains the LMM fixed-effect coefficient that
  `p_val` tests; `attrs['logFC_method']` /
  `n_genes_logFC_mixedlm_sign_discordant` document the dual statistics.
- **`active_score(..., show_plot=True)`** after the `tl` package split used
  `from . import pl` (resolved to missing `scatrans.tl.pl`) and silently
  skipped plotting. Now imports `scatrans.pl` correctly.
- **`PipelineResult` read-only semantics** (dict-subclass edge cases):
  - `copy.copy` / `copy.deepcopy` / pickle via `__copy__` / `__deepcopy__` /
    `__reduce__` (all pickle protocols; no `__getstate__`/`__setstate__` —
    once `__reduce__` is defined those hooks are never called)
  - **`result |= {...}`** no longer silently mutates via C-level
    `dict.__ior__` (raises `TypeError`)
  - **`result | other`** returns a **new mutable** `dict` (left operand
    unchanged)
  - **`result.copy()`** matches **`to_dict()`** (mutable plain `dict`);
    `copy.copy(result)` still returns a read-only `PipelineResult`
- **`run_default_pipeline` `meta` now includes diagnostics**: previously the
  pipeline read `adata.uns["scatrans"]` only for a log branch, then overwrote
  a same-named `meta` variable with just version/organism. `result.meta`
  now merges the nested `diagnostics` block and selected run flags so the
  docstring / CHANGELOG promise matches runtime.
- **Memento DE `.attrs` after reindex**: `_run_memento_de` now re-assigns
  `n_genes_not_returned_by_memento` / `n_genes_missing_pval` after
  `reindex`/`fillna` instead of relying on experimental pandas
  `DataFrame.attrs` propagation (which can silently drop diagnostics).
- **`filter_active_genes` numeric validation**: residual / FDR /
  `effective_gamma_*` / `delta_variance_min` cutoffs now use the same
  `_coerce_numeric_cutoff` path as `active_score_cutoff` / `pval_cutoff` /
  `logfc_cutoff` (clear `ValueError` on non-numeric input; `None` still
  allowed for optional FDR / optional bounds).
- **`_compute_perm_velocity_delta`**: removed no-op if/else that called the
  same estimator on both branches; both tracks still share
  `_compute_velocity_delta` (layers encode the track).
- **`active_score` mode branch**: defensive `else: raise AssertionError` if
  an unhandled `mode` reaches velocity estimation (avoids silent all-NaN
  `gamma_ref` if a new mode is added without wiring).
- **`gseaplot` fallback RES**: removed incorrect min-max scaling of the
  approximate running enrichment score to [0, 1], which destroyed sign and
  contradicted negative NES annotations. Fallback now plots the signed
  weighted KS curve (and warns when precomputed `gsea_details` are missing).
- **`enrich_dotplot` multi-cluster term selection**: no longer applies a
  global positional `.head()` after per-cluster sampling that could drop
  entire later clusters; keeps a per-cluster union as documented.
- **`enrich_vennplot`**: full exclusive-region counts for 3-set diagrams;
  4-set warns that only exclusive counts are labeled (prefer upset).
- **`bias_diagnostic_plot`**: document and annotate that the left-panel
  trend is length-only 1D, while the real Huber correction uses
  length + intron number.
- **Cross-group enrichment p-adjust**: `run_go(..., adjust_across_all=True)`
  and `compare_enrichment(..., adjust_across_clusters=True)` now re-adjust
  with the same `p_adjust_method` as the sub-calls (via kwargs; default
  `fdr_bh`), instead of always hardcoding Benjamini–Hochberg.
- **PyDESeq2 uses `counts=` / `layers['counts']`**: count matrix for DESeq2
  no longer always comes from `.X`. Aligned `counts=` is preferred; when
  `use_pseudobulk` + `pydeseq2` and a `counts` layer exists, pseudobulk
  aggregates that layer into `.X` so log1p `.X` is not mistaken for counts.
- **MixedLM BH NaN safety**: non-finite Wald/LRT p-values are neutral-filled
  *before* `multipletests`; a single NaN no longer collapses all `p_adj` to 1.
- **Built-in `significant` list scale**: uses
  `PSEUDOBULK_FILTER_DEFAULTS` (residual > 0.05, active_score ≥ 5, …) when
  `is_pseudobulk=True`, matching `filter_active_genes(preset="pseudobulk")`.
- **Huber residual for unannotated genes**: after a successful fit, genes
  missing length/intron features get median-centered delta (not residual=0).
- **`velocity_delta_layer` metadata**: cell-level path correctly reports
  `size_factor_normalized_spliced_unspliced` (always normalized).
- **Default mouse gene features**: prefer `mouse_2020A_...` over
  lexicographic `Mus_musculus.GRCm39...`.
- **ORA term-count diagnostics**: `n_terms_tested` /
  `n_terms_size_excluded` replace the misnamed “filtered” counter; legacy
  key kept as alias of tested count.
- **Permutation failure vectors**: failed shuffle replicates return length
  `n_vars` (not `n_obs`).
- **MixedLM false “failed fits”**: no longer discard genes solely because
  statsmodels reports `converged=False` (common when random-effect variance
  is on the boundary / singular `cov_re`). Results are kept when the target
  condition coefficient and Wald p-value are finite; LRT falls back to Wald
  when log-likelihoods are non-finite. Fixes flaky reverse-contrast tests on
  small paired designs.
- **Permutation method note**: metadata now states that γ/DE/residual are
  recomputed under shuffled labels (layers fixed), not that gamma is frozen.
- **`.uns["scatrans"]` merge**: non-sticky keys from prior runs are cleared;
  `None` means “feature off” and removes the key (no sticky `sample_col` etc.).
- **Group label `"1.0"` strings**: `_normalize_group_label` maps CSV-style
  integer floats to `"1"` / `"2"` like true floats.
- **Advanced precomputed Mu/Ms**: after size-factor rewriting U/S, precomputed
  moments are discarded and recomputed (with warning).
- **Memento `counts=` alignment**: matrices go through the same shape/name
  coercion as PyDESeq2; misaligned counts raise clearly.
- **MixedLM `logFC`**: sample-aware (mean of per-sample means) scanpy-style
  log2FC for cutoffs; LMM coefficient stored as `mixedlm_coef` (see Unreleased
  Fixed for the cell-weighted → sample-aware correction).
- **Filter presets**: `effective_gamma_min/max` default to `None` (no silent
  drop when `show_effective_gamma=True`).
- **`differential_expression`**: default `pb_use_total_for_x=False` (do not
  silently DE on U+S total when velocity layers exist).
- **Enrichment `_apply_p_adjust`**: NaN-safe (finite-only multipletests).
- **Empty `significant` UX**: clearer warnings pointing to `filter_active_genes`.
- **Docs**: `docs/api/index.md` no longer claims `generate_gene_features_main`
  is importable from the top-level `scatrans` package.
- **mypy clean** under `mypy src/scatrans --ignore-missing-imports` (union
  narrowing for AnnData/matrices, `Path` vs `str`, heterogeneous
  `load_info` dicts, pl size-legend numerics, etc.).
- **Matrix union-type helpers** in `_utils` (`_dense_expression_matrix`,
  `_matrix_copy`, `_matrix_shape`, `_matrix_sum_axis0/1`,
  `_matrix_row_subset_sum_axis0`, `_as_var_dataframe`) used from `_de`,
  `_velocity`, `tl.active`, and pseudobulk paths so AnnData ``.X`` /
  ``.layers`` wide unions no longer need per-line mypy ignore. Also
  replace ``matplotlib.cm.viridis`` with ``plt.get_cmap("viridis")``.

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