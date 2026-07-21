# API Stability

Stable public API versus implementation detail for the PyPI package and the
published documentation.

scATrans is **0.10.x (Beta)** (`Development Status :: 4 - Beta` in
`pyproject.toml`). A future **1.0.0** release will move to Production/Stable and
apply SemVer (breaking changes only with a major version bump). Until 1.0, minor
versions may refine behavior with deprecation warnings where practical.

**Scientific heuristic defaults** (e.g. `HEURISTIC_FILTER_DEFAULTS` values
such as `logfc_cutoff`, residual/FDR gates) are
**not frozen API**. They may change in a minor release when domain
feedback warrants it; the public guarantee is the **parameter names and
filter semantics**, not the numeric defaults. Always report the installed
`scatrans.__version__` and the cutoffs you used (or
`filter_active_genes(preset=...)`) when publishing.

## Recommended import style

```python
import scatrans as scat

scat.active_score(...)
scat.differential_expression(...)
scat.run_default_pipeline(...)
scat.pl.volcano_plot(...)
scat.qc.unspliced_global(...)
scat.qc.regime_diagnosis(...)
```

Prefer the top-level `scat.*` surface (and `scat.pl` / `scat.qc`) for all
application and paper code.

## Stable public surface

The following are intended to remain importable and callable across
compatible releases (after 1.0: without breaking changes in a minor/patch):

1. **Top-level package** — every name in `scatrans.__all__`, including:
   - primary workflow: `partition_de_by_mechanism` (DE selects → mechanism
     partition; the recommended entry point) and its result `PartitionResult`
     (fields `adata`, `regime`, `gene_table`, `selected`, `programs`,
     `enrichment`, `meta`, parallel to `PipelineResult`). The composite
     `run_default_pipeline(select_by="composite")` path is deprecated as a
     discovery entry (see {doc}`faq`).
   - scoring / DE / pipeline: `active_score`, `active_score_simple`,
     `adaptive_active_score`, `add_adaptive_score`, `adaptive_weight`,
     `labeling_anchor`, `add_abundance_normalized_residual`,
     `annotate_mechanism_class`, `program_mechanism`,
     `nascent_activity_score` (active-transcription detection score; opt-in
     detection columns via `partition_de_by_mechanism(add_nascent_score=True)`,
     decoupled from the mechanism partition),
     `threshold_sensitivity`,
     `differential_expression`, `differential_expression_simple`,
     `run_default_pipeline`, `PipelineResult`, `filter_active_genes`,
     `diagnose_design`, `recommend_workflow`, `WORKFLOW_PRESETS`,
     raw-count helpers (`store_raw_counts`, `restore_raw_counts`;
     `ensure_raw_counts` is **deprecated** — use
     `store_raw_counts(..., mode="auto")`)
   - gene features: `add_gene_features`, `generate_gene_features_from_gtf`,
     `list_available_gene_features`
   - enrichment: `run_enrichment`, `run_go`, `run_kegg`, `run_gsea`,
     `simplify_enrichment`, `compare_enrichment`, and related helpers listed
     in `__all__`
   - version: `scatrans.__version__`

**Scientific maturity (not the same as import stability):** differential
expression, enrichment, and plotting that do **not** depend on
spliced/unspliced layers are suitable for routine use. Nascent-transcription
scoring (`active_score` and its velocity-dependent add-ons) remains
**experimental**; see {doc}`faq` and the README.
2. **`scatrans.pl`** — names in `scatrans.pl.__all__` (plotting helpers).
3. **`scatrans.qc`** — names in `scatrans.qc.__all__`
   (`unspliced_global`, `regime_diagnosis`).
4. **CLI entry points** declared in packaging metadata (e.g.
   `generate-gene-features` → `scatrans.generate_gene_features:main`).

### `PartitionResult`

`partition_de_by_mechanism` returns a dataclass with fields:

`adata`, `regime`, `gene_table`, `selected`, `programs`, `enrichment`, `meta`.

- **`selected` / `gene_table`:** DE membership is only in `selected`; mechanism
  columns live on both when annotation ran. Detection columns
  (`nascent_poisson_z`, `de_reproducible`, …) appear only if
  `add_nascent_score=True`.
- **`regime`:** copy of the reliability pre-flight dict (also under
  `meta["regime"]`).
- **`meta` keys:** always `scatrans_version`, `organism`, `de_source`, `de`,
  `select`, `regime`, `mechanism`, `programs`, `nascent_score`
  (`enabled` / `status` / …). Mechanism is **always** residual-based;
  `nascent_score` never drives `transcription_support` / program directions.

### `PipelineResult`

`run_default_pipeline` returns a **read-only** `dict` subclass
(`isinstance(result, dict)` is `True`). Guaranteed field keys:

`adata`, `significant`, `all_results`, `candidates`, `enrichment`,
`filter_preset`, `backend`, `meta`.

`meta` always includes `scatrans_version`, `organism`, and `select_by`. When
`active_score` ran, it also surfaces the nested `diagnostics` block and
selected run flags from `adata.uns["scatrans"]` (e.g. `use_permutation`,
`gamma_method`, `mode`). On velocity-capable objects the pipeline also
records **`meta["regime"]`** from `scat.qc.regime_diagnosis` (fail-soft if
layers are missing). Optional add-ons record under `meta["bias"]`,
`meta["adaptive"]`, and `meta["mechanism"]` when used.
The full run metadata remains on `result.adata.uns["scatrans"]`.

In-place mutation (`result[k] = …`, `result |= …`, `update` / `pop` / …)
raises `TypeError`. Use `result.to_dict()` or `result.copy()` for a mutable
plain `dict`. Attribute access (`result.candidates`) mirrors the same keys.

## Supported but not path-stable

These imports work and are useful for advanced users, but **import paths
below the package root may move** in a minor release before 1.0 (and only
with a major bump after 1.0):

| Import | Status |
|--------|--------|
| `from scatrans.tl import active_score` | Supported re-export of the public function |
| `from scatrans.enrich import run_enrichment` | Supported re-export |
| `import scatrans.tl.active` | **Implementation module** — not a stability promise |
| `import scatrans.enrich.ora` | **Implementation module** — not a stability promise |
| `scatrans._de`, `scatrans._utils`, `scatrans._*` | **Private** — may change without notice |

If you only need a public callable, import it from `scatrans` or from the
package root (`scatrans.tl` / `scatrans.enrich`), not from leaf modules such
as `scatrans.tl.active`.

## Private / unstable

- Any name starting with `_` (except what is re-exported in a public
  `__all__` by mistake — treat `_` as private).
- Contents of `scatrans/data/` gene-set **files** (bundled data may be updated
  between releases; licensing is separate — see the license page).
- Undocumented keyword arguments marked experimental in docstrings.

## Deprecations

Before removing or renaming a stable symbol after 1.0, scATrans will emit a
`DeprecationWarning` for at least one minor release when practical. Behavior
changes that affect scientific interpretation should be called out in
`CHANGELOG.md`.

## How to stay safe as a user

1. Depend on `scatrans.__all__` / documented functions, not internal modules.
2. Pin a minor version range in papers and production
   (`scatrans>=0.10.7,<0.11` or `==0.10.7` for exact reproducibility).
3. Record `scatrans.__version__` (and backend versions such as PyDESeq2) in
   Methods / session logs.
