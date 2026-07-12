# API stability

This page defines what scATrans treats as **stable public API** versus
**implementation detail**. It applies to the published package on PyPI
and to the documentation on Read the Docs.

scATrans is currently **0.10.x (Beta)** (`Development Status :: 4 - Beta` in
`pyproject.toml`). A future **1.0.0** release will adopt
`Development Status :: 5 - Production/Stable` and treat the contract below
under normal SemVer (breaking changes only with a major version bump).
Until 1.0, minor versions may still refine behaviour with deprecation
warnings where practical.

## Recommended import style

```python
import scatrans as scat

scat.active_score(...)
scat.differential_expression(...)
scat.run_default_pipeline(...)
scat.pl.volcano_plot(...)
scat.qc.unspliced_global(...)
```

Prefer the top-level `scat.*` surface (and `scat.pl` / `scat.qc`) for all
application and paper code.

## Stable public surface

The following are intended to remain importable and callable across
compatible releases (after 1.0: without breaking changes in a minor/patch):

1. **Top-level package** — every name in `scatrans.__all__`, including:
   - scoring / DE / pipeline: `active_score`, `active_score_simple`,
     `differential_expression`, `differential_expression_simple`,
     `run_default_pipeline`, `PipelineResult`, `filter_active_genes`,
     `diagnose_design`, `recommend_workflow`, `WORKFLOW_PRESETS`,
     raw-count helpers (`store_raw_counts`, `ensure_raw_counts`,
     `restore_raw_counts`)
   - gene features: `add_gene_features`, `generate_gene_features_from_gtf`,
     `list_available_gene_features`
   - enrichment: `run_enrichment`, `run_go`, `run_kegg`, `run_gsea`,
     `simplify_enrichment`, `compare_enrichment`, and related helpers listed
     in `__all__`
   - version: `scatrans.__version__`
2. **`scatrans.pl`** — names in `scatrans.pl.__all__` (plotting helpers).
3. **`scatrans.qc`** — names in `scatrans.qc.__all__`.
4. **CLI entry points** declared in packaging metadata (e.g.
   `generate-gene-features` → `scatrans.generate_gene_features:main`).

### `PipelineResult`

`run_default_pipeline` returns a **read-only** `dict` subclass
(`isinstance(result, dict)` is `True`). Guaranteed field keys:

`adata`, `significant`, `all_results`, `candidates`, `enrichment`,
`filter_preset`, `backend`, `meta`.

`meta` always includes `scatrans_version` and `organism`. When
`active_score` ran, it also surfaces the nested `diagnostics` block and
selected run flags from `adata.uns["scatrans"]` (e.g. `use_permutation`,
`gamma_method`, `mode`). The full run metadata remains on
`result.adata.uns["scatrans"]`.

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
`DeprecationWarning` for at least one minor release when practical. Behaviour
changes that affect scientific interpretation should be called out in
`CHANGELOG.md`.

## How to stay safe as a user

1. Depend on `scatrans.__all__` / documented functions, not internal modules.
2. Pin a minor version range in papers and production
   (`scatrans>=0.10.1,<0.11` or `==0.10.1` for exact reproducibility).
3. Record `scatrans.__version__` (and backend versions such as PyDESeq2) in
   Methods / session logs.
