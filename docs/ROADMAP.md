# Roadmap — done vs still open

Working notes from the 0.9.x / 0.10.x post-release review / refactor track
(2026-07). Use this when picking up work later. Not a commitment to
ship dates.

---

## Already done (this track)

Engineering / packaging

- [x] `tl` / `enrich` first-level package split (`scatrans.tl.*`, `scatrans.enrich.*`)
- [x] Public `__all__` / `__dir__` tightened (`pl`, `qc`, `tl`, `enrich`, top-level)
- [x] `PipelineResult` (read-only dict subclass): `__getitem__`, `| =` / `__ior__`,
      `copy` / `deepcopy`, pickle all protocols, `to_dict` / `copy()`
- [x] `show_plot=True` after split (`from .. import pl`)
- [x] Matrix helpers in `_utils` for AnnData `.X` / layers union types
- [x] mypy clean on `src/scatrans --ignore-missing-imports`
- [x] mypy **hard fail** in CI (3.11 core)
- [x] Coverage split: default vs with-slow; `_de` ≥70% gate on default suite
- [x] CI coverage installs optional stacks so reports are not skip-undercounted:
      `pydeseq2`, `memento-de`, `gseapy`, `gtfparse` (same idea as extras
      `pseudobulk` / `memento` / `gsea` / `gene_features`)
- [x] `tests/test_de_core_coverage.py` (scanpy / PyDESeq2 / MixedLM / memento mocks)
- [x] Docs: `docs/api_stability.md` + README / API index links
- [x] Docs fix: `generate-gene-features` not top-level `generate_gene_features_main`
- [x] Lazy-import comments in `tl/pipeline.py`

Out of scope for **this** working tree (release repo is separate)

- [x] Yank PyPI `*.dev*` — done on release side; not required in this dev tree

---

## Still open (prioritized)

### P1 — before calling the package “1.0-ready”

| Item | Why | Notes / acceptance |
|------|-----|--------------------|
| **Real-data smoke** | New abstractions + DE paths need field data | Run `run_default_pipeline` / `active_score` / DE + enrich on 1–2 real h5ad (user-owned data). Record version, backends, pass/fail, odd warnings. |
| **`PipelineResult` “fermentation”** | Review found dunder gaps only under edge tests | Use in notebooks, `joblib`/`ProcessPool`, `copy.deepcopy`, `{**result}` / `|` in real scripts. Prefer keep 0.10.x until no surprises. |
| **Coverage honesty for GSEA / plot** | Default suite hides `@pytest.mark.slow` / `plot` | **Partly done:** CI coverage job installs gseapy/gtfparse; with-slow XML is uploaded. Still optional: badges or docs that show **two** Codecov numbers (default vs with-slow). Confirm plot suite coverage separately if needed. |
| **Raise non-DE “science-adjacent” coverage** | Still thin vs 1.0 narrative | Targets (suggest, not CI-hard yet): `_velocity` ≥65%, `pp_bias` / gene-features ≥60%, overall default suite ≥70–75%. |

### P2 — product / packaging

| Item | Why | Notes |
|------|-----|--------|
| **Wheel size / KEGG bundling** | ~9 MB wheel; KEGG commercial terms | Optional: KEGG not in default wheel; download or extra; GO compress / slim |
| **Dual mouse gene-feature tables** | Confusing defaults | Doc default table + when to switch; or single default + `list_available_gene_features` |
| **CITATION.cff / authors** | Academic citability | Real names/ORCID when ready to publish papers |
| **Optional `plot` extra** | Slim core install | Move heavy viz deps if desired |
| **Further module split** | `tl.active` / `pl` still large | Only if maintainability hurts; keep public imports stable |

### P3 — 1.0 release checklist (do not rush)

- [ ] Freeze public API list = `scatrans.__all__` + `pl.__all__` + `qc.__all__` (already documented in `api_stability.md`)
- [ ] `pyproject.toml` classifier → `Development Status :: 5 - Production/Stable` **only with 1.0.0**
- [ ] SemVer policy in docs (breaking = major after 1.0)
- [ ] CHANGELOG entry for 1.0 with “stable surface” summary
- [ ] Real-data smoke signed off; `_de` gate still green with pydeseq2 present
- [ ] Optional: TestPyPI dry-run of release workflow

**Recommended versioning until then:** stay on **0.10.x** (Beta; current **0.10.1**). Do not label 1.0 until P1 real-data smoke + no open PipelineResult semantics changes.

### P4 — nice-to-have / earlier review notes

- Decision-tree doc: which DE backend when
- Output column dictionary more prominent from README
- Small downloadable demo h5ad (Zenodo) for tutorials
- “How scATrans relates to scanpy / scVelo / gseapy / PyDESeq2”
- `set_nature_style` → `set_publication_style` alias (trademark caution)
- CONTRIBUTING.md / issue templates
- `py.typed` if shipping type information to dependents
- pre-commit mypy flags == CI mypy flags (already close)

---

## How to re-verify quickly

```bash
# Types
mypy src/scatrans --ignore-missing-imports

# Default suite (matches daily CI tests)
pytest -m "not plot and not slow" -q

# _de coverage gate (needs pydeseq2 for fair measure)
pip install "pydeseq2>=0.4.0"
pytest -m "not plot and not slow" --cov=scatrans --cov-report=xml:coverage-default.xml -q
# parse line-rate for *_de.py ≥ 70%

# With slow (GSEA etc.)
pytest -m "not plot" -q
```

GitHub review zip (dev tree):

```bash
python scripts/make_release_zips.py
# or single: scatrans-github-upload-YYYY-MM-DD.zip
```

---

## Context for “continue later”

1. User may bring **real h5ad** for smoke — treat as P1 acceptance, not more refactor.
2. Prefer **not** expanding public API surface before 1.0.
3. If CI `_de` gate fails again, first check **pydeseq2 is installed** on the coverage job (see `.github/workflows/ci.yml`).
4. Release / PyPI publish lives in a **separate** repo from this workspace unless explicitly merged.

Last updated: 2026-07-08 (after P0/P1 + coverage-gate fix).
