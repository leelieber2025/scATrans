# Roadmap

**Sequence:** testing → manuscript → **1.0** and preprint. Stay on **0.10.x
(Beta)** until the manuscript is complete. This is not a calendar commitment.

**Scope:** a method package whose primary workflow is DE selection followed by
mechanism partition (transcription-driven versus stabilization-driven via
`partition_de_by_mechanism` / `PartitionResult`) from spliced/unspliced or
nascent/mature layers—not a multi-omics platform. The lower-level nascent
scoring engine (`active_score`) remains available but is **not** the primary
gene-discovery path.

**Non-goals (until well after 1.0):** LLM/MCP productization, bundling large
numbers of third-party methods, or replacing general DE tools such as scanpy.

---

## Phase map

```text
  Phase T          Phase M           Phase R
  Testing          Manuscript        Release
  (0.10.x)         (0.10.x freeze)   (1.0.0 + preprint)
      │                 │                    │
      ▼                 ▼                    ▼
  CI green +       Paper text,         SemVer stable,
  paper-critical   figures, methods    PyPI 1.0,
  paths covered    locked to a         CITATION + preprint
  + real-data      pinable 0.10.x
  smoke for paper
```

| Phase | Version | Exit criterion |
|-------|---------|----------------|
| **T — Testing** | 0.10.x | Paper-critical paths tested; real-data smoke for manuscript figures recorded; no open API churn planned |
| **M — Manuscript** | **freeze** a 0.10.x pin | Full draft + figures/tables reproducible from that pin; Methods cite version + cutoffs |
| **R — Release** | **1.0.0** | Paper complete → tag 1.0 → preprint (+ Zenodo/software cite as needed) |

Do **not** label 1.0 during drafting. Do **not** expand the public API surface
during Phases T–M unless the paper absolutely requires it.

---

## Phase T — Testing (current priority)

Goal: enough automated + real-data evidence that **Methods and Results are
defensible**, not maximal coverage for its own sake.

### T1 — Paper-critical automated tests (must)

| Area | Why for the paper | Acceptance |
|------|-------------------|------------|
| **Active score core** | Central method | Synthetic / fixture paths for residual and key columns; regression on column names used in figures |
| **DE backends used in paper** | Backend choice appears in Methods | Cover every backend the manuscript reports (e.g. scanpy, PyDESeq2, …) with mocks or small fixtures; CI `_de` gate stays green **with pydeseq2 installed** |
| **Filter presets / cutoffs** | Numbers appear in text | Preset names + semantics stable; document that *numeric* defaults may change pre-1.0 but paper must report installed version + values used |
| **Pipeline / `PipelineResult`** | One-liner workflow | Keys used in tutorials/paper (`candidates`, `enrichment`, `meta`, …) stable; pickle/deepcopy if you serialize results for figures |
| **Enrichment paths in paper** | GO/KEGG/GSEA if shown | ORA (and GSEA if claimed) runnable under CI `slow` or documented offline protocol; universe handling consistent with docs |

### T2 — Real-data smoke for the manuscript (must)

- Run the **exact analysis recipe** planned for figures on 1–2 real h5ad sets
  (user-owned or public with clear license).
- Record: `scatrans.__version__`, extras/backends, command or notebook path,
  pass/fail, unexpected warnings, seed if any.
- Prefer freezing that recipe into a notebook or script checked into docs/
  examples (even if data stay private).

### T3 — Support tests (should, not blockers for starting to write)

| Item | Notes |
|------|--------|
| Raise coverage on science-adjacent modules used in paper (`_velocity`, `pp_bias`, gene features) | Target only modules you **cite or rely on** in Methods |
| Default suite ≥ ~70–75% overall | Nice for confidence; not a paper gate by itself |
| Plot tests | Only if figures depend on `pl` helpers beyond “matplotlib on a DataFrame” |
| GSEA / with-slow honesty | Keep CI with-slow job; dual Codecov badges optional |

### T4 — Deliberately later (after paper, or only if paper needs it)

- KEGG out of default wheel / slim GO (compliance; do before wide commercial use)
- Optional `plot` extra, further `pl` / `tl.active` splits
- `py.typed`, CONTRIBUTING polish, dual Codecov badges
- Expanding public API

### Phase T exit checklist

- [ ] CI default suite green; `_de` gate green with pydeseq2
- [ ] Every analysis step **named in the paper** has a test or a recorded real-data script
- [ ] Real-data smoke log for manuscript figures
- [ ] Decision: **pin version** for Methods (e.g. `scatrans==0.10.x`) and stop feature work

---

## Phase M — Manuscript (while still 0.10.x)

Goal: a complete paper whose computational results are **reproducible from a
pinned Beta**, not from a moving target.

### Must align package ↔ paper

| Doc / artifact | Role in paper |
|----------------|---------------|
| `docs/method.md` | Method description consistency |
| `docs/statistical_guidance.md` | What scores may / may not claim |
| `docs/domain_assumptions.md` | Explicit assumptions (reviewer defense) |
| `docs/api_stability.md` | Import style recommended in code snippets |
| Tutorials / figure notebooks | Supplement or Zenodo bundle |
| `CITATION.cff` | Author: Zhao Li (李钊); ORCID when available |

### Writing discipline

1. **Pin** the package version in Methods (`scatrans` 0.10.x + key extras).
2. Report **cutoffs and presets**, not only “default”.
3. Prefer public API: `import scatrans as scat` / `scat.pl` / `scat.qc`.
4. No public API renames during revision unless a reviewer forces it—and then
   prefer aliases until 1.0.

### Phase M exit checklist

- [ ] Full draft (abstract → discussion) internally complete
- [ ] All figures/tables reproducible from pinned 0.10.x + recorded recipe
- [ ] Limitations match `statistical_guidance` / domain assumptions
- [ ] Software availability statement ready (GitHub, PyPI, license Apache-2.0;
      data licenses for GO/KEGG noted)

---

## Phase R — 1.0 + preprint (only after Phase M)

Order after the manuscript is **done** (ready to post):

1. **Freeze API** — `scatrans.__all__` + `pl.__all__` + `qc.__all__` =
   `api_stability.md`; no surprise renames.
2. **Release 1.0.0**
   - Bump `src/scatrans/_version.py` → `1.0.0`
   - Classifier → `Development Status :: 5 - Production/Stable`
   - CHANGELOG “stable surface” summary
   - SemVer note in docs (breaking = major after 1.0)
   - Optional: TestPyPI dry-run, then PyPI
3. **Preprint** — post with software version **1.0.0** (or “results produced
   with 0.10.x; software released as 1.0.0 with the same public API” if you
   must distinguish analysis pin vs marketing version—prefer **one clear
   story**: analysis redone or confirmed on 1.0 if time allows).
4. **Cite** — update `CITATION.cff` / README with preprint DOI when available.

### 1.0 packaging hygiene (same window, not science scope)

- [ ] KEGG bundling / commercial notice resolved or documented as opt-in
- [ ] Mouse gene-feature default unambiguous in docs
- [ ] ORCID on `CITATION.cff` if available
- [ ] Small public demo h5ad (Zenodo) if not already released with the paper

---

## After 1.0 (backlog only)

Science-led, not platform-led:

- Stronger design diagnostics / backend recommendation narrative
- Multi-sample / multi-condition reporting patterns
- Systematic comparison notebooks (nascent residual vs DE-only / velocity context)
- Species / custom gene-set ergonomics as users demand

Still non-goals unless strategy changes: agent stack, method zoo, GPU-first
rewrite.

---

## Engineering already done (0.9.x / 0.10.x track)

Reference only—do not re-open unless paper-blocking.

- [x] `tl` / `enrich` package split; public `__all__` / `__dir__`
- [x] `PipelineResult` (read-only dict subclass, copy/pickle/`to_dict`)
- [x] Matrix helpers; mypy clean + CI hard-fail (3.11 core)
- [x] Coverage split; `_de` ≥70% gate; optional stacks on coverage job
- [x] `tests/test_de_core_coverage.py`; `docs/api_stability.md`
- [x] Authors in metadata / README: Zhao Li (李钊), leelieber@gmail.com
- [x] Yank PyPI `*.dev*` (release side)

---

## How to re-verify quickly

```bash
# Types
mypy src/scatrans --ignore-missing-imports

# Default suite (daily CI)
pytest -m "not plot and not slow" -q

# _de coverage gate (install pydeseq2 for a fair measure)
pip install "pydeseq2>=0.4.0"
pytest -m "not plot and not slow" --cov=scatrans --cov-report=xml:coverage-default.xml -q

# With slow (GSEA etc.)
pytest -m "not plot" -q
```

```bash
python scripts/make_release_zips.py
```

---

## Working rules

1. **Paper first:** if a change does not help testing confidence or the
   manuscript, it waits until after 1.0.
2. **No API growth** in Phases T–M unless the paper needs a name that already
   exists conceptually—prefer private helpers over new `__all__` entries.
3. **Real h5ad smoke** is manuscript acceptance, not an excuse for more
   refactor.
4. If CI `_de` gate fails, check **pydeseq2 is installed** on the coverage job
   first (`.github/workflows/ci.yml`).

Last updated: 2026-07-12 — strategy: test → paper → 1.0 + preprint.
