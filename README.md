# scATrans

[![PyPI version](https://img.shields.io/pypi/v/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Python versions](https://img.shields.io/pypi/pyversions/scatrans.svg)](https://pypi.org/project/scatrans/)
[![Documentation Status](https://readthedocs.org/projects/scatrans/badge/?version=latest)](https://scatrans.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml/badge.svg)](https://github.com/leelieber2025/scATrans/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

scATrans is a Python toolkit for single-cell differential analysis. It is
primarily designed for datasets that contain spliced/unspliced (or
mature/nascent) RNA layers. In this setting it computes a composite active
transcription score that integrates differential expression with
reference-based excess unspliced RNA to rank genes.

It also supports conventional differential expression workflows (no
velocity data required) using scanpy, PyDESeq2 pseudobulk, linear mixed
models, or optional Memento. Functional enrichment (ORA, GSEA, GO, KEGG)
uses bundled gene sets with consistent universe handling, and a set of
visualization functions is provided.

**📚 Full documentation, tutorials, and the complete API reference are on
Read the Docs: https://scatrans.readthedocs.io**

## Installation

```bash
pip install scatrans

# Optional extras: advanced (scVelo) mode, gene features, pseudobulk DE (PyDESeq2), Memento, GSEA
pip install "scatrans[advanced,gene_features,pseudobulk,memento,gsea]"
```

See [Installation](https://scatrans.readthedocs.io/en/latest/installation.html)
for extras, source installs, and logging setup.

When developing from a git checkout, install the **editable** tree you are
editing (`pip install -e .` from that checkout). A stale editable install
pointing at another path can make `import scatrans` load an older copy.

## Quickstart

```python
import scatrans as scat

# One-liner pipeline: score → filter → GO enrichment
result = scat.run_default_pipeline(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    sample_col="sample",   # optional; auto-selects pseudobulk when >=3 replicates/group
    organism="mouse",
)
print(result["candidates"].head())
print(result["enrichment"].head())
```

See the [Quickstart](https://scatrans.readthedocs.io/en/latest/quickstart.html)
for a complete end-to-end walkthrough, the
[Tutorials](https://scatrans.readthedocs.io/en/latest/tutorials/index.html)
for fully worked, real-data notebooks (with and without RNA-velocity
layers), and the
[User Guide](https://scatrans.readthedocs.io/en/latest/user_guide/index.html)
for DE backends, enrichment, plotting, and advanced options.

## Before reporting results in a paper

`active_score` is a **composite heuristic rank**, not a p-value or FDR on
its own. See
[Statistical Guidance](https://scatrans.readthedocs.io/en/latest/statistical_guidance.html)
for what each output column means, safe vs. unsafe uses, and a reporting
checklist before you cite scATrans results in a manuscript or supplement.

## API stability

scATrans is **Beta (0.10.x)**. Prefer `import scatrans as scat` and the names
in `scatrans.__all__` / `scat.pl` / `scat.qc`. Leaf modules such as
`scatrans.tl.active` are **implementation detail** and may move. Full
contract: [API stability](https://scatrans.readthedocs.io/en/latest/api_stability.html)
(source: `docs/api_stability.md`).


## License

Software (Python source) is licensed under [Apache License 2.0](LICENSE).
Bundled gene-set data (GO, KEGG) carries its own licensing terms — see
[License](https://scatrans.readthedocs.io/en/latest/license.html) before
commercial use.
