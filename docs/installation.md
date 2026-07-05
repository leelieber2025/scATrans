# Installation

```bash
# Basic installation
pip install scatrans

# With support for scVelo-based advanced mode and the gene feature generation CLI
pip install "scatrans[advanced,gene_features]" gseapy

# With support for pseudobulk differential expression using PyDESeq2
pip install "scatrans[pseudobulk]"

# Optional: Memento (Cell 2024) as an additional cell-level DE backend
pip install "scatrans[memento]"
```

The package ships precomputed gene feature tables (gene length + intron
number) for both mouse and human. These are used for optional bias
correction in `active_score`. You can also supply custom tables (e.g. from
your own GTF) — see {doc}`user_guide/gene_features`.

## Install from source

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev]"
```

## Logging

The package logs under the name `scatrans`. You can control verbosity with:

```python
import logging
logging.getLogger("scatrans").setLevel(logging.INFO)
```

## Quick data quality check

Before analysis, inspect the global unspliced fraction:

```python
import scatrans as scat
ufrac = scat.qc.unspliced_global(adata)   # logs INFO + WARNING if > 50%
```

`active_score` automatically runs this check and records the value in
diagnostics.
