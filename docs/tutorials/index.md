# Tutorials

Worked examples on a real dataset: endothelial cells (EC) from mouse spinal
cord, comparing uninjured controls (**UN**, 3 replicates) against spinal
cord injury (**SCI**, 3 replicates). See {doc}`../references` for the full
citation and GEO accession.

```{toctree}
:maxdepth: 1

t_ec_active_transcription
t_ec_standalone_de_enrichment
```

- **{doc}`t_ec_active_transcription`** — for data **with** spliced/unspliced
  (or mature/nascent) layers: the full `active_score` workflow (heuristic,
  pseudobulk, permutation, gamma robustness, bias correction, advanced mode)
  end to end on real SCI vs. UN data.
- **{doc}`t_ec_standalone_de_enrichment`** — for data **without**
  spliced/unspliced layers: `differential_expression` across backends, plus
  a full tour of enrichment methods (ORA, KEGG, GO all-ontology, GSEA,
  redundancy reduction) and the plotting gallery.

## Reproducing these notebooks locally

Both notebooks ship with their outputs already rendered, so you can read
them here without running anything. To re-run them yourself:

```bash
git clone https://github.com/leelieber2025/scATrans.git
cd scATrans
pip install -e ".[dev,advanced,pseudobulk,gene_features,memento]" gseapy
jupyter lab docs/tutorials/t_ec_active_transcription.ipynb
```

`EC.h5ad` is included at the repository root, which is where both notebooks
load it from (`sc.read_h5ad("../../EC.h5ad")`).
