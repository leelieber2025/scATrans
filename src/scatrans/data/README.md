# Data Folder

This folder contains two kinds of data:

1. Gene feature files (for bias correction in active transcription scoring)
2. Bundled gene sets for functional enrichment (GO / KEGG)

---

## 1. Gene feature files (bias correction)

- `mouse_2020A_gene_features.parquet`
- `Mus_musculus.GRCm39.115_gene_features.parquet`

These provide `gene_length` + `intron_number` for mouse. Used by `add_gene_features()` / `active_score()`.

You can generate similar files for other species using the CLI or `pp_bias.generate_gene_features_from_gtf`.

---

## 2. Bundled enrichment gene sets (scATrans / clusterProfiler-derived)

We provide GO and KEGG gene sets extracted from clusterProfiler (using the same TERM2GENE mappings), named with the package name (`*_scATrans.gmt`).

These are stored as standard GMT files so they are compatible with the existing `run_enrichment` / `run_kegg` API. `run_kegg()` now defaults to the bundled scATrans version.

The original Enrichr collections (via gseapy) have many historical versions (different years). You can select any specific version by passing the exact name together with `gene_set_source="enrichr"`.

### Why use these instead of (or in addition to) Enrichr/gseapy sets?

- **Consistency**: Same gene ↔ term associations as clusterProfiler → more comparable results.
- Works great together with the improved background/universe handling we added (conservative intersect by default, `force_universe`, detailed `universe_info` in attrs).
- You can still use the original Enrichr libraries (e.g. `"GO_Biological_Process_2023"`, `"KEGG_2026"`) — the loader tries bundled sets first, then falls back to gseapy.

### Usage

```python
import scatrans as scat

# List what is bundled
print(scat.list_bundled_gene_sets())

# Default — package's bundled scATrans version (only need organism for KEGG)
kegg = scat.run_kegg(my_genes, organism="mouse")

# For GO, base name is enough — gets bundled version automatically
res = scat.run_enrichment(
    gene_list=my_genes,
    gene_sets="GO_Biological_Process_2023",
    background=background_genes,
)

# To use original Enrichr or a specific historical version: just write the full name
classic = scat.run_enrichment(
    gene_list=my_genes,
    gene_sets="GO_Biological_Process_2021",  # or 2023, 2019, 2018...
)

classic_kegg = scat.run_kegg(
    my_genes, organism="mouse",
    kegg_library="KEGG_2021",   # or KEGG_2019, KEGG_2016, etc.
)
```

### How to add / name your files

Drop your GMT files directly into this `data/` directory (they will be packaged automatically).

Recommended naming convention (uses the package name for clarity):
- `GO_Biological_Process_scATrans.gmt`
- `GO_Cellular_Component_scATrans.gmt`
- `GO_Molecular_Function_scATrans.gmt`
- `KEGG_scATrans.gmt`

The loader accepts the name **with or without the `.gmt` extension**.

After adding files, users (and you) can discover them with:
```python
scat.list_bundled_gene_sets(verbose=True)
```

### File format

Standard GMT (tab-separated):
```
TERM_ID<TAB>description<TAB>GENE1<TAB>GENE2<TAB>...
```

This is the same format Enrichr/gseapy GMTs use, so any existing GMT tooling will work.
