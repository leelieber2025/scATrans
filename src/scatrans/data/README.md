# Data Folder

This folder contains two kinds of data:

1. Gene feature files (for bias correction in active transcription scoring)
2. Bundled gene sets for functional enrichment (GO / KEGG)

---

## 1. Gene feature files (bias correction)

- `mouse_2020A_gene_features.parquet`
- `Mus_musculus.GRCm39.115_gene_features.parquet`
- `human_GRCh38_2024A_gene_features.parquet`

These provide `gene_length` + `intron_number`. `add_gene_features(..., organism="mouse")` picks a mouse table; `organism="human"` picks the human table. Used by `add_gene_features()` / `active_score()`.

You can generate similar files for other species or custom annotations using:

CLI (recommended):
```bash
pip install "scatrans[gene_features]"
generate-gene-features --gtf /path/to/genes.gtf \
                       --output human_GRCh38_gene_features.parquet \
                       --organism human
```

Programmatic:
```python
import scatrans as scat
df = scat.generate_gene_features_from_gtf(
    "/path/to/genes.gtf",
    output_name="my_features.parquet",
    organism="human"
)
```

Then attach to your data:
```python
adata = scat.add_gene_features(adata, gene_features_path="my_features.parquet")
```

The resulting table must have columns `gene_name`, `gene_length`, `intron_number` (and optionally `gene_type`).

---

## 2. Bundled enrichment gene sets (scATrans / clusterProfiler-derived)

We provide GO and KEGG gene sets extracted from clusterProfiler (using the same TERM2GENE mappings), named with the package name (`*_scATrans.gmt`).

These are stored as standard GMT files so they are compatible with the existing `run_enrichment` / `run_kegg` API.

`run_kegg()` (and `run_enrichment` with base GO names) now defaults to the organism-specific built-in libraries (Hs_*/Mm_*_2026.txt) when you only specify `organism`. The old _scATrans names are still supported for backward compatibility.

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

# Default — organism-specific built-in library (Hs_*/Mm_*_2026)
# You only need to give organism for KEGG (or a base GO name).
kegg = scat.run_kegg(my_genes, organism="mouse")

# GO base name — automatically resolved to the correct Hs/Mm 2026 built-in
# (only the BP + KEGG 2026 files are bundled; other GO branches fall back to gseapy/Enrichr if named explicitly)
res = scat.run_enrichment(
    gene_list=my_genes,
    gene_sets="GO_Biological_Process",  # or "GO_BP"
    # Best: pass adata= after store_raw_counts(...) for automatic correct universe
    adata=adata,
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

### License and redistribution (read before commercial use)

Bundled enrichment files have **separate** license terms from the scATrans **source code**
(Apache-2.0). See **[DATA_LICENSES.md](DATA_LICENSES.md)** in this folder.

Summary:

- **Gene feature `.parquet` files** — Apache-2.0 (project).
- **GO BP `*_GO_Biological_Process_2026.txt`** — GO / Bioconductor-derived; cite GO and
  Bioconductor annotation packages in Methods.
- **KEGG `Hs_KEGG_2026.txt` / `Mm_KEGG_2026.txt`** — **not Apache-2.0**. KEGG pathway
  mappings are subject to [KEGG legal terms](https://www.kegg.jp/kegg/legal.html).
  Academic use is typically allowed with attribution; **commercial use requires a KEGG
  license**. Do not redistribute these files outside academic fair-use without checking
  KEGG policy. For commercial pipelines, use runtime Enrichr/gseapy KEGG libraries or
  your own licensed gene-set file instead of the bundled `.txt` files.

### Provenance of bundled `*_2026.txt` libraries (for Methods / reproducibility)

The four default organism-specific libraries shipped with scATrans 0.10.x are:

| File | Species | Terms (approx.) | Source pipeline |
|------|---------|-----------------|-----------------|
| `Hs_GO_Biological_Process_2026.txt` | Human | 14,208 | clusterProfiler GO map (`org.Hs.eg.db` + `GO.db`) |
| `Mm_GO_Biological_Process_2026.txt` | Mouse | 14,956 | clusterProfiler GO map (`org.Mm.eg.db` + `GO.db`) |
| `Hs_KEGG_2026.txt` | Human | 222 | clusterProfiler KEGG cache (organism `hsa`) — **KEGG license, not Apache-2.0** |
| `Mm_KEGG_2026.txt` | Mouse | 218 | clusterProfiler KEGG cache (organism `mmu`) — **KEGG license, not Apache-2.0** |

- **Extracted:** June 2026 (`extracted_date: 2026-06` in enrichment attrs).
- **Gene ID type:** gene symbols (same convention as clusterProfiler ORA output).
- **Format:** GMT-like tab-separated `.txt` (term, optional empty description column, genes).
- **Runtime access:** `run_enrichment` / `run_go` / `run_kegg` write `gene_set_info.provenance` into result `.attrs` when a bundled library is loaded.

**Suggested Methods sentence:**

> Functional enrichment used scATrans bundled GO Biological Process and KEGG gene sets (clusterProfiler-derived mappings, June 2026 release; `Hs_*/Mm_*_2026` libraries), with Benjamini–Hochberg FDR across tested terms (`p_adjust_method='fdr_bh'`).

For Enrichr historical versions, cite the exact library name passed to `gene_sets=` (e.g. `GO_Biological_Process_2021`) and note `actual_source: gseapy` in `gene_set_info`.
