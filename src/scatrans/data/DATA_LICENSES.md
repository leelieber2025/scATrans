# Bundled data — license and redistribution notices

The **scATrans source code** is licensed under [Apache License 2.0](../../../LICENSE).
Files in this directory are **not** all covered by Apache-2.0. Each category below
has its own terms. By using bundled enrichment gene sets you agree to comply with
the applicable third-party licenses.

---

## 1. Gene feature parquet files

| File | License | Notes |
|------|---------|-------|
| `mouse_2020A_gene_features.parquet` | Apache-2.0 (project) | Derived by scATrans maintainers from public GTF annotations |
| `Mus_musculus.GRCm39.115_gene_features.parquet` | Apache-2.0 (project) | Same |
| `human_GRCh38_2024A_gene_features.parquet` | Apache-2.0 (project) | Same |

You may redistribute these tables under the project Apache-2.0 license.

---

## 2. GO Biological Process gene sets (`Hs_GO_*`, `Mm_GO_*`)

| File | Provenance | License |
|------|------------|---------|
| `Hs_GO_Biological_Process_2026.txt` | clusterProfiler GO map via `org.Hs.eg.db` + `GO.db` (Bioconductor) | Bioconductor / GO terms — see below |
| `Mm_GO_Biological_Process_2026.txt` | clusterProfiler GO map via `org.Mm.eg.db` + `GO.db` (Bioconductor) | Same |

- **Gene Ontology**: [GO citation / reuse policy](https://geneontology.org/docs/go-citation-policy/)
  ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) for ontology content;
  annotation sets may carry additional terms).
- **Bioconductor annotation packages** (`org.*.eg.db`, `GO.db`): subject to Bioconductor and upstream data-provider terms.

These files are **derived mappings** (gene symbol ↔ GO BP term) extracted for offline
enrichment. They are provided for academic research convenience. If you redistribute
a modified PyPI wheel or repackage these files, retain this notice and cite GO/Bioconductor
in your Methods section.

---

## 3. KEGG pathway gene sets (`Hs_KEGG_2026.txt`, `Mm_KEGG_2026.txt`)

| File | Provenance | License |
|------|------------|---------|
| `Hs_KEGG_2026.txt` | clusterProfiler KEGG cache (organism `hsa`) | **KEGG — not Apache-2.0** |
| `Mm_KEGG_2026.txt` | clusterProfiler KEGG cache (organism `mmu`) | **KEGG — not Apache-2.0** |

**Important:** KEGG pathway–gene mappings are subject to the
[KEGG FTP / academic use policy](https://www.kegg.jp/kegg/legal.html).

- **Academic use** (universities, non-profit research): generally permitted with
  attribution; confirm current KEGG terms before publication or redistribution.
- **Commercial use** (for-profit companies, commercial services/products): requires a
  **separate license agreement with KEGG** (Pathway Solutions / Kanehisa Laboratories).
  **Do not** assume the Apache-2.0 license on scATrans extends to KEGG data inside
  the wheel.

Bundling these files in the PyPI package is intended as a convenience for academic
offline ORA. If your use case is commercial, either:

1. Obtain a KEGG commercial license and verify redistribution is allowed, or
2. Do **not** use the bundled KEGG files — call `run_kegg(..., kegg_library="KEGG_2021")`
   (or another version) so gseapy/Enrichr fetches data at runtime under your own
   account/terms, or supply your own licensed gene-set file via `gene_sets=`.

scATrans maintainers make **no warranty** that redistribution of KEGG-derived mappings
in a Python package satisfies KEGG policy for your jurisdiction or use case.

---

## 4. Recommended citations / Methods text

- **GO**: Gene Ontology Consortium; cite the GO release used (see `gene_set_info` attrs).
- **KEGG**: Kanehisa M, Goto S; *KEGG: Kyoto Encyclopedia of Genes and Genomes* —
  cite [kegg.jp](https://www.kegg.jp/) and state whether use is academic or licensed commercial.
- **scATrans software**: see `CITATION.cff` in the repository root.

---

## 5. Maintainer note (future packaging)

A future release may move KEGG libraries to **optional runtime download** (not shipped
in the default wheel) to reduce redistribution risk. GO BP sets may remain bundled with
expanded provenance. Watch the CHANGELOG before upgrading.