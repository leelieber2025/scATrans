# License

**Software (Python source):** [Apache License 2.0](https://github.com/leelieber2025/scATrans/blob/main/LICENSE).

**Bundled data** in `src/scatrans/data/` is **not** entirely Apache-2.0:

| Data | License |
|------|---------|
| Gene feature `.parquet` files | Apache-2.0 (project) |
| GO Biological Process `*_2026.txt` | GO / Bioconductor-derived — see [DATA_LICENSES.md](https://github.com/leelieber2025/scATrans/blob/main/src/scatrans/data/DATA_LICENSES.md) |
| KEGG `Hs_KEGG_2026.txt` / `Mm_KEGG_2026.txt` | **[KEGG terms](https://www.kegg.jp/kegg/legal.html)** — academic use with attribution; commercial use requires a separate KEGG license. Not redistributable under Apache-2.0. |

If you ship a product or offer a commercial service, review
[DATA_LICENSES.md](https://github.com/leelieber2025/scATrans/blob/main/src/scatrans/data/DATA_LICENSES.md)
before using bundled KEGG mappings. You can avoid bundled KEGG files by
passing an Enrichr/gseapy library name (e.g. `kegg_library="KEGG_2021"`) or
your own `gene_sets` file.
