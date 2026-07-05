# References & Data Sources

## Tutorial dataset

The worked-example tutorials use `EC.h5ad`, a subset of endothelial cells
(EC) extracted from a public single-nucleus RNA-seq dataset of adult mouse
spinal cord, comparing **uninjured controls (UN, n=3 replicates)** against
**spinal cord injury (SCI, n=3 replicates)**:

> Squair, J.W., Gautier, M., Kathe, C., Anderson, M.A., James, N.D., Hutson,
> T.H., Hudelle, R., Qaiser, T., Matson, K.J.E., Barraud, Q., Levine, A.J.,
> La Manno, G., Skinnider, M.A., Courtine, G. (2021).
> **Confronting false discoveries in single-cell differential expression.**
> *Nature Communications* 12, 5692.
> DOI: [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)

- Raw data: GEO accession
  [GSE165003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165003)
  (6 samples: 3× UN, 3× SCI).

This paper is a particularly fitting choice for scATrans tutorials: it is
itself about the statistical pitfalls of treating cells as independent
replicates in single-cell differential expression — the same
pseudoreplication concern that motivates scATrans's pseudobulk /
mixed-model / permutation options and the {doc}`statistical_guidance` page.

## Methods and libraries scATrans builds on

| Component | Reference |
|-----------|-----------|
| scanpy (preprocessing, `rank_genes_groups` DE backends) | Wolf, F.A., Angerer, P., Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19, 15. DOI: [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0). [github.com/scverse/scanpy](https://github.com/scverse/scanpy) |
| PyDESeq2 (`pseudobulk_de_backend="pydeseq2"`) | Muzellec, B., Teleńczuk, M., Cabeli, V., Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39(9), btad547. DOI: [10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547). [github.com/scverse/PyDESeq2](https://github.com/scverse/PyDESeq2) |
| scVelo (`mode="advanced"` moments smoothing) | Bergen, V., Lange, M., Peidli, S., Wolf, F.A., Theis, F.J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology* 38, 1408–1414. DOI: [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3). [github.com/theislab/scvelo](https://github.com/theislab/scvelo) |
| Memento (`use_memento_de=True`) | Kim, M.C., Gate, R., Lee, D.S., Tolopko, A., Lu, A., Gordon, E., Shifrut, E., Garcia-Nieto, P.E., Marson, A., Ntranos, V., Ye, C.J. (2024). Method of moments framework for differential expression analysis of single-cell RNA sequencing data. *Cell* 187(22), 6393–6410.e16. DOI: [10.1016/j.cell.2024.09.044](https://doi.org/10.1016/j.cell.2024.09.044). [github.com/yelabucsf/scrna-parameter-estimation](https://github.com/yelabucsf/scrna-parameter-estimation) |
| GSEApy (Enrichr access, `run_gsea` prerank engine) | [github.com/zqfang/GSEApy](https://github.com/zqfang/GSEApy), docs at [gseapy.readthedocs.io](https://gseapy.readthedocs.io/en/latest/) |
| ggVolcano-style volcano plots (`volcano_plot(style="ggvolcano"/"gradual")`) | [github.com/BioSenior/ggVolcano](https://github.com/BioSenior/ggVolcano) |
| PathwayDenester (`simplify_enrichment(method="pathway_denester")`) | [github.com/Helmy-Lab/PathwayDenester](https://github.com/Helmy-Lab/PathwayDenester) |

## Bundled gene set data provenance

GO and KEGG gene sets bundled with scATrans have their own provenance and
licensing terms (GO is CC BY 4.0-derived; KEGG requires a commercial license
for non-academic redistribution). See
[`src/scatrans/data/DATA_LICENSES.md`](https://github.com/leelieber2025/scATrans/blob/main/src/scatrans/data/DATA_LICENSES.md)
and {doc}`license` for details.

## Citing scATrans

If you use scATrans in your research, please cite it using the metadata in
[`CITATION.cff`](https://github.com/leelieber2025/scATrans/blob/main/CITATION.cff)
(also exposed as GitHub's "Cite this repository" button).
