# References & Data Sources

This page lists **external methods and packages** scATrans depends on or
implements wrappers for, with **citation metadata** and **verified project
links**. Prefer citing the original method papers (and scATrans via
`CITATION.cff`) in your Methods section.

---

## Tutorial datasets

### EC (`EC.h5ad` — UN vs SCI)

The endothelial-cell tutorials use `EC.h5ad`, a subset of endothelial cells
(EC) extracted from a public single-nucleus RNA-seq dataset of adult mouse
spinal cord, comparing **uninjured controls (UN, n=3 replicates)** against
**spinal cord injury (SCI, n=3 replicates)**:

> Squair, J.W., Gautier, M., Kathe, C., Anderson, M.A., James, N.D., Hutson,
> T.H., Hudelle, R., Qaiser, T., Matson, K.J.E., Barraud, Q., Levine, A.J.,
> La Manno, G., Skinnider, M.A., Courtine, G. (2021).
> **Confronting false discoveries in single-cell differential expression.**
> *Nature Communications* 12, 5692.
> DOI: [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2)
> · PubMed: [34531461](https://pubmed.ncbi.nlm.nih.gov/34531461/)

- Raw data: GEO accession
  [GSE165003](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165003)
  (6 samples: 3× UN, 3× SCI).

This paper is a particularly fitting choice for scATrans tutorials: it is
itself about the statistical pitfalls of treating cells as independent
replicates in single-cell differential expression — the same
pseudoreplication concern that motivates scATrans's pseudobulk /
mixed-model / permutation options and the {doc}`statistical_guidance` page.

### GA (`GA_test.h5ad` — GA vs Ctrl)

The higher-power active-transcription tutorial uses `GA_test.h5ad`, a bundled
mouse scRNA-seq test object (GA vs. Ctrl, three individuals per group, with
`spliced` / `unspliced` layers). The underlying experiment is public
single-cell RNA-seq of cells from a mouse calvarial bone defect after TrkA
agonism with gambogic amide (**GA**) versus vehicle control:

> Li, Z., Xing, X., Du, B., Zhou, M., Chen, A.Z., Archer, M., Rao, C.,
> Zhu, M., Cherief, M., James, A.W. (2026).
> **Boosting Sensory Nerve-to-Bone Interactions Enhances Hedgehog Mediated
> Calvarial Bone Repair.**
> *Advanced Science* e75389.
> DOI: [10.1002/advs.75389](https://doi.org/10.1002/advs.75389)
> · Publisher:
> [Wiley Advanced Science](https://advanced.onlinelibrary.wiley.com/doi/abs/10.1002/advs.75389)
> · PubMed: [42003742](https://pubmed.ncbi.nlm.nih.gov/42003742/)

- Raw data: GEO accession
  [GSE266598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266598)
  (NCBI GEO).

---

## Differential expression & single-cell backends

| Component | Citation | Links |
|-----------|----------|-------|
| **scanpy** (preprocessing, `rank_genes_groups` DE) | Wolf, F.A., Angerer, P., Theis, F.J. (2018). SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 19, 15. DOI: [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0) | [GitHub](https://github.com/scverse/scanpy) |
| **DESeq2** (statistical model behind PyDESeq2) | Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15, 550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8) · PubMed: [25516281](https://pubmed.ncbi.nlm.nih.gov/25516281/) | — |
| **PyDESeq2** (`pseudobulk_de_backend="pydeseq2"`) | Muzellec, B., Teleńczuk, M., Cabeli, V., Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* 39(9), btad547. DOI: [10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547) · PubMed: [37669147](https://pubmed.ncbi.nlm.nih.gov/37669147/) | [GitHub](https://github.com/scverse/PyDESeq2) |
| **Memento** (`use_memento_de=True`) | Kim, M.C., et al. (2024). Method of moments framework for differential expression analysis of single-cell RNA sequencing data. *Cell* 187(22), 6393–6410.e16. DOI: [10.1016/j.cell.2024.09.044](https://doi.org/10.1016/j.cell.2024.09.044) | [GitHub](https://github.com/yelabucsf/scrna-parameter-estimation) |
| **statsmodels MixedLM** (`use_mixed_model=True`) | Seabold, S., Perktold, J. (2010). Statsmodels: Econometric and statistical modeling with Python. *Proc. 9th Python in Science Conf.* | [MixedLM docs](https://www.statsmodels.org/stable/mixed_linear.html) · [GitHub](https://github.com/statsmodels/statsmodels) |

---

## RNA velocity / moments (advanced track)

| Component | Citation | Links |
|-----------|----------|-------|
| **RNA velocity** (concept) | La Manno, G., et al. (2018). RNA velocity of single cells. *Nature* 560, 494–498. DOI: [10.1038/s41586-018-0414-6](https://doi.org/10.1038/s41586-018-0414-6) | — |
| **scVelo** (`mode="advanced"` moments smoothing) | Bergen, V., Lange, M., Peidli, S., Wolf, F.A., Theis, F.J. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. *Nature Biotechnology* 38, 1408–1414. DOI: [10.1038/s41587-020-0591-3](https://doi.org/10.1038/s41587-020-0591-3) | [GitHub](https://github.com/theislab/scvelo) |

scATrans’s default track is **reference-gamma unspliced excess** (group contrast),
not a full dynamical velocity fit. Do not cite scVelo as the engine of
`mode="heuristic"` results.

---

## Enrichment (ORA / GSEA / gene sets)

| Component | Citation | Links |
|-----------|----------|-------|
| **GSEA** (preranked enrichment concept) | Subramanian, A., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102(43), 15545–15550. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102) · PubMed: [16199517](https://pubmed.ncbi.nlm.nih.gov/16199517/) | — |
| **GSEApy** (`run_gsea`, Enrichr access) | Fang, Z., Liu, X., Peltz, G. (2023). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics* 39(1), btac757. DOI: [10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757) · PubMed: [36426870](https://pubmed.ncbi.nlm.nih.gov/36426870/) | [GitHub](https://github.com/zqfang/GSEApy) · [Docs](https://gseapy.readthedocs.io/en/latest/) |
| **Enrichr** (gene-set libraries via gseapy) | Kuleshov, M.V., et al. (2016). Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Research* 44(W1), W90–W97. DOI: [10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377) · PubMed: [27141961](https://pubmed.ncbi.nlm.nih.gov/27141961/) | [maayanlab.cloud/Enrichr](https://maayanlab.cloud/Enrichr/) |
| **clusterProfiler** (ORA / GO workflow inspiration) | Yu, G., Wang, L.G., Han, Y., He, Q.Y. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16(5), 284–287. DOI: [10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118) · PubMed: [22455463](https://pubmed.ncbi.nlm.nih.gov/22455463/) | — |
| **Gene Ontology** (bundled BP sets) | Gene Ontology Consortium. Cite the GO release used; see [GO citation policy](https://geneontology.org/docs/go-citation-policy/) (CC BY 4.0 for ontology content). | [geneontology.org](https://geneontology.org/) |
| **KEGG** (bundled pathway sets) | Kanehisa, M., Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research* 28(1), 27–30. DOI: [10.1093/nar/28.1.27](https://doi.org/10.1093/nar/28.1.27) · PubMed: [10592173](https://pubmed.ncbi.nlm.nih.gov/10592173/) | [KEGG legal terms](https://www.kegg.jp/kegg/legal.html) · [kegg.jp](https://www.kegg.jp/) |
| **PathwayDenester** (`simplify_enrichment(method="pathway_denester")`) | Helmy Lab PathwayDenester (repository). | [GitHub](https://github.com/Helmy-Lab/PathwayDenester) |

---

## Multiple testing, robust regression, permutation p-values

| Component | Citation | Links |
|-----------|----------|-------|
| **Benjamini–Hochberg FDR** (`p.adjust` / `multipletests` FDR_bh) | Benjamini, Y., Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J. R. Stat. Soc. B* 57(1), 289–300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x) | — |
| **Permutation *p* = (1+exceed)/(n+1)** | Phipson, B., Smyth, G.K. (2010). Permutation *P*-values should never be zero: calculating exact *P*-values when permutations are randomly drawn. *Stat. Appl. Genet. Mol. Biol.* 9, Article 39. DOI: [10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585) · PubMed: [21044043](https://pubmed.ncbi.nlm.nih.gov/21044043/) | — |
| **Huber regression** (`sklearn.linear_model.HuberRegressor` bias correction) | Huber, P.J. (1964). Robust estimation of a location parameter. *Ann. Math. Statist.* 35(1), 73–101. Implementation: scikit-learn. | [HuberRegressor docs](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.HuberRegressor.html) |

---

## Visualization utilities

| Component | Citation / source | Links |
|-----------|-------------------|-------|
| **ggVolcano-style** volcano (`style="ggvolcano"` / `"gradual"`) | BioSenior/ggVolcano (R package inspiration for layout). | [GitHub](https://github.com/BioSenior/ggVolcano) |

---

## Bundled gene set data provenance

GO and KEGG gene sets bundled with scATrans have their own provenance and
licensing terms (GO is CC BY 4.0-derived; KEGG requires a commercial license
for non-academic redistribution). See:

- Local: `src/scatrans/data/DATA_LICENSES.md`
- GitHub: [DATA_LICENSES.md](https://github.com/leelieber2025/scATrans/blob/main/src/scatrans/data/DATA_LICENSES.md)
- {doc}`license`

---

## Citing scATrans

If you use scATrans in your research, please cite it using the metadata in
[`CITATION.cff`](https://github.com/leelieber2025/scATrans/blob/main/CITATION.cff)
(also exposed as GitHub's "Cite this repository" button).

Also cite the **backends and databases you actually used** (e.g. PyDESeq2 +
DESeq2; GSEApy + GSEA + Enrichr; GO/KEGG) from the tables above.
