# User Guide

The recommended entry point is {func}`~scatrans.partition_de_by_mechanism`: a
standard DE test **selects** the changed genes, then scATrans **partitions** them
into transcription- vs stabilization-driven — with a per-gene soft annotation and a
decisive program-level call. See {doc}`workflow` for the full walkthrough.

Under the hood it composes lower-level building blocks documented here:
`active_score` (differential expression, reference-gamma unspliced-excess, optional
bias correction, composite scoring + diagnostics), `filter_active_genes` for custom
gene lists (including `select_by="de"` for **DE selects, proxy annotates**), then
functional enrichment and plotting. The complete results table is `all_results`;
diagnostics live under `adata_res.uns["scatrans"]["diagnostics"]`.

:::{important}
Gene-list membership comes from **DE**; the spliced/unspliced signal only
**annotates** mechanism (a low-confidence per-gene hint — conclude at the program
level, and heed the reliability pre-flight). The composite
`run_default_pipeline(select_by="composite")` ranking is **deprecated**. DE,
enrichment, and plotting paths that do not depend on velocity layers are stable —
see {doc}`standalone_de` and the note on the {doc}`../index`.
:::

```{toctree}
:maxdepth: 2

workflow
enrichment
plotting
advanced
gene_features
standalone_de
```
