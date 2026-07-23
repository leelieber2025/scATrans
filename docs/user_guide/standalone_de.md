# Standalone Differential Expression

Use this path when nascent (spliced/unspliced) layers are unavailable or when
only conventional DE is required.

```python
import scatrans as scat

# Early (right after load + basic QC, before HVG/normalize/log).
# sidecar=True (default) snapshots the full-gene counts into .uns so they
# survive later HVG/cell subsetting.
scat.store_raw_counts(adata, layer="counts")

# Standard count AnnData (no spliced/unspliced layers required)
adata, de_results = scat.differential_expression(
    adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    # de_method="t-test_overestim_var",   # or "wilcoxon", etc. (default)
    # use_memento_de=True,                # optional: use the integrated Memento (Cell 2024) backend
    # memento_capture_rate=0.07,
)

# Then use the same downstream tools as with active_score results
candidates = scat.filter_active_genes(de_results, padj_cutoff=0.05, logfc_cutoff=0.3)  # upregulated (default)
# downregulated: logfc_direction="down"
# both: logfc_direction="both"
# DE-only defaults (padj<0.05 & |log2FC|>1 when cutoffs omitted):
# candidates = scat.filter_active_genes(de_results, select_by="de")

# After scat.store_raw_counts(adata) early in the workflow,
# just pass adata= here. It auto-supplies the full measured gene list as background/universe.
enrich = scat.run_enrichment(
    candidates.index.tolist(),
    gene_sets="GO_Biological_Process",  # auto → correct Hs/Mm 2026 bundled
    adata=adata,
)
scat.pl.volcano_plot(de_results)
scat.pl.enrich_dotplot(enrich)
```

`differential_expression` supports the same flexible backends as
`active_score` (scanpy methods, PyDESeq2 pseudobulk, mixed models, and
optionally Memento as a method-of-moments estimator). **Do not enable
`use_mixed_model` and `use_memento_de` together** — they are mutually
exclusive. With MixedLM, reported `logFC` is sample-mean-of-means log2FC
(see {doc}`advanced`). The returned table is directly compatible
with `filter_active_genes`, enrichment functions, and all `scat.pl.*`
plotting helpers.

Example script: `examples/memento_de_example.py`.

## Raw counts requirement

Count-based backends (Memento, PyDESeq2) expect raw integer counts. The following
pattern leaves unsuitable data:

```python
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

This leaves `adata.X` as log-transformed HVGs only, which is unsuitable.

**Early in the workflow:**

```python
import scatrans as scat

# Before HVG + normalize + log1p
scat.store_raw_counts(adata, layer="counts")   # saves raw counts to layers["counts"] + .uns snapshot

# Then normal Scanpy preprocessing
sc.pp.highly_variable_genes(adata, ...)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Now safe (count-based backends read the snapshot / counts layer automatically)
adata, de_res = scat.differential_expression(adata, use_memento_de=True, ...)

# To run DE / enrichment on the full pre-HVG gene set:
adata_full = scat.restore_raw_counts(adata, full_genes=True)
```

`store_raw_counts(adata, mode="auto")` is idempotent and will also try to recover
counts from `adata.raw` when `.X` is already normalized/log-transformed. The
functions emit clear warnings when they detect this situation.

**Note on `anndata.concat()` and `de_preprocess="auto"`:** `ad.concat()`
drops `.uns` by default, including the `uns['log1p']` marker that scATrans
uses to detect already-log-normalized data. This is common when combining
multiple samples for a case-vs-control comparison — each sample may be
correctly `normalize_total` + `log1p`'d before concatenation, but the marker
is gone afterward. `de_preprocess="auto"` still guards against double-log1p
in this case via heuristics on `.X`, but for certainty, either re-set the
marker after concatenating (`combined.uns["log1p"] = {"base": None}`) or
pass `de_preprocess="none"` explicitly when you know `.X` is already
log-normalized.
