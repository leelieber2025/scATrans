# API Reference

```python
import scatrans as scat
```

Public surface: top-level names in `scatrans.__all__`, plus `scat.pl` and
`scat.qc` (scanpy-style submodules). Other internal modules are not part of the
stable API. Narrative usage: {doc}`../user_guide/index`. Output semantics:
{doc}`../statistical_guidance`. Import stability: {doc}`../api_stability`.

## Core: scoring, filtering, DE

`partition_de_by_mechanism` is the recommended primary entry point (returns
{class}`~scatrans.PartitionResult`). `active_score` is the lower-level residual
and scoring engine; `differential_expression` is the DE path without nascent
layers. Scope and deprecated composite ranking: {doc}`../faq`.

**Optional detection columns (`partition_de_by_mechanism` only):**

| Parameter | Default | Effect |
|-----------|---------|--------|
| `add_nascent_score` | `False` | When `True`, append **detection** columns from {func}`~scatrans.nascent_activity_score` (`nascent_poisson_z`, `dlog_*`, `de_reproducible` / `de_repro_frac`) to the gene table. **Decoupled from mechanism:** the transcription-vs-stabilization call always uses the induction-normalized residual. Fail-soft (error recorded in `meta["nascent_score"]`). |

In the table below "both" means `active_score` and `differential_expression` (the
low-level scorers). The convenience entry points `partition_de_by_mechanism` /
`active_score_simple` / `run_default_pipeline` instead default
`target_group`/`reference_group` to `"Disease"`/`"Control"`.

| Parameter | Applies to | Default | Options |
|-----------|-----------|---------|---------|
| `groupby` | both | `"condition"` | any `obs` column holding group labels |
| `target_group` / `reference_group` | both | `None` (required for the low-level scorers) | must be set explicitly to two values in `adata.obs[groupby]`; tutorials often use e.g. `"GA"` / `"Ctrl"` |
| `use_pseudobulk` + `sample_col` | both | `False` / `None` | aggregate to per-replicate pseudobulk before DE (needs `sample_col`) |
| `pseudobulk_de_backend` | both | `"pydeseq2"` | `"pydeseq2"` (count-based DESeq2) or `"scanpy"` (rank_genes_groups on aggregated profiles) |
| `de_method` | both | `"t-test_overestim_var"` | any `scanpy.tl.rank_genes_groups` method, e.g. `"wilcoxon"` |
| `use_mixed_model` + `sample_col` | both | `False` | cell-level LMM with sample random intercept; needs ≥4 samples/group; `logFC` = sample-mean-of-means (not LMM coef); **incompatible with** `use_memento_de` |
| `use_memento_de` | both | `False` | method-of-moments cell-level DE (raw integer counts required); **incompatible with** `use_mixed_model` |
| `use_permutation` + `n_perm` + `perm_de_backend` | `active_score` | `False` / `100` / `"same"` | permutation FDR on unspliced excess; `perm_de_backend="fast"` trades accuracy for speed |
| `gamma_method` | `active_score` | `"heuristic_shrink"` | `"heuristic_shrink"`, `"robust_median"`, `"empirical_bayes"` (hierarchical, recommended for small reference groups), `"raw"` |
| `bias_correction` | `active_score` | `"huber_length_intron"` | `"huber_length_intron"` or `"none"` |
| `mode` | `active_score` | `"heuristic"` | `"heuristic"` or `"advanced"` (scVelo moments smoothing) |

Always call `recommend_workflow(...)` (or let `active_score_simple` /
`run_default_pipeline` call it for you) before picking these by hand — it
inspects cell/sample counts and suggests a preset.

```{eval-rst}
.. currentmodule:: scatrans

.. autosummary::
   :toctree: generated/
   :nosignatures:

   partition_de_by_mechanism
   PartitionResult
   active_score
   active_score_simple
   adaptive_active_score
   add_adaptive_score
   adaptive_weight
   labeling_anchor
   add_abundance_normalized_residual
   annotate_mechanism_class
   program_mechanism
   nascent_activity_score
   threshold_sensitivity
   differential_expression
   differential_expression_simple
   diagnose_design
   recommend_workflow
   run_default_pipeline
   filter_active_genes
   store_raw_counts
   ensure_raw_counts
   restore_raw_counts
```

**`PartitionResult` fields:** `adata`, `regime`, `gene_table`, `selected`,
`programs`, `enrichment`, `meta`. Useful `meta` keys include `de_source`,
`select` (cutoffs + `n_selected`), `regime`, `mechanism`, `programs`, and
`nascent_score` (`enabled` / `status` when `add_nascent_score` was used).
Default `run_go_enrichment=False` on partition (unlike
`run_default_pipeline`, which defaults to `True`).

`filter_active_genes(results_df, preset=..., select_by=..., logfc_direction=...,
return_mask=...)` accepts `preset="heuristic"` (standard cutoffs), `"pseudobulk"`
(looser, post-aggregation), `"significant"` (replays the built-in strict mask;
requires `use_permutation=True` upstream), or `"permissive"`; or pass explicit
`*_cutoff` kwargs instead of a preset. The **default is `preset=None`**
(permissive — only explicitly-passed cutoffs apply); `select_by="de"` additionally
applies DE defaults (`padj<0.05`, `logFC>1`) when no cutoffs are given.

**`select_by`** (default `"composite"`):

| Value | Membership decided by | Proxy columns |
|-------|----------------------|---------------|
| `"composite"` | DE gates **and** nascent/composite gates (prior behavior) | Participate in selection |
| `"de"` | **DE only** (`p_adj` / `logFC`; defaults `padj < 0.05` and `|log2FC| > 1` when no cutoffs given) | Annotation only — never gate membership. Sorted by `p_adj` then `logFC`. Incompatible with `preset="significant"` |

### Pipeline add-ons (`run_default_pipeline`)

Optional kwargs (all default off / prior behavior):

| Kwarg | Effect | `meta` key |
|-------|--------|------------|
| `bias_method="abundance"` / `"abundance_length"` | Adds `unspliced_excess_residual_abnorm` | `meta["bias"]` |
| `adaptive_weighting=True` | Adds `adaptive_score` / `adaptive_score_pct`; `adaptive_anchor=` selects reliability anchor (`"de"` or `labeling_anchor()` / callable) | `meta["adaptive"]` |
| `select_by="de"` | Candidates from DE gates only (proxy annotates) | `meta["select_by"]` |
| `annotate_mechanism=True` | Adds `transcription_support` / `mechanism_class` / `mechanism_confidence` (confidence scaled by regime reliability) | `meta["mechanism"]` |

`run_default_pipeline` **always** records `meta["regime"]` from
`scat.qc.regime_diagnosis` (cheap, fail-soft pre-flight from the global
unspliced fraction). Add-ons fail soft when columns are missing; invalid
`bias_method` raises.

### Post-hoc ranking & mechanism helpers (additive; core `active_score` unchanged)

| Function | Adds / returns | When to use |
|----------|----------------|-------------|
| `add_adaptive_score` / `adaptive_active_score` | `adaptive_score`, `adaptive_score_pct` | Reliability-weighted nascent leg; `anchor=` (`"de"` or `labeling_anchor(...)`) chooses the induced-gene set used for reliability AUC |
| `labeling_anchor` | callable for `anchor=` | Metabolic-labeling truth (e.g. `new_log2fc`) instead of DE-induced genes |
| `add_abundance_normalized_residual` | `unspliced_excess_residual_abnorm` | Demote abundance / nuclear-retention artifacts (e.g. *MALAT1*) |
| `nascent_activity_score` | `nascent_poisson_z`, `dlog_unspliced` / `dlog_spliced`, `de_reproducible` / `de_repro_frac` | Pseudobulk variance-stabilized nascent **detection** score + spliced-side DE-reproducibility flag; opt-in on partition via `add_nascent_score=True` (**not** the mechanism residual) |
| `annotate_mechanism_class` | `transcription_support`, `mechanism_class`, `mechanism_confidence` | Low-confidence per-gene transcription vs stabilization label (**annotation only**) |
| `program_mechanism` | program-level DataFrame | Threshold-free gene-set pooling of support (stronger than per-gene) |
| `threshold_sensitivity` | padj×logFC grid table | Report DE-list robustness instead of defending one cutoff |

See {doc}`../user_guide/advanced` and {doc}`../statistical_guidance`.

## Gene features

Used for optional bias correction inside `active_score` (length + intron
count regressed out of the raw unspliced-excess delta).

`WORKFLOW_PRESETS` is a public constant dict with keys `"explore"`, `"report"`,
`"pseudobulk_report"`, and `"nascent_focus"`. `recommend_workflow()` selects one
and returns its `active_score_kwargs` as `suggested_kwargs`; you apply them via
`scat.active_score(..., **rec["suggested_kwargs"])` (they are **not** a parameter of
`active_score_simple` / `run_default_pipeline`, which instead take `filter_preset`).
It is not autosummary'd here because it is data, not a callable.

The `generate-gene-features` console script maps to
`scatrans.generate_gene_features:main` (not re-exported on the top-level
`scatrans` package). Programmatic use is via `generate_gene_features_from_gtf()`
below, or `from scatrans.generate_gene_features import main` for the CLI
entrypoint function.

| Parameter | Function | Notes |
|-----------|----------|-------|
| `organism` | `add_gene_features` | `"mouse"` (default) or `"human"`, uses the bundled table |
| `gene_features_path` | `add_gene_features` | supply your own parquet instead of the bundled table |
| `gtf_path` / `organism` | `generate_gene_features_from_gtf` | build a custom table from a 10x/GENCODE GTF (`pip install "scatrans[gene_features]"`) |

```{eval-rst}
.. autosummary::
   :toctree: generated/
   :nosignatures:

   add_gene_features
   generate_gene_features_from_gtf
   list_available_gene_features
```

## Functional enrichment

All ORA-style functions (`run_enrichment`, `run_kegg`, `run_go`) share the
same universe-handling and gene-set-resolution machinery; `run_gsea` takes a
ranked list instead of a candidate gene list.

| Parameter | Applies to | Default | Options |
|-----------|-----------|---------|---------|
| `gene_sets` | `run_enrichment`, `run_gsea` | — | base name (e.g. `"GO_Biological_Process"`, auto-resolved per `organism`) or a full versioned Enrichr name (e.g. `"...2021"`) for a historical library |
| `organism` | all ORA/GSEA functions | — | `"mouse"` or `"human"` |
| `adata` | all ORA/GSEA functions | `None` | pass the object you called `store_raw_counts` on to auto-supply the measured-gene background/universe |
| `gene_set_source` | all ORA/GSEA functions | auto-detected | `"scatrans"` (bundled) or `"enrichr"` (gseapy) to force a source |
| `ontology` | `run_go` | `"BP"` | `"BP"`, `"CC"`, `"MF"`, or `"ALL"` (with `adjust_across_all=True` for unified correction) |
| `method` | `simplify_enrichment` | `"jaccard"` | `"jaccard"` (fast, overlap-based) or `"pathway_denester"` (nested-pathway test) |

```{eval-rst}
.. autosummary::
   :toctree: generated/
   :nosignatures:

   run_enrichment
   run_kegg
   run_go
   run_gsea
   simplify_enrichment
   save_enrichment_report
   expand_enrichment_genes
   list_bundled_gene_sets
   compare_enrichment
   extract_gene_lists
   concat_compare_results
```

## Quality control (`scat.qc`)

| Function | Returns | When to use |
|----------|---------|-------------|
| `qc.unspliced_global` | float in [0, 1] | Raw global unspliced fraction; logs a warning if high |
| `qc.regime_diagnosis` | dict: `unspliced_fraction`, `reliability` [0, 1], `regime` (`ok` / `low_unspliced` / `high_unspliced`), `basis`, `message` | Pre-flight proxy **data-quality** reliability (U-shaped map of fraction → reliability). Pass `reliability` into `annotate_mechanism_class`; `partition_de_by_mechanism` always runs it; `run_default_pipeline` stores it in `meta["regime"]` and uses it when `annotate_mechanism=True` |

**Scope:** `regime_diagnosis` is the **data-quality / gamma** half of a regime
check (too little nascent signal, or too much ≈ nuclear/gDNA → gamma mis-fit).
It does **not** yet distinguish dynamic vs steady-state transcription (that
needs a velocity-magnitude signal such as `velocity_length`, pending
validation). High `reliability` means the proxy is not clearly corrupted, not
that it outperforms DE.

```{eval-rst}
.. autosummary::
   :toctree: generated/
   :nosignatures:

   qc.unspliced_global
   qc.regime_diagnosis
```

## Plotting (`scat.pl`)

Every `scat.pl.*` function accepts `ax=`/`axes=` (multi-panel embedding),
`save_path=` (vector/300 dpi export), `show=`, and `use_style=`; most return
`(fig, ax)`.

| Parameter | Applies to | Notes |
|-----------|-----------|-------|
| `style` | `volcano_plot` | `"auto"` (legacy `active_score` colormap), `"ggvolcano"` (3-color classic), `"gradual"` (FDR gradient) |
| `context` | major plotters | `"notebook"` defaults vs `"paper"` (larger figsize/fonts, dpi=300); aliases include `"print"` / `"publication"` |
| `x` / `size_by` / `color_by` | `enrich_dotplot`, `enrich_barplot` | `x="GeneRatio"` (ORA) or `"NES"` (GSEA, auto-detected); `color_by` defaults to adjusted p-value |
| `show_terms` | `enrich_dotplot`, `compare_dotplot` | `int` (top N), `"auto"` (significance + count heuristic), or an explicit term list |
| `top_n` / `label_genes` / `label_repel` | `comet_plot`, `volcano_plot` | control auto-labeling; `label_genes=[...]` adds manual labels; `label_repel=False` skips adjustText |
| `s` / `point_scale` / `min_size` / `max_size` | comet / volcano / volcano_3d | `s=` forces fixed point size; otherwise score-based sizing with hard bounds |
| `use_style` / `set_style()` | all | opt-in publication style (off by default so notebooks are not surprised by global `rcParams` changes) |

`compare_dotplot` is the clusterProfiler ``compareCluster``-style multi-group
grid (groups on the x-axis, terms on the y-axis). Use it on the long table from
`compare_enrichment` / `concat_compare_results`. Prefer `enrich_dotplot` for a
single contrast or when faceting with `facet_by_cluster=True`.

`gene_upsetplot` is the **gene-level** UpSet (companion to the term-level
`enrich_upsetplot`): it shows how genes overlap across several DE results or
gene lists. Feed it either a `{name: de_df}` mapping (filtered internally) or a
pre-built membership matrix from `build_gene_membership`; in the default
`direction="separate"` mode each DE result contributes a `name::up` and
`name::down` set, so common-up and common-down genes appear as their own
intersection columns. `common_genes(membership, direction="up"|"down")` pulls
those intersection genes back out as a list ready for `run_enrichment`. Colors
are fully customizable (`set_color`, `intersection_color`, `dot_color`,
`inactive_color`, `line_color`; the intersection/dot colors also accept a
per-column list to highlight specific intersections). The same color parameters
were added to `enrich_upsetplot`; `bias_diagnostic_plot`
(`raw_color`/`corrected_color`/`trend_color`) and `gamma_shrinkage_plot`
(`cmap`/`color`) are now recolorable too.

```{eval-rst}
.. autosummary::
   :toctree: generated/
   :nosignatures:

   pl.comet_plot
   pl.volcano_plot
   pl.volcano_3d
   pl.bias_diagnostic_plot
   pl.enrich_dotplot
   pl.compare_dotplot
   pl.enrich_barplot
   pl.enrich_upsetplot
   pl.gene_upsetplot
   pl.build_gene_membership
   pl.common_genes
   pl.enrich_vennplot
   pl.gseaplot
   pl.active_score_rankplot
   pl.active_genes_heatmap
   pl.velocity_phase_portraits
   pl.gamma_shrinkage_plot
   pl.set_style
   pl.set_nature_style
   pl.style_context
   pl.figure_export_context
   pl.save_all_figures
```
