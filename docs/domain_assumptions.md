# Domain Assumptions

Domain conventions enforced by the implementation. Use this page to check
expectations against behavior.

Related: {doc}`statistical_guidance` (reporting), {doc}`faq` (scope).

---

## 1. What “active” means

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **Active transcription ranking is upregulation-oriented** | Built-in `significant` requires positive `logFC` (and positive residual). Default `filter_active_genes` direction is `"up"`. | `logfc_direction="down"` / `"both"` on pure DE tables; or use `differential_expression` alone |
| **The unspliced-excess residual is independent of DE direction and DE significance** | Positive unspliced excess can occur even if DE is weak, untested, or `p_adj` filled to 1 (e.g. PyDESeq2 independent filtering). Ranking by the residual is not a DE-significant list. | `filter_active_genes(..., select_by="de")` or explicit `padj_cutoff` / `logfc_cutoff`, or `significant` conjunction |
| **`select_by="de"` decouples membership from the proxy** | DE gates define the list; nascent gates are skipped; sorting is by `p_adj` then `logFC`. | Default `select_by="composite"` keeps prior proxy-aware filtering |
| **`qc.regime_diagnosis` reliability is U-shaped in unspliced fraction** | Full reliability in a normal band (~10–45%); degrades at low (noise) and high (nuclear/gDNA, gamma mis-fit) extremes. Scales `mechanism_confidence` in `partition_de_by_mechanism` and when `annotate_mechanism=True`. | This is data-quality only — not dynamic vs steady-state |
| **Detection ≠ mechanism** | `nascent_poisson_z` / `add_nascent_score=True` is an absolute nascent-increase **detection** score (induction-coupled). Mechanism labels always use the induction-normalized residual. | Do not pass `nascent_poisson_z` as `residual_col` for `annotate_mechanism_class` if you want transcription-vs-stabilization |
| **Residual is one-sided (positive excess)** | Negative excess does not contribute to the residual soft-scale. Permutation FDR on residual is one-sided for positive excess. | Do not interpret low residual FDR as “repression” |
| **The residual is a within-run relative magnitude** | The residual soft-scale `λ ≈ median(positive x) / ln(2)` is estimated **from the genes in that run**. Residual magnitudes are comparable within one run, not as absolute cross-run units. | For cross-subset or cross-dataset claims use transportable quantities: `logFC`, `p_adj`, or re-run on a **shared gene universe**. Inspect `diagnostics["scoring"]`. |

---

## 2. p-values and cutoffs

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **Parameters named `pval_cutoff` filter adjusted p when available** | Prefer `padj_cutoff` (enrichment, `filter_active_genes`, `extract_gene_lists`). Legacy `pval_cutoff` still maps to **BH/FDR-adjusted** columns (`p_adj` / `p.adjust`), not raw `p_val`. | Always pass adjusted columns, or recompute BH yourself |
| **Heuristic filter defaults live in one code dict** | `HEURISTIC_FILTER_DEFAULTS` / `PSEUDOBULK_FILTER_DEFAULTS` set `logfc_cutoff`, residual, score, and FDR gates. Documented table: {doc}`statistical_guidance`. | Override cutoffs explicitly; do not hard-code stale numbers in papers without checking the installed version |
| **`extract_gene_lists` prefers adjusted p** | If only raw p exists, a **warning** is emitted and the cutoff is applied to raw p (inflates false positives). | Provide `p_adj` / `padj` |
| **MixedLM: `p_adj` tests `mixedlm_coef`, not sample-aware `logFC`** | Significant / s3 direction use `mixedlm_coef > 0`. Sign discordance emits a **warning** and increments `n_genes_logFC_mixedlm_sign_discordant`. | Inspect `mixedlm_coef` vs `logFC` before claiming effect direction |
| **Built-in `significant` is empty without permutation** | Default `use_permutation=False` → empty `significant` by design (needs residual FDR). | `use_permutation=True` or exploratory `filter_active_genes` |

---

## 3. Bias correction and gene features

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **Huber fit uses only `gene_length > 0` (and finite intron ≥ 0)** | Length 0 / NaN (GTF missing exons, unmapped symbols) are **not** fit covariates — avoids `log1p(0)` leverage. Expressed genes without length get median-centered residuals. | Fix GTF / feature table mapping; check `diagnostics["n_genes_with_valid_features"]` |
| **Missing GTF length is “no usable length”, not a real 0 bp gene** | Feature generator stores non-positive length as NaN. | Custom tables should use NaN for unknown length |
| **GTF exons need `transcript_id` for `intron_number`** | Slim GFF3 conversions fail with a clear error. | Full GENCODE/Ensembl GTF or prebuilt parquet |

---

## 4. Enrichment

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **GSEA needs signed ranks** | Auto-pick prefers `logFC` / t-stat-like columns; one-sided non-negative score columns are **not** auto-selected (non-negative → one-sided NES). Forcing one warns. | Pass `all_results["logFC"]` or `score_column="logFC"` |
| **Gene IDs must match gene-set universe (case included)** | Mapping rate &lt; 20% warns; 0% → empty GSEA/ORA with `reason="no_ranked_genes_mapped"`. Duplicate IDs keep the entry with largest absolute score. Enrichr is usually UPPERCASE. | `gene_case="upper"` for mouse-style symbols |
| **ORA uses the same mapping-rate gate** | Low overlap prints input vs gene-set examples. | Fix organism / symbol type / case |

---

## 5. Design, state, and packaging

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **`diagnose_design` runs inside `active_score` when `sample_col` or pseudobulk is set** | Warnings are logged and stored under `adata.uns["scatrans"]["diagnostics"]["design"]`. | Call `diagnose_design` / `recommend_workflow` yourself before scoring |
| **`copy_input=True` (default) isolates caller AnnData** | Writes to labels, layers, `.var` stay on the working copy. | `copy_input=False` only when you accept mutation |
| **Wheel installs only the importable package + data** | No tests/docs/h5ad; no shadow `tl.py`/`enrich.py`. sdist includes tests/docs for development. | Editable install from a clean checkout for development |

---

## 6. What automated tests can and cannot catch

| Layer | CI coverage (examples) | Typically needs domain review |
|-------|------------------------|--------------------------------|
| Numerics / indexing / state | Sparse vs dense, `n_jobs`, seed, order shuffle, `copy_input` | — |
| Null calibration | Type I of permutation p @ 0.05 ∈ [0.03, 0.07] | Choice of α, FDR interpretation under sparse alternatives |
| Ranking semantics | Planted up in top-N, down not; s3 gated | Whether residual-only mid-ranks are desirable for your paper |
| Huber leverage | length 0/NaN not in `n_genes_used_for_fit` | Annotation quality thresholds |
| Scale / λ | Same raw logFC maps to different soft-scale under different gene backgrounds | Cross-dataset absolute score comparisons |
| Naming / API contracts | `padj_cutoff` vs raw p warnings; Seurat `avg_log2FC` | User mental model of “pval” |
| Assumption lock-in | `tests/test_domain_assumptions_verified.py` + `tests/test_statistical_guards.py` | New product semantics not yet encoded |

Bugs that “do exactly what the code says” but violate domain intent (wrong direction gate, sentinel 0 length, silent raw-p cutoffs, treating a within-run residual magnitude as absolute) are **semantic**. They are listed above so they stay explicit.

---

## Checklist before interpreting top hits

1. Did I want **DE-significant** genes, or **nascent-excess** rank?  
2. If MixedLM: do `logFC` and `mixedlm_coef` agree in sign?  
3. Is residual high because of biology, or length/annotation holes?  
4. For enrichment: signed metric + symbol case matching?  
5. For claims: `p_adj` / `unspliced_excess_fdr`, not the nascent residual alone.
