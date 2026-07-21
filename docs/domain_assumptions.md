# Domain Assumptions

Domain conventions enforced by the implementation. Use this page to check
expectations against behavior.

Related: {doc}`statistical_guidance` (reporting), {doc}`faq` (scope).

---

## 1. What “active” means

| Assumption | Implication | If you need something else |
|------------|-------------|----------------------------|
| **Active transcription ranking is upregulation-oriented** | Built-in `significant` requires positive `logFC` (and positive residual). Default `filter_active_genes` direction is `"up"`. | `logfc_direction="down"` / `"both"` on pure DE tables; or use `differential_expression` alone |
| **Composite `active_score` DE legs (s1 logFC, s3 −log p) are upregulation-gated** | Negative logFC → s1 = 0 and s3 = 0 (`mixedlm_coef > 0` when MixedLM). Strongly down genes do not get mid-rank scores from p-values alone. | Pass residual-only mode: `ranking_mode="nascent_excess"` |
| **Residual leg s2 is independent of DE direction and DE significance** | Positive unspliced excess can score even if DE is weak, untested, or `p_adj` filled to 1 (e.g. PyDESeq2 independent filtering). **Top-N by `active_score` is not a DE-significant list.** | `filter_active_genes(..., select_by="de")` or explicit `padj_cutoff` / `logfc_cutoff`, or `significant` conjunction |
| **`select_by="de"` decouples membership from the proxy** | DE gates define the list; nascent / composite gates are skipped; sorting is by `p_adj` then `logFC`. | Default `select_by="composite"` keeps prior proxy-aware filtering |
| **`qc.regime_diagnosis` reliability is U-shaped in unspliced fraction** | Full reliability in a normal band (~10–45%); degrades at low (noise) and high (nuclear/gDNA, gamma mis-fit) extremes. Scales `mechanism_confidence` in `partition_de_by_mechanism` and when `annotate_mechanism=True`. | This is data-quality only — not dynamic vs steady-state |
| **Detection ≠ mechanism** | `nascent_poisson_z` / `add_nascent_score=True` is an absolute nascent-increase **detection** score (induction-coupled). Mechanism labels always use the induction-normalized residual. | Do not pass `nascent_poisson_z` as `residual_col` for `annotate_mechanism_class` if you want transcription-vs-stabilization |
| **Residual is one-sided (positive excess)** | Negative excess does not contribute to s2 soft-scale. Permutation FDR on residual is one-sided for positive excess. | Do not interpret low residual FDR as “repression” |
| **`active_score` 0–100 is within-run relative (λ is data-adaptive)** | Each leg uses soft-scale `1 − exp(−x/λ)` with `λ ≈ median(positive x) / ln(2)` estimated **from the genes in that run** (plus floors). The same gene with the same raw logFC/residual/p can map to **different** scores if the background gene set or subset changes (other genes shift the median → λ → all soft-scaled values). | Compare **ranks within one `active_score` call** only. For cross-subset or cross-dataset claims use transportable quantities: `logFC`, `p_adj`, residual magnitude, or re-run on a **shared gene universe** with the understanding scores are still not absolute. Inspect `diagnostics["scoring"]` (`lambda_fc`, `lambda_res`, `lambda_pval`). |

### High-risk misuse of the 0–100 scale

These are **incorrect** uses of `active_score` numbers (code will not stop you):

1. **“Gene X scores 80 in B cells but 50 in T cells”** as absolute activity across subsets (each subset re-estimates λ).
2. **Cross-dataset / cross-experiment numeric comparison** of scores.
3. **HVG filter before/after** comparison of the same gene’s score (gene set change → λ change → rescaling).
4. Treating score cutoffs (e.g. 55) as universal thresholds across designs (heuristic presets assume a typical single-run scale).

Safe: rank genes **within** one contrast + one gene universe; report DE/`unspliced_excess_fdr` for claims.

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
| **GSEA needs signed ranks** | Auto-pick prefers `logFC` / t-stat-like columns; `active_score` is **not** auto-selected (non-negative → one-sided NES). Forcing it warns. | Pass `all_results["logFC"]` or `score_column="logFC"` |
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

Bugs that “do exactly what the code says” but violate domain intent (wrong direction gate, sentinel 0 length, silent raw-p cutoffs, treating 0–100 as absolute) are **semantic**. They are listed above so they stay explicit.

---

## Checklist before interpreting top hits

1. Did I want **DE-significant** genes, or **composite / nascent-excess** rank?  
2. Am I treating 0–100 as absolute across subsets/datasets? (**Don't** — λ is data-adaptive.)  
3. If MixedLM: do `logFC` and `mixedlm_coef` agree in sign?  
4. Is residual high because of biology, or length/annotation holes?  
5. For enrichment: signed metric + symbol case matching?  
6. For claims: `p_adj` / `unspliced_excess_fdr`, not `active_score` alone.
