"""Low-cost statistical guards for CI.

These three tests are intentionally small but block entire *classes* of bugs:

1. **Null calibration** — random labels → unadjusted permutation Type I @ α=0.05
   in a narrow band (guards permutation estimator + one-sided direction + FDR plumbing).
2. **Ground-truth ranking** — planted up / down / flat genes: up in top-N, down not
   (guards active_score direction gating / composite semantics — bug class #1).
3. **Huber leverage** — gene_length 0 / NaN never enter the fit count
   (guards valid_feat / GTF-zero leverage — bug class #4).
"""

from __future__ import annotations

import anndata as ad
import numpy as np

import scatrans as scat


def _null_velocity_adata(
    n_cells: int = 80,
    n_genes: int = 100,
    *,
    seed: int = 1,
) -> ad.AnnData:
    """i.i.d. counts + random group labels (global null, no true DE / excess)."""
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(5, 0.45, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(X)
    labels = np.array(["A"] * (n_cells // 2) + ["B"] * (n_cells - n_cells // 2))
    adata.obs["condition"] = rng.permutation(labels)
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = (X * 0.4 + rng.poisson(0.3, size=X.shape)).astype(float)
    adata.var["gene_length"] = rng.integers(800, 4000, n_genes)
    adata.var["intron_number"] = rng.integers(0, 8, n_genes)
    adata.var_names = [f"g{i}" for i in range(n_genes)]
    return adata


def test_null_label_permutation_type_i_error_calibrated():
    """Random labels: permutation p-values @0.05 stay near nominal α.

    With n_perm=100 the discrete Phipson–Smyth p-values are (1+k)/101.
    Under a global null the empirical Type I rate of unadjusted p < 0.05
    should sit near 0.05 (we allow [0.03, 0.07] for finite-sample noise).

    BH-FDR under a *global* null is intentionally conservative (often ~0);
    we only assert it is not *inflated* above 0.07.
    """
    adata = _null_velocity_adata(seed=1)
    ad_out, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="A",
        reference_group="B",
        use_permutation=True,
        n_perm=100,
        n_jobs=2,
        show_plot=False,
        min_total_counts=1,
        bias_correction="none",
        ranking_mode="nascent_excess",  # residual-focused; still runs full null machinery
    )
    meta = ad_out.uns["scatrans"]
    assert meta.get("use_fdr_for_significance") is True, meta.get("perm_disabled_reason")

    for col in ("unspliced_excess_pval", "active_score_pval"):
        assert col in allr.columns
        p = allr[col].to_numpy(dtype=float)
        p = p[np.isfinite(p)]
        assert len(p) >= 50
        type_i = float(np.mean(p < 0.05))
        assert 0.03 <= type_i <= 0.07, f"{col} Type I @0.05 = {type_i:.4f} outside [0.03, 0.07]"

    for col in ("unspliced_excess_fdr", "active_score_fdr"):
        if col not in allr.columns:
            continue
        fdr = allr[col].to_numpy(dtype=float)
        fdr = fdr[np.isfinite(fdr)]
        type_i_fdr = float(np.mean(fdr < 0.05))
        assert type_i_fdr <= 0.07, f"{col} FDR Type I inflated: {type_i_fdr:.4f}"


def test_ground_truth_up_in_topn_down_not():
    """Planted DE: strong up ranks in top-N; strong down does not (bug #1 class).

    Isolates DE legs (weight_unspliced=0) so residual noise cannot rescue down genes
    via s2. Direction gating on s1/s3 must keep downregulation out of the top ranks.
    """
    rng = np.random.default_rng(42)
    n_cells, n_genes = 100, 40
    n_t = n_cells // 2
    X = rng.negative_binomial(5, 0.4, size=(n_cells, n_genes)).astype(float)

    # Gene 0: strong upregulation in target (Disease)
    X[:n_t, 0] = rng.negative_binomial(20, 0.3, size=n_t).astype(float) + 25
    X[n_t:, 0] = rng.negative_binomial(3, 0.5, size=n_cells - n_t).astype(float) + 1
    # Gene 1: strong downregulation in target
    X[:n_t, 1] = rng.negative_binomial(3, 0.5, size=n_t).astype(float) + 1
    X[n_t:, 1] = rng.negative_binomial(20, 0.3, size=n_cells - n_t).astype(float) + 25
    # Gene 2: flat (null-ish)
    X[:, 2] = rng.negative_binomial(5, 0.4, size=n_cells).astype(float) + 5

    adata = ad.AnnData(X)
    adata.obs["condition"] = ["Disease"] * n_t + ["Control"] * (n_cells - n_t)
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.35 + 0.5
    adata.var["gene_length"] = rng.integers(900, 3500, n_genes)
    adata.var["intron_number"] = rng.integers(0, 6, n_genes)
    adata.var_names = [f"G{i}" for i in range(n_genes)]

    _, _, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        weight_fc=1.0,
        weight_unspliced=0.0,
        weight_pval=1.0,
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
        bias_correction="none",
    )

    assert float(allr.loc["G0", "logFC"]) > 0.5
    assert float(allr.loc["G1", "logFC"]) < -0.5

    top_n = 8
    top = list(allr.nlargest(top_n, "active_score").index)
    assert "G0" in top, f"planted up-gene G0 not in top-{top_n}: {top}"
    assert "G1" not in top, f"planted down-gene G1 incorrectly in top-{top_n}: {top}"
    # Down gene must not outrank the planted up gene on composite score
    assert float(allr.loc["G1", "active_score"]) < float(allr.loc["G0", "active_score"])


def test_huber_excludes_zero_and_nan_length_from_fit():
    """gene_length 0 / NaN must not count toward n_genes_used_for_fit (bug #4 class)."""
    from scatrans._utils import _fit_huber_bias_correction

    rng = np.random.default_rng(0)
    n_good, n_zero, n_nan = 60, 25, 15
    n = n_good + n_zero + n_nan
    gene_length = np.concatenate(
        [
            rng.uniform(600, 4000, n_good),
            np.zeros(n_zero),
            np.full(n_nan, np.nan),
        ]
    )
    intron_number = rng.integers(0, 10, n).astype(float)
    # Put wild y on bad lengths so a leak into the fit would bias the slope hard
    delta = np.zeros(n)
    delta[:n_good] = 0.04 * np.log1p(gene_length[:n_good]) + rng.normal(0, 0.02, n_good)
    delta[n_good:] = rng.normal(40.0, 5.0, n_zero + n_nan)

    # Same rule as active_score valid_feat
    valid_feat = (
        np.isfinite(gene_length)
        & np.isfinite(intron_number)
        & (gene_length > 0)
        & (intron_number >= 0)
    )
    assert int(valid_feat.sum()) == n_good
    assert not valid_feat[n_good:].any()

    valid_expr = np.ones(n, dtype=bool)
    X_features = np.column_stack(
        [np.log1p(gene_length[valid_feat]), np.log1p(intron_number[valid_feat])]
    )
    _residual, info = _fit_huber_bias_correction(
        delta,
        gene_length,
        intron_number,
        total_us_for_weights=np.full(n, 100.0),
        valid_feat=valid_feat,
        valid_expr=valid_expr,
        X_features=X_features,
        bias_correction="huber_length_intron",
        min_fit_obs=30,
    )
    assert info.get("bias_corrected") is True
    n_used = int(info["n_genes_used_for_fit"])
    assert n_used == n_good, f"fit used {n_used} genes; expected only {n_good} with length>0"
    # Sanity: recovered length coefficient still near planted trend
    coef = float(info["coef_gene_length"])
    assert abs(coef - 0.04) < 0.025


def test_active_score_lambda_is_data_adaptive_not_cross_run_comparable():
    """Same gene raw signals → different active_score when gene background changes.

    Soft-scale λ = median(positive)/ln2 from the gene vector in that run. Enlarging
    the gene set with many large logFC values inflates λ_fc and rescales scores.
    """
    from scatrans._utils import _get_exponential_scale_lambda, _soft_scale

    # Shared "biology" for gene 0
    logfc_core = np.array([1.0, 0.5, 0.2, 0.1])
    # Background A: modest positives → smaller λ
    bg_a = np.concatenate([logfc_core, np.full(20, 0.15)])
    # Background B: many large positives → larger λ → same logFC=1 maps lower soft-scale
    bg_b = np.concatenate([logfc_core, np.full(20, 3.0)])
    lam_a = _get_exponential_scale_lambda(bg_a)
    lam_b = _get_exponential_scale_lambda(bg_b)
    assert lam_b > lam_a * 1.5
    s_a = float(_soft_scale(np.array([1.0]), lam_a)[0])
    s_b = float(_soft_scale(np.array([1.0]), lam_b)[0])
    assert s_a > s_b + 0.05  # same raw logFC, different score contribution

    # End-to-end: diagnostics expose lambdas
    rng = np.random.default_rng(0)
    n_cells, n_genes = 60, 40
    X = rng.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(X)
    adata.obs["condition"] = ["Disease"] * 30 + ["Control"] * 30
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.4
    adata.var["gene_length"] = rng.integers(800, 3000, n_genes)
    adata.var["intron_number"] = rng.integers(0, 5, n_genes)
    ad_out, _, _ = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
        bias_correction="none",
    )
    scoring = ad_out.uns["scatrans"]["diagnostics"]["scoring"]
    assert "lambda_fc" in scoring and scoring["lambda_fc"] > 0
    assert "scale_note" in scoring and "within" in scoring["scale_note"].lower()


def test_active_score_records_valid_feat_excluding_bad_lengths():
    """End-to-end: active_score diagnostics n_genes_with_valid_features ignores length≤0/NaN."""
    rng = np.random.default_rng(3)
    n_cells, n_genes = 60, 30
    X = rng.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(X)
    adata.obs["condition"] = ["Disease"] * 30 + ["Control"] * 30
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.4
    gl = rng.integers(800, 3000, n_genes).astype(float)
    gl[0] = 0.0
    gl[1] = np.nan
    gl[2] = -1.0  # non-positive → excluded like 0
    adata.var["gene_length"] = gl
    adata.var["intron_number"] = rng.integers(0, 5, n_genes)
    adata.var_names = [f"g{i}" for i in range(n_genes)]

    ad_out, _, _ = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=False,
        show_plot=False,
        n_jobs=1,
        min_total_counts=1,
    )
    n_valid = int(ad_out.uns["scatrans"]["diagnostics"]["n_genes_with_valid_features"])
    expected = int(np.sum(np.isfinite(gl) & (gl > 0)))
    assert n_valid == expected
    assert n_valid == n_genes - 3
