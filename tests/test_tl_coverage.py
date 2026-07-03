"""Tests for tl module: DE, raw counts, design diagnosis, permutation backends."""

import warnings

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scatrans as scat


def test_differential_expression_delta_variance_warning(adata_de_only, caplog):
    import logging

    with caplog.at_level(logging.WARNING):
        scat.differential_expression(
            adata_de_only,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            de_method="wilcoxon",
            use_delta_variance_pval=True,
            delta_var_pval_cutoff=0.01,
        )
    assert any("use_delta_variance_pval=True is not enforced" in r.message for r in caplog.records)


def test_differential_expression_basic(adata_de_only):
    ad, res = scat.differential_expression(
        adata_de_only,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        de_method="wilcoxon",
    )
    assert "logFC" in res.columns
    assert "p_adj" in res.columns
    assert ad.uns.get("scatrans", {}).get("mode") == "differential_expression"


def test_differential_expression_subset(adata_de_only):
    ad = adata_de_only.copy()
    # Keep both conditions within the subset
    ad.obs["batch"] = ["B1" if i % 2 == 0 else "B2" for i in range(ad.n_obs)]
    _, res = scat.differential_expression(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        subset_col="batch",
        subset_values=["B1", "B2"],
        de_method="wilcoxon",
    )
    assert len(res) == ad.n_vars


def test_diagnose_design_basic(adata_basic):
    diag = scat.diagnose_design(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        sample_col="sample",
    )
    assert diag["n_cells_target"] == 60
    assert diag["n_cells_reference"] == 60
    assert diag["n_samples_target"] >= 3
    assert "warnings" in diag
    assert diag["suggested_preset"] in ("heuristic", "pseudobulk", None)


def test_diagnose_design_small_cells_warning(adata_basic):
    ad = adata_basic[:30].copy()
    ad.obs["condition"] = ["Disease"] * 15 + ["Control"] * 15
    diag = scat.diagnose_design(ad, "condition", "Disease", "Control")
    assert any("small" in w.lower() for w in diag["warnings"])


def test_store_and_restore_raw_counts(adata_de_only):
    ad = adata_de_only.copy()
    X_raw = np.asarray(ad.X, dtype=float).copy()
    scat.store_raw_counts(ad, layer="counts", save_raw=False)
    assert "counts" in ad.layers
    assert "raw_gene_list" in ad.uns.get("scatrans", {})
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    restored = scat.restore_raw_counts(ad, layer="counts", inplace=False)
    assert np.allclose(np.asarray(restored.X), X_raw)
    assert "log1p" not in restored.uns


def test_reconcile_log1p_marker_detects_stale_metadata(adata_de_only):
    from scatrans._utils import _reconcile_log1p_marker, _x_looks_log_normalized

    ad = adata_de_only.copy()
    scat.store_raw_counts(ad, layer="counts", save_raw=False)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.X = ad.layers["counts"].copy()
    assert "log1p" in ad.uns
    assert not _x_looks_log_normalized(ad.X)

    assert _reconcile_log1p_marker(ad) is False
    assert "log1p" not in ad.uns


def test_de_preprocess_auto_after_restore_raw_counts(adata_de_only):
    ad = adata_de_only.copy()
    scat.store_raw_counts(ad, layer="counts", save_raw=False)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    restored = scat.restore_raw_counts(ad, layer="counts", inplace=False)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        _, res = scat.differential_expression(
            restored,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            de_method="t-test_overestim_var",
            de_preprocess="auto",
        )

    raw_count_msgs = [
        w.message
        for w in caught
        if "raw count" in str(w.message).lower() or "logarithmize" in str(w.message).lower()
    ]
    assert not raw_count_msgs
    assert "logFC" in res.columns


def test_de_preprocess_auto_strips_stale_log1p_marker_inplace(adata_de_only):
    ad = adata_de_only.copy()
    scat.store_raw_counts(ad, layer="counts", save_raw=False)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    ad.X = ad.layers["counts"].copy()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        scat.differential_expression(
            ad,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            de_method="t-test_overestim_var",
            de_preprocess="auto",
            copy_input=False,
        )

    raw_count_msgs = [
        w.message
        for w in caught
        if "raw count" in str(w.message).lower() or "logarithmize" in str(w.message).lower()
    ]
    assert not raw_count_msgs
    assert "log1p" in ad.uns


def test_ensure_raw_counts_from_raw_attr(adata_de_only):
    ad = adata_de_only.copy()
    ad.raw = ad.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    scat.ensure_raw_counts(ad)
    assert "counts" in ad.layers


@pytest.mark.slow
def test_perm_de_backend_default_is_same(adata_basic):
    res, _, _ = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=True,
        n_perm=4,
        de_method="wilcoxon",
        n_jobs=1,
        show_plot=False,
    )
    meta = res.uns["scatrans"]
    assert meta.get("perm_de_backend") == "same"


@pytest.mark.slow
def test_perm_de_backend_fast_still_works(adata_basic):
    res, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=True,
        perm_de_backend="fast",
        n_perm=4,
        n_jobs=1,
        show_plot=False,
    )
    assert res.uns["scatrans"].get("perm_de_backend") == "fast"
    assert "unspliced_excess_fdr" in allr.columns


def test_active_score_show_effective_gamma(adata_basic):
    res, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        show_effective_gamma=True,
        use_permutation=False,
        show_plot=False,
    )
    assert "effective_gamma" in allr.columns


def test_active_score_bias_correction_none(adata_basic):
    res, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        bias_correction="none",
        use_permutation=False,
        show_plot=False,
    )
    assert "unspliced_excess_residual" in allr.columns


def test_active_score_gamma_robust_median(adata_basic):
    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="robust_median",
        use_permutation=False,
        show_plot=False,
    )
    assert "unspliced_excess_residual" in allr.columns


def test_active_score_gamma_empirical_bayes(adata_basic):
    adata_res, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="empirical_bayes",
        use_permutation=False,
        show_plot=False,
    )
    vel = adata_res.uns["scatrans"]["diagnostics"]["velocity"]
    assert vel["gamma_method"] == "empirical_bayes"
    assert "gamma_prior_mean" in vel
    assert "shrinkage_summary" in vel
    assert "gamma_shrinkage_weight" in allr.columns


@pytest.mark.slow
def test_gamma_methods_ranking_correlation(adata_basic):
    from scipy.stats import spearmanr

    _, _, allr_heur = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="heuristic_shrink",
        use_permutation=False,
        show_plot=False,
    )
    _, _, allr_eb = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="empirical_bayes",
        use_permutation=False,
        show_plot=False,
    )
    merged = allr_heur[["active_score", "unspliced_excess_residual"]].join(
        allr_eb[["active_score", "unspliced_excess_residual"]],
        lsuffix="_heur",
        rsuffix="_eb",
    )
    rho_score, _ = spearmanr(merged["active_score_heur"], merged["active_score_eb"])
    rho_res, _ = spearmanr(
        merged["unspliced_excess_residual_heur"], merged["unspliced_excess_residual_eb"]
    )
    assert rho_score > 0.5
    assert rho_res > 0.5


def test_empirical_bayes_small_reference(adata_small_reference):
    adata_res, _, allr = scat.active_score(
        adata_small_reference,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="empirical_bayes",
        use_permutation=False,
        show_plot=False,
    )
    vel = adata_res.uns["scatrans"]["diagnostics"]["velocity"]
    assert vel["n_genes_used_for_prior"] >= 10
    assert allr["gamma_shrinkage_weight"].between(0, 1).all()


@pytest.mark.plot
def test_gamma_shrinkage_plot(adata_basic):
    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        gamma_method="empirical_bayes",
        use_permutation=False,
        show_plot=False,
    )
    fig, ax = scat.pl.gamma_shrinkage_plot(allr, show=False)
    assert fig is not None and ax is not None


def test_active_score_prioritize_velocity(adata_basic):
    with pytest.warns(DeprecationWarning, match="prioritize_velocity"):
        _, _, allr = scat.active_score(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            prioritize_velocity=True,
            use_permutation=False,
            show_plot=False,
        )
    assert "active_score" in allr.columns


def test_filter_active_genes_return_mask_and_inplace(adata_mixed_small):
    _, _, allr = scat.active_score_simple(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        show_plot=False,
    )
    mask = scat.filter_active_genes(allr, preset="permissive", return_mask=True)
    assert isinstance(mask, pd.Series)
    copy_df = allr.copy()
    out = scat.filter_active_genes(copy_df, preset="permissive", inplace=True)
    assert out is copy_df
    assert len(copy_df) <= len(allr)


def test_differential_expression_pseudobulk_on_de_only(adata_de_only):
    """Critical path: use_pseudobulk=True must not require spliced/unspliced layers (Bug 1 fix)."""
    ad, res = scat.differential_expression(
        adata_de_only,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        pseudobulk_de_backend="scanpy",
        sample_col="sample",
        de_method="wilcoxon",
        min_cells=1,
        min_counts=1,
    )
    assert "logFC" in res.columns
    assert "p_adj" in res.columns
    assert ad.n_obs < adata_de_only.n_obs  # pseudobulked (after filtering)
    # When no velocity layers, use_total_for_x should have been ignored and X used for agg
    assert "pb_x_source" in ad.obs.columns


def test_prepare_log_normalized_expression_gap(adata_de_only):
    """Cover the former 5<mx<=20 gap path (Bug 3). Data after mild log-like transform should not be re-log1p'ed."""
    from scatrans._utils import _prepare_log_normalized_expression

    ad = adata_de_only.copy()
    # Simulate "already somewhat logged" data with max ~10 (scran-like or light log)
    rng = np.random.default_rng(42)
    mat = rng.normal(loc=2.0, scale=1.5, size=ad.shape)
    mat = np.clip(mat, -1, 12)  # ensure some range covering gap, allow mild neg for realism
    ad.X = mat
    # No log1p in uns, not integer counts
    X_out = _prepare_log_normalized_expression(ad)
    # Should return roughly same scale (no additional log1p applied)
    assert X_out.max() < 15
    assert np.allclose(X_out.mean(), mat.mean(), atol=1.0)
