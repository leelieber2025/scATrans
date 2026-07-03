"""Tests for DE backends and raw-count utilities (boost _de / _utils / tl coverage)."""

import importlib.util

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scatrans as scat
from scatrans._de import _run_mixedlm_de


@pytest.mark.slow
def test_mixedlm_logfc_respects_target_reference_order():
    """Target alphabetically before reference must still yield target-minus-reference logFC."""
    np.random.seed(0)
    n_cells, n_genes = 40, 5
    X = np.zeros((n_cells, n_genes))
    X[:20] = 100.0
    X[20:] = 1.0
    X += np.random.randn(n_cells, n_genes) * 5.0
    X = np.maximum(X, 0.0)
    obs = pd.DataFrame(
        {
            "condition": ["A"] * 20 + ["Z"] * 20,
            "sample": (["s1"] * 10 + ["s2"] * 10) * 2,
        }
    )
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]),
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    a_vs_z = _run_mixedlm_de(
        adata,
        groupby="condition",
        target_group="A",
        reference_group="Z",
        sample_col="sample",
        n_jobs=1,
    )
    z_vs_a = _run_mixedlm_de(
        adata,
        groupby="condition",
        target_group="Z",
        reference_group="A",
        sample_col="sample",
        n_jobs=1,
    )
    # Exclude degenerate fits (logFC=0, p_val=1) — MixedLM can fail on ~constant genes.
    fitted_a = a_vs_z["p_val"] < 1.0
    fitted_z = z_vs_a["p_val"] < 1.0
    assert fitted_a.sum() >= 4
    assert fitted_z.sum() >= 4
    assert np.all(a_vs_z.loc[fitted_a, "logFC"] > 0)
    assert np.all(z_vs_a.loc[fitted_z, "logFC"] < 0)
    assert np.allclose(a_vs_z["logFC"], -z_vs_a["logFC"])


@pytest.mark.slow
def test_differential_expression_mixed_model(adata_mixed_small):
    ad, res = scat.differential_expression(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_mixed_model=True,
        sample_col="sample",
        n_jobs=1,
    )
    assert "delta_variance" in res.columns


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
@pytest.mark.slow
def test_differential_expression_pseudobulk_pydeseq2(adata_pb):
    scat.store_raw_counts(adata_pb, layer="counts")
    ad, res = scat.differential_expression(
        adata_pb,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        sample_col="sample",
        pseudobulk_de_backend="pydeseq2",
        de_preprocess="none",
    )
    assert "logFC" in res.columns
    assert ad.uns["scatrans"]["use_pseudobulk"] is True


@pytest.mark.slow
def test_pseudobulk_categorical_groupby_column(adata_pb):
    """Regression: Categorical obs columns must not break pb_key string concat."""
    ad = adata_pb.copy()
    ad.obs["sample"] = pd.Categorical(ad.obs["sample"].astype(str))
    ad.obs["condition"] = pd.Categorical(ad.obs["condition"].astype(str))
    _, res = scat.differential_expression(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        sample_col="sample",
        pseudobulk_de_backend="scanpy",
        de_method="wilcoxon",
        min_cells=1,
        min_counts=1,
    )
    assert len(res) == adata_pb.n_vars


def test_differential_expression_pseudobulk_scanpy(adata_pb):
    ad, res = scat.differential_expression(
        adata_pb,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        sample_col="sample",
        pseudobulk_de_backend="scanpy",
        de_method="wilcoxon",
    )
    assert len(res) == adata_pb.n_vars


def test_store_raw_counts_mature_nascent_layers():
    np.random.seed(3)
    X = np.random.negative_binomial(4, 0.4, size=(30, 40)).astype(float)
    ad = sc.AnnData(X)
    ad.layers["mature"] = X.copy()
    ad.layers["nascent"] = X * 0.4
    scat.store_raw_counts(ad, layer="counts")
    assert "raw_mature" in ad.layers or "raw_spliced" in ad.layers


def test_active_score_subset_and_gene_type(adata_basic):
    ad = adata_basic.copy()
    ad.var["gene_type"] = "protein_coding"
    _, _, allr = scat.active_score(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        subset_col="sample",
        subset_values=["s0", "s1", "s2", "s3"],
        gene_type_filter="protein_coding",
        use_permutation=False,
        show_plot=False,
    )
    assert len(allr) > 0


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
@pytest.mark.slow
def test_active_score_pseudobulk_pydeseq2(adata_pb):
    scat.store_raw_counts(adata_pb, layer="counts")
    res, _, allr = scat.active_score(
        adata_pb,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        sample_col="sample",
        de_preprocess="none",
        use_permutation=False,
        show_plot=False,
    )
    assert "active_score" in allr.columns
    assert res.uns["scatrans"]["use_pseudobulk"] is True


@pytest.mark.plot
def test_active_score_show_plot_comet(adata_basic):
    """Exercise show_plot=True path (Agg backend)."""
    import matplotlib.pyplot as plt

    _, _, _ = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=False,
        show_plot=True,
    )
    plt.close("all")


@pytest.mark.slow
def test_run_go_bp_smoke():
    res = scat.run_go(
        ["FakeGene"],
        ontology="BP",
        organism="mouse",
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    assert isinstance(res, pd.DataFrame)


def test_wilcoxon_logfc_is_log2_scale_no_secondary_div():
    """Regression test for wilcoxon logFC handling.

    scanpy always returns logfoldchanges on log2( (expm1(m_t)+e)/(expm1(m_r)+e) ) scale
    independent of the rank test method. The old secondary raw_lfc / np.log(2) for
    wilcoxon was incorrect and systematically shrunk values.
    """
    np.random.seed(0)
    # Construct data with a clear fold difference
    n = 80
    X = np.random.poisson(lam=2, size=(n, 3)).astype(float)
    X[: n // 2, 0] += 30  # strong up in first group for gene 0
    ad = sc.AnnData(X)
    ad.obs["group"] = ["A"] * (n // 2) + ["B"] * (n // 2)
    ad.var_names = ["g0", "g1", "g2"]
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # Run via the internal wrapper using scanpy path with wilcoxon and t-test
    from scatrans._de import _run_de_wrapper

    res_w = _run_de_wrapper(
        ad, groupby="group", target_group="A", reference_group="B", de_method="wilcoxon"
    )
    res_t = _run_de_wrapper(
        ad, groupby="group", target_group="A", reference_group="B", de_method="t-test"
    )

    lfc_w = res_w.loc["g0", "logFC"]
    lfc_t = res_t.loc["g0", "logFC"]
    # They must match (within float tol); old code would have made wilcox ~1.44x smaller
    assert np.allclose(lfc_w, lfc_t, rtol=1e-5, atol=1e-5)
    # Sanity: positive and reasonably large
    assert lfc_w > 1.0
    # Ensure no one is doing natural-log scale either (log2( (30+2)/2 ) ~ log2(16)~4, actual will be lower due to norm)
    # Just check not hugely inflated or deflated. The equality is the key assertion.


def test_validation_errors(adata_basic):
    with pytest.raises(ValueError, match="must be different"):
        scat.active_score(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Disease",
            show_plot=False,
        )
    with pytest.raises(ValueError, match="sample_col"):
        scat.differential_expression(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_mixed_model=True,
        )


def test_inactive_backend_params_not_validated(adata_basic):
    """Memento / mixed-model options are only validated when the feature is enabled."""
    # Invalid memento settings must not block default heuristic path
    scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_memento_de=False,
        memento_num_boot=50,
        memento_capture_rate=1.5,
        show_plot=False,
        use_permutation=False,
    )
    _, _ = scat.differential_expression(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_memento_de=False,
        memento_capture_rate=1.5,
    )
    # Invalid mixed_model_pval is ignored unless use_mixed_model=True
    scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_mixed_model=False,
        mixed_model_pval="invalid",
        show_plot=False,
        use_permutation=False,
    )
