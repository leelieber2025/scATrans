"""Numerical edge cases and cross-backend consistency for DE helpers."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from statsmodels.stats.multitest import multipletests

from scatrans._de import (
    _run_de_wrapper,
    _run_memento_de,
    _run_mixedlm_de,
    _validate_de_result,
)
from scatrans._permutation import run_permutation_test


def test_validate_de_result_raises_on_missing_columns():
    bad = pd.DataFrame({"logFC": [0.1], "p_val": [0.05]}, index=["G1"])
    with pytest.raises(RuntimeError, match="missing columns"):
        _validate_de_result(bad, backend="test")


def test_validate_de_result_raises_on_all_nan_column():
    bad = pd.DataFrame(
        {"logFC": [np.nan], "p_val": [0.05], "p_adj": [0.1]},
        index=["G1"],
    )
    with pytest.raises(RuntimeError, match="no finite values"):
        _validate_de_result(bad, backend="test")


def test_scanpy_backends_agree_direction_on_toy():
    """Wilcoxon and t-test should agree on up/down sign for a clear DE gene."""
    np.random.seed(1)
    n = 60
    X = np.random.poisson(lam=2, size=(n, 4)).astype(float)
    X[: n // 2, 0] += 40
    adata = ad.AnnData(X, obs=pd.DataFrame({"group": ["T"] * (n // 2) + ["R"] * (n // 2)}))
    adata.var_names = [f"g{i}" for i in range(4)]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    res_w = _run_de_wrapper(
        adata, groupby="group", target_group="T", reference_group="R", de_method="wilcoxon"
    )
    res_t = _run_de_wrapper(
        adata, groupby="group", target_group="T", reference_group="R", de_method="t-test"
    )
    assert res_w.loc["g0", "logFC"] > 0
    assert res_t.loc["g0", "logFC"] > 0
    assert ((res_w["logFC"] > 0) == (res_t["logFC"] > 0)).all()


@pytest.mark.slow
def test_mixedlm_constant_gene_counted_as_failed():
    """Near-constant genes must increment n_genes_failed_fit, not silently pass as significant."""
    np.random.seed(2)
    n_cells = 24
    X = np.random.exponential(1.0, size=(n_cells, 3))
    X[:, 2] = 5.0  # constant column
    obs = pd.DataFrame(
        {
            "condition": ["A"] * 12 + ["B"] * 12,
            "sample": (["s1"] * 6 + ["s2"] * 6) * 2,
        }
    )
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=["v0", "v1", "const"]))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    res = _run_mixedlm_de(
        adata,
        groupby="condition",
        target_group="A",
        reference_group="B",
        sample_col="sample",
        n_jobs=1,
    )
    assert res.attrs.get("n_genes_failed_fit", 0) >= 1
    assert res.attrs.get("failed_fit_rate", 0) > 0
    assert res.loc["const", "p_val"] == 1.0
    assert res.loc["const", "logFC"] == 0.0


def test_de_wrapper_finite_output_on_all_zero_gene():
    """All-zero expression should return finite schema columns, not NaN explosion."""
    np.random.seed(3)
    X = np.random.poisson(3, size=(40, 3)).astype(float)
    X[:, 1] = 0.0
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"group": ["T"] * 20 + ["R"] * 20}),
        var=pd.DataFrame(index=["up", "zero", "other"]),
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    res = _run_de_wrapper(
        adata, groupby="group", target_group="T", reference_group="R", de_method="wilcoxon"
    )
    for col in ("logFC", "p_val", "p_adj"):
        assert col in res.columns
        assert np.isfinite(res[col].to_numpy()).all()


def test_memento_counts_ann_data_obs_names_alignment():
    """counts=<AnnData> must align by obs_names, not positional index."""
    pytest.importorskip("memento")
    np.random.seed(4)
    n = 24
    X = np.random.poisson(8, size=(n, 6)).astype(float)
    obs = pd.DataFrame({"condition": ["T"] * 12 + ["R"] * 12})
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(6)]))
    counts_ad = ad.AnnData(X, obs=obs.copy(), var=adata.var.copy())
    counts_ad.obs_names = counts_ad.obs_names[::-1]

    res = _run_memento_de(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        counts=counts_ad,
        num_boot=50,
        n_cpus=1,
    )
    assert "logFC" in res.columns
    assert len(res) == adata.n_vars


def test_memento_counts_ann_data_missing_obs_raises():
    pytest.importorskip("memento")
    np.random.seed(5)
    n = 20
    X = np.random.poisson(4, size=(n, 4)).astype(float)
    obs = pd.DataFrame({"condition": ["T"] * 10 + ["R"] * 10})
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(4)]))
    partial = ad.AnnData(X[:8], obs=obs.iloc[:8].copy(), var=adata.var.copy())

    with pytest.raises(ValueError, match="missing"):
        _run_memento_de(
            adata,
            groupby="condition",
            target_group="T",
            reference_group="R",
            counts=partial,
            num_boot=20,
            n_cpus=1,
        )


def test_run_permutation_test_use_fdr_based_on_n_success(monkeypatch):
    """use_fdr must reflect successful shuffles, not requested n_perm."""
    import scatrans._permutation as perm_mod

    np.random.seed(6)
    n_cells, n_genes = 30, 8
    X = np.random.poisson(3, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"group": ["T"] * 15 + ["R"] * 15}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )
    adata.var["gene_length"] = np.linspace(800, 4000, n_genes)
    adata.var["intron_number"] = np.arange(n_genes)
    uns = X.copy()
    spl = X * 0.5
    real_score = np.linspace(10, 90, n_genes)
    real_residual = np.linspace(0.1, 2.0, n_genes)
    valid_feat = np.ones(n_genes, dtype=bool)
    valid_expr = np.ones(n_genes, dtype=bool)

    orig_task = perm_mod._single_permutation_task
    call_n = {"i": 0}

    def flaky_task(*args, **kwargs):
        call_n["i"] += 1
        if call_n["i"] <= 80:
            return np.full(n_genes, np.nan), np.full(n_genes, np.nan)
        return orig_task(*args, **kwargs)

    monkeypatch.setattr(perm_mod, "_single_permutation_task", flaky_task)

    *_, use_fdr, reason = run_permutation_test(
        n_perm=150,
        effective_n_jobs=1,
        random_seed=0,
        obs_labels=adata.obs["group"].to_numpy(),
        target_group="T",
        reference_group="R",
        adata=adata,
        X_features=None,
        valid_feat=valid_feat,
        velocity_layer_for_perm_uns=uns,
        velocity_layer_for_perm_spl=spl,
        total_us_raw=uns.sum(axis=0) + spl.sum(axis=0),
        min_total_counts=1,
        weight_fc=1.0,
        weight_unspliced=1.0,
        weight_pval=1.0,
        lambda_fc=1.0,
        lambda_res=1.0,
        lambda_pval=1.0,
        is_pseudobulk=False,
        perm_pb_backend="scanpy",
        perm_de_method="t-test_overestim_var",
        prior_weight=5.0,
        gamma_method="heuristic_shrink",
        de_preprocess="none",
        strict_pydeseq2_counts=True,
        real_score=real_score,
        real_residual=real_residual,
        valid_expr=valid_expr,
    )
    assert use_fdr is False
    assert reason == "small_permutation_space"


def test_memento_bh_excludes_nan_pvals_from_denominator():
    """Defensive: NaN p-values must not inflate BH denominator if they ever appear."""
    pvals_raw = pd.Series([0.01, 0.04, np.nan, 0.2], index=["g0", "g1", "g2", "g3"])
    valid = pvals_raw.notna()
    p_adj_clean = pd.Series(1.0, index=pvals_raw.index)
    if valid.sum() > 0:
        p_adj_clean.loc[valid] = multipletests(pvals_raw[valid].values, method="fdr_bh")[1]
    p_adj_polluted = multipletests(pvals_raw.fillna(1.0).values, method="fdr_bh")[1]
    assert p_adj_clean.loc["g0"] < p_adj_polluted[0]
