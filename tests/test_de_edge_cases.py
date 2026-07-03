"""Numerical edge cases and cross-backend consistency for DE helpers."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

from scatrans._de import (
    _run_de_wrapper,
    _run_mixedlm_de,
    _validate_de_result,
)


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
