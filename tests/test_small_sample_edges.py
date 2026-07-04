"""Extreme small-sample and zero-count boundary tests.

Guards against divide-by-zero / empty-array crashes on toy designs that users
sometimes hit during subsetting or QC (1–2 cells per group, all-zero genes).
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import scatrans as scat


def _tiny_velocity_adata(
    n_target: int,
    n_ref: int,
    *,
    X: np.ndarray | None = None,
    var_names: list[str] | None = None,
    sample_col: str | None = None,
) -> ad.AnnData:
    """Minimal AnnData with spliced/unspliced layers and gene features."""
    n_cells = n_target + n_ref
    n_genes = X.shape[1] if X is not None else 3
    if X is None:
        X = np.random.poisson(3, size=(n_cells, n_genes)).astype(float)
    if var_names is None:
        var_names = [f"g{i}" for i in range(n_genes)]
    obs: dict[str, list] = {
        "condition": ["T"] * n_target + ["R"] * n_ref,
    }
    if sample_col:
        obs[sample_col] = [f"s{i}" for i in range(n_cells)]
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(obs),
        var=pd.DataFrame(index=var_names),
    )
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.5
    adata.var["gene_length"] = np.linspace(800, 4000, n_genes)
    adata.var["intron_number"] = np.arange(n_genes)
    return adata


@pytest.mark.parametrize(
    "runner",
    [
        "active_score",
        "differential_expression",
    ],
)
def test_one_cell_per_group_raises_clear_error(runner):
    """Single cell per group cannot support scanpy DE — fail fast with guidance."""
    adata = _tiny_velocity_adata(1, 1, X=np.array([[5.0, 0.0, 2.0], [1.0, 0.0, 3.0]]))
    with pytest.raises(ValueError, match="at least 2 cells per group"):
        if runner == "active_score":
            scat.active_score(
                adata,
                groupby="condition",
                target_group="T",
                reference_group="R",
                use_permutation=False,
                min_cells=1,
                min_total_counts=0,
                show_plot=False,
            )
        else:
            scat.differential_expression(
                adata,
                groupby="condition",
                target_group="T",
                reference_group="R",
                min_cells=1,
            )


def test_two_cells_per_group_all_zero_gene_returns_finite_outputs():
    """All-zero gene column must not produce NaN/inf in the full active_score table."""
    X = np.array([[5.0, 0.0, 2.0], [6.0, 0.0, 1.0], [1.0, 0.0, 3.0], [2.0, 0.0, 4.0]])
    adata = _tiny_velocity_adata(2, 2, X=X, var_names=["expressed", "all_zero", "other"])

    _, significant, all_results = scat.active_score(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        use_permutation=False,
        min_cells=1,
        min_total_counts=0,
        show_plot=False,
    )

    numeric = all_results.select_dtypes(include=[np.number])
    assert np.isfinite(numeric.to_numpy()).all()
    assert "all_zero" in all_results.index
    row = all_results.loc["all_zero"]
    assert np.isfinite(float(row["logFC"]))
    assert np.isfinite(float(row["p_adj"]))
    assert np.isfinite(float(row["active_score"]))
    filt = scat.filter_active_genes(all_results, preset="permissive")
    assert isinstance(filt, pd.DataFrame)
    assert isinstance(significant, pd.DataFrame)


def test_two_cells_per_group_all_zero_velocity_layers():
    """Gene with zero spliced+unspliced in every cell should yield neutral velocity terms."""
    X = np.array([[1.0, 0.0, 3.0], [2.0, 0.0, 4.0], [3.0, 0.0, 5.0], [4.0, 0.0, 6.0]])
    adata = _tiny_velocity_adata(2, 2, X=X, var_names=["nz", "all_zero", "nz2"])
    adata.layers["spliced"][:, 1] = 0.0
    adata.layers["unspliced"][:, 1] = 0.0

    _, _, all_results = scat.active_score(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        use_permutation=False,
        min_cells=1,
        min_total_counts=0,
        show_plot=False,
    )

    row = all_results.loc["all_zero"]
    assert abs(float(row["unspliced_excess_residual"])) < 1e-6
    assert float(row["active_score"]) == 0.0
    assert np.isfinite(float(row["logFC"]))


def test_two_cells_per_group_entire_matrix_zero_counts():
    """Fully zero count matrix should complete without crash (DE may be degenerate)."""
    X = np.zeros((4, 3), dtype=float)
    adata = _tiny_velocity_adata(2, 2, X=X)

    _, _, all_results = scat.active_score(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        use_permutation=False,
        min_cells=1,
        min_total_counts=0,
        show_plot=False,
    )

    assert len(all_results) == adata.n_vars
    assert np.isfinite(all_results.select_dtypes(include=[np.number]).to_numpy()).all()


def test_diagnose_and_recommend_workflow_survive_one_cell_per_group():
    """Pre-flight helpers must not crash on degenerate cell counts."""
    adata = _tiny_velocity_adata(1, 1, X=np.array([[5.0], [1.0]]), sample_col="sample")

    diag = scat.diagnose_design(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        sample_col="sample",
        copy_input=False,
    )
    assert diag["n_cells_target"] == 1
    assert diag["n_cells_reference"] == 1
    assert isinstance(diag["warnings"], list)

    rec = scat.recommend_workflow(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        sample_col="sample",
    )
    assert rec["use_permutation"] is False
    assert "suggested_kwargs" in rec


def test_pseudobulk_one_cell_per_sample_scanpy_backend():
    """Pseudobulk with 1 cell/sample is tiny but must not crash on scanpy path."""
    X = np.random.poisson(5, size=(4, 6)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(
            {
                "condition": ["T", "T", "R", "R"],
                "sample": ["s1", "s2", "s3", "s4"],
            }
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(6)]),
    )
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.5

    _, _, all_results = scat.active_score(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        use_pseudobulk=True,
        sample_col="sample",
        pseudobulk_de_backend="scanpy",
        de_method="wilcoxon",
        min_cells=1,
        min_counts=0,
        use_permutation=False,
        show_plot=False,
    )
    for col in ("logFC", "p_val", "p_adj", "active_score", "unspliced_excess_residual"):
        assert col in all_results.columns
        assert np.isfinite(all_results[col].to_numpy()).all(), f"non-finite in {col}"


@pytest.mark.slow
def test_two_cells_per_group_permutation_completes():
    """Permutation on 2 cells/group is low-power but must not raise."""
    adata = _tiny_velocity_adata(2, 2)

    _, _, all_results = scat.active_score(
        adata,
        groupby="condition",
        target_group="T",
        reference_group="R",
        use_permutation=True,
        n_perm=10,
        n_jobs=1,
        min_cells=1,
        min_total_counts=0,
        show_plot=False,
    )
    assert "unspliced_excess_pval" in all_results.columns
    assert np.isfinite(all_results["unspliced_excess_pval"].to_numpy()).all()
