"""Regression tests for confirmed bug fixes (2026-07 audit)."""

from __future__ import annotations

import logging
import sys
from unittest.mock import patch

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scatrans as scat
from scatrans._de import _run_de_wrapper, _validate_de_result
from scatrans._utils import _is_integer_counts_like
from scatrans.tl import _materialize_if_view, _select_var


def test_validate_de_result_before_fillna_in_run_de_wrapper(adata_de_only):
    """fillna must not mask all-NaN backend output (_run_de_wrapper end-to-end)."""
    adata = adata_de_only.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_genes = adata.n_vars
    nan_df = pd.DataFrame(
        {
            "names": adata.var_names.to_list(),
            "logfoldchanges": [np.nan] * n_genes,
            "pvals": [np.nan] * n_genes,
            "pvals_adj": [np.nan] * n_genes,
        }
    )
    with (
        patch("scatrans._de.sc.tl.rank_genes_groups"),
        patch("scatrans._de.sc.get.rank_genes_groups_df", return_value=nan_df),
        pytest.raises(RuntimeError, match="no finite values"),
    ):
        _run_de_wrapper(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            de_method="wilcoxon",
        )


def test_validate_de_result_raises_before_fillna_helper():
    """Direct helper still rejects all-NaN columns (unit-level guard)."""
    bad = pd.DataFrame(
        {"logFC": [np.nan], "p_val": [np.nan], "p_adj": [np.nan]},
        index=["G1"],
    )
    with pytest.raises(RuntimeError, match="no finite values"):
        _validate_de_result(bad, backend="test")


def test_gene_type_filter_copy_input_false_materializes_var_view():
    """_select_var + copy_input=False must be materialized before .var writes."""
    rng = np.random.default_rng(0)
    X = rng.poisson(3, size=(20, 8)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"condition": ["A"] * 10 + ["B"] * 10}),
        var=pd.DataFrame(
            {
                "gene_type": ["protein_coding"] * 5 + ["lncRNA"] * 3,
            },
            index=[f"G{i}" for i in range(8)],
        ),
    )
    mask = adata.var["gene_type"] == "protein_coding"
    sub = _select_var(adata, mask, copy_input=False)
    sub = _materialize_if_view(sub)
    assert not getattr(sub, "is_view", False)
    sub.var["marker"] = 1
    assert "marker" in sub.var.columns


def test_core_modules_have_null_handler():
    """tl/_de/_velocity/_permutation/_utils loggers must attach NullHandler."""
    modules = [
        "scatrans.tl",
        "scatrans._de",
        "scatrans._velocity",
        "scatrans._permutation",
        "scatrans._utils",
    ]
    for name in modules:
        logger = logging.getLogger(name)
        assert any(isinstance(h, logging.NullHandler) for h in logger.handlers), name


def test_run_kegg_accepts_hs_and_mm_shorthand():
    """organism='hs'/'mm' must work in run_kegg like other enrich functions."""
    genes = ["Tp53", "Gapdh", "Actb", "Mdm2", "Cdkn1a"]
    scat.run_kegg(genes, organism="hs", verbose=False)
    scat.run_kegg(genes, organism="mm", verbose=False)


def test_de_auto_after_scale_yields_nonzero_logfc(adata_de_only):
    """normalize→log1p→scale standard scanpy path must not zero-out logFC."""
    adata = adata_de_only.copy()
    adata.X[40:, :15] *= 2.5
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    _, results = scat.differential_expression(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        de_preprocess="auto",
    )
    assert (results["logFC"].abs() > 1e-6).any()


def test_add_gene_features_organism_case_insensitive():
    """organism='Human' must resolve like 'human', not silently fall back to mouse."""
    adata = ad.AnnData(
        np.ones((4, 4)),
        var=pd.DataFrame(index=["TP53", "GAPDH", "ACTB", "MALAT1"]),
    )
    out_lower = scat.add_gene_features(adata.copy(), organism="human")
    out_title = scat.add_gene_features(adata.copy(), organism="Human")
    assert out_lower.var["gene_length"].notna().all()
    assert out_title.var["gene_length"].notna().all()
    pd.testing.assert_series_equal(
        out_lower.var["gene_length"],
        out_title.var["gene_length"],
        check_names=False,
    )


def test_is_integer_counts_like_detects_stride_contamination():
    """Stride subsample must hit contamination at index 0 on large matrices."""
    rng = np.random.default_rng(0)
    mat = rng.integers(0, 12, size=500_000).astype(float)
    mat[0] = 0.5
    assert not _is_integer_counts_like(mat)


def test_diagnose_design_kb_python_layers(adata_mature_nascent):
    """mature/nascent (kb_python) must populate unspliced_global_fraction."""
    diag = scat.diagnose_design(
        adata_mature_nascent,
        groupby="condition",
        target_group="GA",
        reference_group="Ctrl",
    )
    assert diag["unspliced_global_fraction"] is not None
    assert 0.0 <= diag["unspliced_global_fraction"] <= 1.0


@pytest.mark.plot
def test_comet_plot_style_context_restored_on_validation_error():
    """Simulate IPython traceback retention: rcParams must restore after ValueError."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = pd.DataFrame(
        {
            "logFC": [1.0],
            "active_score": [1.0],
            "unspliced_excess_residual": [0.2],
        },
        index=["G1"],
    )
    before = float(mpl.rcParams["figure.dpi"])
    try:
        scat.pl.comet_plot(df, use_style=True, top_n=-1, show=False)
    except ValueError:
        sys.last_traceback = sys.exc_info()[2]
    else:
        pytest.fail("expected ValueError for top_n=-1")
    assert float(mpl.rcParams["figure.dpi"]) == before
    plt.rcdefaults()
