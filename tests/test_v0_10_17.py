"""Regression tests for the 0.10.17 bugfix release.

Default pytest suite (no plot/slow markers). Plotters use the Agg backend.
"""

from __future__ import annotations

import inspect

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scatrans as scat
from scatrans._utils import _categorical_or_unique
from scatrans.tl._common import (
    COMPARE_LOGFC_CUTOFF,
    HEURISTIC_FILTER_DEFAULTS,
    PARTITION_LOGFC_CUTOFF,
)
from scatrans.tl.filter import filter_active_genes
from scatrans.tl.pipeline import run_default_pipeline


def test_version_is_0_10_17():
    assert scat.__version__ == "0.10.17"


def test_logfc_cutoffs_remain_distinct():
    """Partition / heuristic / compare cutoffs differ on purpose."""
    assert PARTITION_LOGFC_CUTOFF == 1.0
    assert HEURISTIC_FILTER_DEFAULTS["logfc_cutoff"] == 0.35
    assert COMPARE_LOGFC_CUTOFF == 0.25
    assert (
        len(
            {
                PARTITION_LOGFC_CUTOFF,
                HEURISTIC_FILTER_DEFAULTS["logfc_cutoff"],
                COMPARE_LOGFC_CUTOFF,
            }
        )
        == 3
    )
    sig = inspect.signature(scat.partition_de_by_mechanism)
    assert sig.parameters["logfc_cutoff"].default == PARTITION_LOGFC_CUTOFF
    sig_c = inspect.signature(scat.compare_de_across_groups)
    assert sig_c.parameters["logfc_cutoff"].default == COMPARE_LOGFC_CUTOFF


# ---------------------------------------------------------------------------
# A. object-dtype grouping must not crash; categoricals keep unused-aware order
# ---------------------------------------------------------------------------
def _mechanism_table(*, categorical: bool, extra_unused: bool = False) -> pd.DataFrame:
    classes = ["transcription-driven", "stabilization-driven", "ambiguous"] * 4
    df = pd.DataFrame(
        {
            "active_score": np.linspace(10, 90, len(classes)),
            "mechanism_class": classes,
            "cell_type": (["T"] * 6) + (["B"] * 6),
        }
    )
    if categorical:
        cats = [
            "transcription-driven",
            "stabilization-driven",
            "ambiguous",
            "unclassified_down",
        ]
        if extra_unused:
            df["mechanism_class"] = pd.Categorical(df["mechanism_class"], categories=cats)
            df["cell_type"] = pd.Categorical(df["cell_type"], categories=["T", "B", "NK"])
        else:
            df["mechanism_class"] = pd.Categorical(df["mechanism_class"], categories=cats[:3])
    else:
        df["mechanism_class"] = df["mechanism_class"].astype(object)
        df["cell_type"] = df["cell_type"].astype(object)
    return df


def test_group_stat_plot_object_dtype_groupby():
    df = _mechanism_table(categorical=False)
    fig, ax = scat.pl.group_stat_plot(
        df, value_col="active_score", groupby="mechanism_class", show=False
    )
    assert fig is not None and ax is not None
    labels = [t.get_text() for t in ax.get_xticklabels()]
    assert "transcription-driven" in labels
    plt.close(fig)


def test_group_stat_plot_categorical_unused_order():
    df = _mechanism_table(categorical=True, extra_unused=True)
    fig, ax = scat.pl.group_stat_plot(
        df, value_col="active_score", groupby="mechanism_class", show=False
    )
    labels = [t.get_text() for t in ax.get_xticklabels()]
    # unused "unclassified_down" is omitted; remaining follow categorical order
    assert labels == ["transcription-driven", "stabilization-driven", "ambiguous"]
    plt.close(fig)


def test_composition_barplot_object_dtype_groupby_and_hue():
    df = _mechanism_table(categorical=False)
    fig, ax, table = scat.pl.composition_barplot(
        df, groupby="cell_type", hue="mechanism_class", show=False, return_table=True
    )
    assert fig is not None
    assert "T" in list(table.index) and "B" in list(table.index)
    plt.close(fig)


def test_composition_barplot_categorical_unused_order():
    df = _mechanism_table(categorical=True, extra_unused=True)
    fig, ax, table = scat.pl.composition_barplot(
        df, groupby="cell_type", hue="mechanism_class", show=False, return_table=True
    )
    assert list(table.index) == ["T", "B"]  # unused NK dropped
    assert list(table.columns) == [
        "transcription-driven",
        "stabilization-driven",
        "ambiguous",
    ]
    plt.close(fig)


def test_categorical_or_unique_object_vs_categorical():
    obj = pd.Series(["b", "a", "b"], dtype=object)
    assert _categorical_or_unique(obj) == ["b", "a"]
    cat = pd.Categorical(["b", "a", "b"], categories=["a", "b", "c"])
    assert _categorical_or_unique(pd.Series(cat)) == ["a", "b"]
    assert _categorical_or_unique(pd.Series(cat), drop_unused=False) == ["a", "b", "c"]


# ---------------------------------------------------------------------------
# B. compare_de_across_groups object-dtype split_by
# ---------------------------------------------------------------------------
def test_compare_de_across_groups_object_dtype_split_by(adata_basic):
    adata = adata_basic.copy()
    adata.obs["cell_type"] = np.array(["T"] * 60 + ["B"] * 60, dtype=object)
    # Huge min_cells skips DE but still exercises group listing (the old hasattr line).
    result = scat.compare_de_across_groups(
        adata,
        split_by="cell_type",
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        min_cells=10**9,
    )
    assert set(result.failed) == {"T", "B"}


# ---------------------------------------------------------------------------
# C. assign_pseudo_replicates balanced assignment
# ---------------------------------------------------------------------------
def _two_condition_adata(n_per: int):
    import anndata as ad

    n = n_per * 2
    adata = ad.AnnData(np.ones((n, 3)))
    adata.obs["condition"] = ["A"] * n_per + ["B"] * n_per
    return adata


def test_assign_pseudo_replicates_balanced_and_reproducible():
    adata = _two_condition_adata(8)
    scat.assign_pseudo_replicates(adata, groupby="condition", n_replicates=3, random_state=0)
    col = adata.obs["pseudo_replicate"]
    assert isinstance(col.dtype, pd.CategoricalDtype)
    for cond in ("A", "B"):
        labs = col.loc[adata.obs["condition"] == cond]
        assert labs.nunique() == 3
        assert set(labs) == {f"{cond}_rep1", f"{cond}_rep2", f"{cond}_rep3"}

    a1 = _two_condition_adata(8)
    a2 = _two_condition_adata(8)
    scat.assign_pseudo_replicates(a1, groupby="condition", n_replicates=3, random_state=7)
    scat.assign_pseudo_replicates(a2, groupby="condition", n_replicates=3, random_state=7)
    assert list(a1.obs["pseudo_replicate"]) == list(a2.obs["pseudo_replicate"])


def test_assign_pseudo_replicates_too_few_cells_raises():
    adata = _two_condition_adata(2)
    with pytest.raises(ValueError, match="n_replicates"):
        scat.assign_pseudo_replicates(adata, groupby="condition", n_replicates=3)


# ---------------------------------------------------------------------------
# D. select_by default is "de"
# ---------------------------------------------------------------------------
def test_filter_and_pipeline_default_select_by_is_de():
    assert inspect.signature(filter_active_genes).parameters["select_by"].default == "de"
    assert inspect.signature(run_default_pipeline).parameters["select_by"].default == "de"

    df = pd.DataFrame(
        {
            "logFC": [3.0, 0.1, -2.0],
            "p_adj": [1e-6, 0.8, 1e-5],
            "active_score": [1.0, 90.0, 80.0],
            "unspliced_excess_residual": [-2.0, 5.0, 4.0],
        },
        index=["de_up_bad_proxy", "proxy_only", "de_down"],
    )
    default = filter_active_genes(df)
    explicit_de = filter_active_genes(df, select_by="de")
    assert list(default.index) == list(explicit_de.index)
    assert "de_up_bad_proxy" in default.index
    assert "proxy_only" not in default.index
    comp = filter_active_genes(df, preset="heuristic", select_by="composite")
    assert "de_up_bad_proxy" not in comp.index


def test_run_default_pipeline_default_select_by_de(adata_basic):
    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
    )
    assert res.meta["select_by"] == "de"


# ---------------------------------------------------------------------------
# E. permutation: one-NaN score vector still counts toward n_success
# ---------------------------------------------------------------------------
def test_permutation_single_nan_score_counts_toward_n_success(monkeypatch):
    import anndata as ad

    import scatrans._permutation as perm_mod
    from scatrans._permutation import run_permutation_test

    n_cells, n_genes = 20, 5
    rng = np.random.default_rng(0)
    X = rng.poisson(3, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"group": ["T"] * 10 + ["R"] * 10}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )
    adata.var["gene_length"] = np.linspace(800, 4000, n_genes)
    adata.var["intron_number"] = np.arange(n_genes)
    uns = X.copy()
    spl = X * 0.5
    real_score = np.full(n_genes, 50.0)
    real_residual = np.full(n_genes, 1.0)
    n_perm = 5

    def one_nan_task(*args, **kwargs):
        s = np.zeros(n_genes)
        s[0] = np.nan
        return s, np.zeros(n_genes)

    monkeypatch.setattr(perm_mod, "_single_permutation_task", one_nan_task)

    pvals, *_rest, use_fdr, reason = run_permutation_test(
        n_perm=n_perm,
        effective_n_jobs=1,
        random_seed=0,
        obs_labels=adata.obs["group"].to_numpy(),
        target_group="T",
        reference_group="R",
        adata=adata,
        X_features=None,
        valid_feat=np.ones(n_genes, dtype=bool),
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
        valid_expr=np.ones(n_genes, dtype=bool),
    )
    assert reason != "permutation_shuffle_failed"
    # scores are 0 (or NaN) vs observed 50 → no exceedances; p = 1/(n_success+1)
    expected = 1.0 / (n_perm + 1.0)
    assert np.allclose(pvals, expected)
    assert use_fdr is False


# ---------------------------------------------------------------------------
# F. MixedLM warning no longer claims logFC=0 fill
# ---------------------------------------------------------------------------
def test_mixedlm_warning_does_not_claim_logfc_zero_fill():
    from scatrans._de import _run_mixedlm_de

    src = inspect.getsource(_run_mixedlm_de)
    assert "mixedlm_coef=0" in src
    assert "logFC is still" in src or "from sample means" in src
    assert "received neutral values (logFC=0" not in src


# ---------------------------------------------------------------------------
# H. de_summary_barplot smoke
# ---------------------------------------------------------------------------
def test_de_summary_barplot_smoke():
    summary = pd.DataFrame({"up": [5, 3, 1], "down": [2, 4, 0]}, index=["T", "B", "NK"])
    fig, ax = scat.pl.de_summary_barplot(summary, mode="stacked", show=False)
    assert fig is not None and ax is not None
    plt.close(fig)
    fig2, ax2 = scat.pl.de_summary_barplot(summary, mode="diverging", show=False)
    plt.close(fig2)


def test_include_groups_compat_exon_union():
    """pp_bias exon union retries without include_groups on older pandas."""
    import scatrans.pp_bias as pb

    src = inspect.getsource(pb.generate_gene_features_from_gtf)
    assert "include_groups=False" in src
    assert "except TypeError" in src
