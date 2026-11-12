"""Tests for the review-fixes: rcParams leak, empty-input save, NaN coloring, returns."""

from __future__ import annotations

import warnings

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scatrans as scat
from scatrans.tl.filter import _sort_values_inplace


def _toy_adata_with_layers(n_obs=12, n_vars=4, with_nan_group=True):
    rng = np.random.default_rng(0)
    groups = np.array(["T"] * 4 + ["B"] * 4 + [None] * 4, dtype=object)
    if not with_nan_group:
        groups = np.array(["T"] * 6 + ["B"] * 6)
    adata = ad.AnnData(
        rng.poisson(3, size=(n_obs, n_vars)).astype(float),
        obs=pd.DataFrame({"celltype": groups}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_vars)]),
    )
    adata.layers["spliced"] = rng.random((n_obs, n_vars))
    adata.layers["unspliced"] = rng.random((n_obs, n_vars))
    adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"])
    return adata


def _mini_results(n=20):
    rng = np.random.default_rng(1)
    return pd.DataFrame(
        {
            "active_score": rng.normal(1.0, 0.5, n),
            "logFC": rng.normal(0.0, 1.0, n),
            "p_adj": rng.uniform(1e-6, 0.4, n),
            "unspliced_excess_residual": rng.normal(0.5, 1.0, n),
            "unspliced_excess_delta": rng.normal(0.2, 0.8, n),
            "gene_length": rng.integers(500, 4000, n),
            "intron_number": rng.integers(0, 15, n),
        },
        index=[f"g{i}" for i in range(n)],
    )


# ---------------------------------------------------------------------------
# Severe: use_style must not leak rcParams
# ---------------------------------------------------------------------------
@pytest.mark.plot
def test_use_style_does_not_leak_rcparams(tmp_path):
    key = "axes.linewidth"
    before = mpl.rcParams[key]
    results = _mini_results()

    fig, _ = scat.pl.bias_diagnostic_plot(results, show=False, use_style=True)
    plt.close(fig)
    assert mpl.rcParams[key] == before

    fig, _ = scat.pl.active_score_rankplot(results, top_n=5, show=False, use_style=True)
    plt.close(fig)
    assert mpl.rcParams[key] == before

    adata = _toy_adata_with_layers(with_nan_group=False)
    fig, _ = scat.pl.active_genes_heatmap(
        adata, genes=list(adata.var_names[:3]), groupby="celltype", show=False, use_style=True
    )
    if isinstance(fig, plt.Figure):
        plt.close(fig)
    assert mpl.rcParams[key] == before

    fig, _ = scat.pl.velocity_phase_portraits(
        adata, genes=list(adata.var_names[:2]), groupby="celltype", show=False, use_style=True
    )
    plt.close(fig)
    assert mpl.rcParams[key] == before


# ---------------------------------------------------------------------------
# Severe: empty-input paths still save (and would show)
# ---------------------------------------------------------------------------
@pytest.mark.plot
def test_gseaplot_empty_writes_save_path(tmp_path):
    dest = tmp_path / "gsea_empty.png"
    fig, ax = scat.pl.gseaplot(pd.Series(dtype=float), show=False, save_path=str(dest))
    plt.close(fig)
    assert dest.is_file() and dest.stat().st_size > 0
    assert fig is not None and ax is not None


@pytest.mark.plot
def test_active_genes_heatmap_empty_writes_save_path(tmp_path):
    dest = tmp_path / "heatmap_empty.png"
    adata = _toy_adata_with_layers(with_nan_group=False)
    fig, ax = scat.pl.active_genes_heatmap(adata, genes=None, show=False, save_path=str(dest))
    plt.close(fig)
    assert dest.is_file() and dest.stat().st_size > 0
    assert isinstance(fig, plt.Figure)


@pytest.mark.plot
def test_velocity_phase_portraits_empty_writes_save_path(tmp_path):
    dest = tmp_path / "phase_empty.png"
    adata = _toy_adata_with_layers(with_nan_group=False)
    fig, axes = scat.pl.velocity_phase_portraits(adata, genes=[], show=False, save_path=str(dest))
    plt.close(fig)
    assert dest.is_file() and dest.stat().st_size > 0
    assert isinstance(fig, plt.Figure)


# ---------------------------------------------------------------------------
# Severe: heatmap returns (Figure, axes-or-dict)
# ---------------------------------------------------------------------------
@pytest.mark.plot
def test_active_genes_heatmap_returns_figure(tmp_path):
    dest = tmp_path / "heatmap.png"
    adata = _toy_adata_with_layers(with_nan_group=False)
    fig, ax_or_dict = scat.pl.active_genes_heatmap(
        adata,
        genes=list(adata.var_names[:3]),
        groupby="celltype",
        show=False,
        save_path=str(dest),
    )
    assert isinstance(fig, plt.Figure)
    assert hasattr(fig, "savefig")
    assert dest.is_file() and dest.stat().st_size > 0
    plt.close(fig)


# ---------------------------------------------------------------------------
# Severe: NaN / category coloring in velocity_phase_portraits
# ---------------------------------------------------------------------------
@pytest.mark.plot
def test_velocity_phase_portraits_nan_groups_distinct_colors():
    adata = _toy_adata_with_layers(with_nan_group=True)
    # Distinct uns colors so we can also verify the helper is honored
    adata.uns["celltype_colors"] = ["#1f77b4", "#ff7f0e"]  # T, B

    fig, axes = scat.pl.velocity_phase_portraits(
        adata, genes=["g0"], groupby="celltype", show=False
    )
    fig.canvas.draw()
    ax = axes[0]
    collections = [c for c in ax.collections if hasattr(c, "get_facecolors")]
    plt.close(fig)
    assert collections, "expected scatter collections on the phase portrait"

    # Gather facecolors together with the observation mask they represent
    # by re-checking unique RGB tuples.
    colors = []
    for coll in collections:
        fc = np.asarray(coll.get_facecolors())
        if fc.size == 0:
            continue
        # unique rows
        rounded = np.round(fc[:, :3], 4)
        for row in np.unique(rounded, axis=0):
            colors.append(tuple(row))
    uniq = set(colors)
    # T, B, and NaN (gray) must not collapse to 1 or 2 colors
    assert len(uniq) >= 3, f"expected >=3 distinct colors (T, B, NA), got {uniq}"


@pytest.mark.plot
def test_colors_from_uns_prefers_stored_colors():
    adata = _toy_adata_with_layers(with_nan_group=False)
    adata.uns["celltype_colors"] = ["#000000", "#ffffff"]
    cols = scat.pl.colors_from_uns(adata, "celltype", categories=["T", "B"])
    assert cols[:2] == ["#000000", "#ffffff"]


# ---------------------------------------------------------------------------
# Medium helpers covered by unit-style tests
# ---------------------------------------------------------------------------
def test_palette_warns_when_cycling():
    with pytest.warns(UserWarning, match="will cycle"):
        cols = scat.pl.palette("okabe_ito", n=20)
    assert len(cols) == 20


def test_sort_values_inplace_preserves_identity_and_order():
    df = pd.DataFrame({"a": [3, 1, 2], "b": ["x", "y", "z"]})
    df.attrs["keep_me"] = 42
    original = df
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _sort_values_inplace(df, "a", ascending=True)
    assert original is df
    assert list(df["a"]) == [1, 2, 3]
    assert list(df["b"]) == ["y", "z", "x"]
    assert df.attrs["keep_me"] == 42


@pytest.mark.plot
def test_de_summary_barplot_accepts_n_up_n_down_aliases():
    summary = pd.DataFrame({"group": ["T", "B"], "n_up": [10, 3], "n_down": [2, 8]})
    fig, ax = scat.pl.de_summary_barplot(summary, show=False)
    assert fig is not None
    plt.close(fig)


@pytest.mark.plot
def test_de_summary_barplot_bad_columns_raise():
    summary = pd.DataFrame({"foo": [1], "bar": [2]})
    with pytest.raises(ValueError, match="n_up"):
        scat.pl.de_summary_barplot(summary, show=False)
