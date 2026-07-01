"""Tests for scatrans.pl plotting helpers (headless Agg backend)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import pytest

import scatrans as scat

pytestmark = pytest.mark.plot


def test_volcano_plot(results_df):
    fig, ax = scat.pl.volcano_plot(results_df, top_n=5, show=False)
    plt.close(fig)


def test_bias_diagnostic_plot(results_df):
    fig, ax = scat.pl.bias_diagnostic_plot(results_df, show=False)
    plt.close(fig)


def test_enrich_barplot(enrich_df):
    fig, ax = scat.pl.enrich_barplot(enrich_df, top_n=2, show=False)
    plt.close(fig)


def test_velocity_phase_portraits(adata_basic):
    genes = adata_basic.var_names[:3].tolist()
    fig, axes = scat.pl.velocity_phase_portraits(
        adata_basic, genes, groupby="condition", show=False
    )
    plt.close(fig)


def test_style_context():
    with scat.pl.style_context(linewidth=0.5):
        pass
    scat.pl.set_nature_style()


def test_enrich_upsetplot(enrich_df):
    try:
        fig, ax = scat.pl.enrich_upsetplot(enrich_df, show=False)
        plt.close(fig)
    except ImportError:
        pytest.skip("upsetplot not installed")


def test_enrich_vennplot(enrich_df):
    fig, ax = scat.pl.enrich_vennplot(enrich_df, show=False)
    plt.close(fig)


def test_volcano_3d_empty_df():
    empty = pd.DataFrame(columns=["logFC", "p_adj", "active_score"])
    fig, ax = scat.pl.volcano_3d(empty, show=False)
    plt.close(fig)
