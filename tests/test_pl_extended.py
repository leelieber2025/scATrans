"""Additional plotting tests."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scatrans as scat

pytestmark = pytest.mark.plot


def test_comet_plot_point_scale_and_s(results_df):
    fig, ax = scat.pl.comet_plot(results_df, top_n=4, s=2, show=False)
    plt.close(fig)
    fig2, ax2 = scat.pl.comet_plot(results_df, top_n=4, point_scale=0.2, min_size=1, show=False)
    plt.close(fig2)


def test_volcano_plot_label_genes(results_df):
    genes = results_df.index[:2].tolist()
    fig, ax = scat.pl.volcano_plot(results_df, label_genes=genes, top_n=3, show=False)
    plt.close(fig)


def test_gseaplot_with_stored_curves():
    ranked = pd.Series(np.linspace(1, 0, 20), index=[f"G{i}" for i in range(20)])
    gsea_res = pd.DataFrame(
        {
            "Term": ["TERM1"],
            "NES": [1.5],
            "pvalue": [0.01],
            "p.adjust": [0.05],
        }
    )
    gsea_res.attrs["gsea_info"] = {
        "rank_metric": ranked,
        "enrichment_score": np.linspace(0, 0.5, 20),
        "hit_indices": [2, 5, 9],
    }
    fig, ax = scat.pl.gseaplot(ranked, gsea_res, term="TERM1", show=False)
    plt.close(fig)


def test_active_score_rankplot_empty():
    empty = pd.DataFrame(columns=["active_score", "logFC", "p_adj"])
    fig, ax = scat.pl.active_score_rankplot(empty, show=False)
    plt.close(fig)


def test_enrich_dotplot_cluster_col(enrich_df):
    df = enrich_df.copy()
    df["cluster"] = ["A", "A", "B"]
    fig, ax = scat.pl.enrich_dotplot(df, top_n=3, show=False)
    plt.close(fig)
