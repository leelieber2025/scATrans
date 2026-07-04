"""Additional plotting tests."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

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


def test_gseaplot_dataframe_first_arg_type_error():
    gsea_res = pd.DataFrame({"Term": ["TERM1"], "NES": [1.5], "pvalue": [0.01]})
    with pytest.raises(TypeError, match="gsea_result"):
        scat.pl.gseaplot(gsea_res, show=False)


def test_gseaplot_with_stored_curves():
    # gseapy-sorted ranking (descending scores)
    ranking = pd.Series(
        {f"G{i}": float(20 - i) for i in range(20)},
        index=[f"G{i}" for i in range(20)],
    ).sort_values(ascending=False)
    # Caller passes a deliberately shuffled series (common when ranked from all_results)
    shuffled = ranking.sample(frac=1, random_state=0)
    res_curve = np.linspace(0.0, 1.0, len(ranking))
    gsea_res = pd.DataFrame(
        {
            "Term": ["TERM1"],
            "NES": [1.5],
            "pvalue": [0.01],
            "p.adjust": [0.05],
        }
    )
    gsea_res.attrs["gsea_details"] = {
        "TERM1": {
            "RES": res_curve.tolist(),
            "hits": [2, 5, 9],
            "nes": 1.5,
            "pval": 0.01,
            "fdr": 0.05,
        }
    }
    gsea_res.attrs["ranking"] = ranking.to_dict()
    fig, ax = scat.pl.gseaplot(shuffled, gsea_res, term="TERM1", show=False)
    bar_heights = [p.get_height() for p in fig.axes[2].patches]
    assert np.allclose(bar_heights, ranking.values)
    assert len(fig.axes[0].lines[0].get_ydata()) == len(ranking)
    plt.close(fig)


@pytest.mark.plot
def test_active_genes_heatmap_default_show():
    import anndata as ad

    X = np.random.poisson(3, size=(20, 5)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"group": ["A"] * 10 + ["B"] * 10}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(5)]),
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    result = scat.pl.active_genes_heatmap(
        adata, genes=adata.var_names[:3].tolist(), groupby="group", show=True
    )
    assert result is not None
    fig, _ax = result
    assert fig is not None
    if isinstance(fig, plt.Figure):
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


def test_bias_diagnostic_plot_rejects_invalid_axes(results_df):
    with pytest.raises(ValueError, match="exactly two matplotlib Axes"):
        scat.pl.bias_diagnostic_plot(results_df, axes=(None, None), show=False)


def test_comet_plot_missing_columns_placeholder():
    df = pd.DataFrame({"logFC": [1.0], "active_score": [50.0]})
    fig, ax = scat.pl.comet_plot(df, show=False)
    assert fig is not None and ax is not None
    plt.close(fig)


def test_volcano_plot_missing_columns_placeholder():
    df = pd.DataFrame({"logFC": [1.0]})
    fig, ax = scat.pl.volcano_plot(df, show=False)
    assert fig is not None and ax is not None
    plt.close(fig)


def test_enrich_dotplot_auto_without_padj_column():
    """show_terms='auto' must not crash when p.adjust column is absent."""
    df = pd.DataFrame(
        {
            "Term": ["T1", "T2", "T3"],
            "Count": [5, 3, 2],
            "GeneRatio": [0.1, 0.2, 0.05],
        }
    )
    fig, ax = scat.pl.enrich_dotplot(df, show_terms="auto", top_n=2, show=False)
    plt.close(fig)
