"""Default (non-@plot) smoke tests for scatrans.pl.

These run in the default pytest selection so pl coverage does not depend solely
on the optional ``plot`` marker job. All figures use the Agg backend.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import scatrans as scat


@pytest.fixture
def mini_results() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    n = 40
    idx = [f"G{i}" for i in range(n)]
    return pd.DataFrame(
        {
            "logFC": rng.normal(0.5, 1.0, n),
            "p_adj": rng.uniform(1e-4, 0.2, n),
            "active_score": rng.uniform(10, 90, n),
            "unspliced_excess_residual": rng.normal(0.5, 1.0, n),
            "unspliced_excess_delta": rng.normal(0.2, 0.8, n),
            "gene_length": rng.integers(500, 4000, n),
            "intron_number": rng.integers(0, 15, n),
        },
        index=idx,
    )


@pytest.fixture
def mini_enrich() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Term": ["pathway_a", "pathway_b", "pathway_c"],
            "Description": ["A", "B", "C"],
            "Count": [8, 5, 3],
            "GeneRatio": [0.2, 0.1, 0.05],
            "BgRatio": [0.05, 0.04, 0.03],
            "pvalue": [1e-4, 1e-3, 0.02],
            "p.adjust": [1e-3, 0.01, 0.05],
            "Genes": ["G0/G1/G2", "G3/G4", "G5"],
        }
    )


def test_smoke_volcano_and_comet(mini_results):
    fig, ax = scat.pl.volcano_plot(mini_results, top_n=5, show=False)
    assert fig is not None and ax is not None
    plt.close(fig)
    fig2, ax2 = scat.pl.comet_plot(mini_results, top_n=5, show=False)
    assert fig2 is not None
    plt.close(fig2)


def test_smoke_rankplot_bias_bar(mini_results, mini_enrich):
    fig, ax = scat.pl.active_score_rankplot(mini_results, top_n=8, show=False)
    plt.close(fig)
    fig2, ax2 = scat.pl.bias_diagnostic_plot(mini_results, show=False)
    plt.close(fig2)
    fig3, ax3 = scat.pl.enrich_barplot(mini_enrich, top_n=3, show=False)
    plt.close(fig3)
    fig4, ax4 = scat.pl.enrich_dotplot(mini_enrich, top_n=3, show=False)
    plt.close(fig4)


def test_enrich_barplot_labels_when_description_empty():
    """Regression: bundled GO/KEGG libraries ship an empty ``Description`` column
    (the readable name and ID both live in ``Term``). Empty strings pass
    ``notna()``, so enrich_barplot must not select ``Description`` and render blank
    bar labels — it should fall back to ``Term`` like enrich_dotplot does.
    """
    df = pd.DataFrame(
        {
            "Term": [
                "myeloid leukocyte activation (GO:0002274)",
                "lymphocyte proliferation (GO:0046651)",
                "leukocyte degranulation (GO:0043299)",
            ],
            "Description": ["", "", ""],  # matches real bundled ORA output
            "Count": [40, 37, 30],
            "GeneRatio": [0.057, 0.052, 0.043],
            "BgRatio": [0.01, 0.01, 0.008],
            "pvalue": [5e-17, 2e-13, 1e-10],
            "p.adjust": [5e-13, 4e-10, 1e-7],
            "Genes": ["G0/G1", "G2/G3", "G4"],
        }
    )
    fig, ax = scat.pl.enrich_barplot(df, top_n=3, show=False)
    labels = [t.get_text().strip() for t in ax.get_yticklabels()]
    plt.close(fig)
    assert labels, "enrich_barplot produced no y tick labels"
    assert all(labels), f"enrich_barplot rendered blank labels from empty Description: {labels}"


def test_smoke_volcano_3d_and_styles(mini_results):
    fig, ax = scat.pl.volcano_3d(mini_results, show=False)
    assert fig is not None
    plt.close(fig)
    # empty input must not crash (placeholder figure)
    empty = pd.DataFrame(columns=["logFC", "p_adj", "active_score"])
    fig_e, ax_e = scat.pl.volcano_3d(empty, show=False)
    assert fig_e is not None
    plt.close(fig_e)
    scat.pl.set_style()
    scat.pl.set_nature_style()
    with scat.pl.style_context(linewidth=0.5):
        fig2, _ = scat.pl.volcano_plot(mini_results, top_n=3, show=False, use_style=False)
        plt.close(fig2)


def test_smoke_phase_portraits(adata_basic):
    genes = list(adata_basic.var_names[:3])
    fig, axes = scat.pl.velocity_phase_portraits(
        adata_basic, genes, groupby="condition", show=False
    )
    plt.close(fig)


def test_smoke_gamma_shrinkage_if_columns(mini_results):
    df = mini_results.copy()
    # Optional columns some gamma plots expect — exercise graceful path
    df["effective_gamma"] = 0.3
    df["gamma_raw"] = 0.4
    try:
        out = scat.pl.gamma_shrinkage_plot(df, show=False)
        if out is not None:
            fig = out[0] if isinstance(out, tuple) else out
            plt.close(fig)
    except (TypeError, ValueError, KeyError):
        # Signature/column requirements may vary; smoke must not fail the suite
        # if optional diagnostic plot is picky — still covered by @plot tests.
        pytest.skip("gamma_shrinkage_plot requires additional columns/signature")


def test_pl_all_symbols_exist():
    for name in scat.pl.__all__:
        assert hasattr(scat.pl, name), name


def test_active_score_show_plot_imports_scatrans_pl(adata_basic):
    """Regression: after tl package split, show_plot must use ``from .. import pl``.

    Runs in the default suite (not only ``-m plot``) so wrong relative imports
    fail CI immediately.
    """
    plt.close("all")
    before = len(plt.get_fignums())
    scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=False,
        show_plot=True,
    )
    assert len(plt.get_fignums()) > before
    plt.close("all")


def test_smoke_volcano_styles_and_export(mini_results, tmp_path):
    fig, ax = scat.pl.volcano_plot(mini_results, style="ggvolcano", top_n=4, show=False, title="")
    assert fig is not None
    plt.close(fig)
    fig2, ax2 = scat.pl.volcano_plot(mini_results, style="gradual", top_n=3, show=False)
    plt.close(fig2)
    fig3, ax3 = scat.pl.comet_plot(mini_results, top_n=4, s=2, show=False)
    with scat.pl.figure_export_context(str(tmp_path / "figs"), fmt="png", dpi=80) as export:
        path = export.save(fig3, "comet")
        assert path.endswith(".png")
    paths = scat.pl.save_all_figures(
        {"volcano": fig2}, str(tmp_path / "batch"), fmt="png", dpi=80, close=True
    )
    assert len(paths) >= 1
    plt.close("all")


def test_smoke_compare_enrich_venn_upset(mini_enrich):
    df = mini_enrich.copy()
    df["Cluster"] = ["A", "A", "B"]
    fig, ax = scat.pl.enrich_vennplot(df, show=False)
    plt.close(fig)
    try:
        fig2, ax2 = scat.pl.enrich_upsetplot(df, show=False)
        plt.close(fig2)
    except ImportError:
        pytest.skip("upsetplot optional dependency not installed")


def _three_cluster_enrich() -> pd.DataFrame:
    """Three groups with known overlaps: up∩down∩shared = 3, etc."""
    def mock(cluster, terms):
        return pd.DataFrame(
            {"Term": terms, "Cluster": cluster, "p.adjust": [0.001] * len(terms),
             "Count": [10] * len(terms)}
        )
    up = [f"UP_{i}" for i in range(7)] + [f"S_{i}" for i in range(4)]
    down = [f"DN_{i}" for i in range(6)] + [f"S_{i}" for i in range(1, 5)]
    shared = [f"S_{i}" for i in range(6)] + [f"MID_{i}" for i in range(4)]
    return pd.concat([mock("up", up), mock("down", down), mock("shared", shared)],
                     ignore_index=True)


def test_enrich_vennplot_labels_all_three_set_regions():
    """Regression: a 3-set Venn must label all 7 exclusive regions, not just the
    three singletons (a 0.9.x bug only drew the pairwise-of-first-two counts).
    """
    fig, ax = scat.pl.enrich_vennplot(_three_cluster_enrich(), pval_cutoff=0.05, show=False)
    nums = [t.get_text() for t in ax.texts if t.get_text().strip().lstrip("-").isdigit()]
    plt.close(fig)
    assert len(nums) == 7, f"expected 7 region counts for a 3-set Venn, got {len(nums)}: {nums}"


def test_compare_dotplot_grid_layout():
    """compare_dotplot must lay groups on the x-axis and terms on the y-axis, with a
    dot at each enriched (group, term) cell (clusterProfiler compareCluster grid).
    """
    df = _three_cluster_enrich()  # up / down / shared, with shared S_* terms
    fig, ax = scat.pl.compare_dotplot(df, top_n=4, show=False)
    xlabels = [t.get_text().split("\n")[0] for t in ax.get_xticklabels()]
    n_dots = len(ax.collections[0].get_offsets()) if ax.collections else 0
    plt.close(fig)
    # x-axis carries the group names (order-independent)
    assert set(xlabels) == {"up", "down", "shared"}, xlabels
    # a shared term enriched in multiple groups yields multiple dots -> more dots than rows
    assert n_dots >= len(set(xlabels)), f"expected a grid of dots, got {n_dots}"


def test_enrich_upsetplot_matrix_has_no_duplicate_row_labels():
    """Regression: the UpSet matrix panel must not repeat the group names — they
    are shown once on the aligned Set-size panel; repeating them collided with the
    set-size value labels (e.g. "10" + "down" -> "10down").
    """
    fig, ax_mat = scat.pl.enrich_upsetplot(_three_cluster_enrich(), pval_cutoff=0.05, show=False)
    mat_labels = [t.get_text() for t in ax_mat.get_yticklabels() if t.get_text().strip()]
    plt.close(fig)
    assert mat_labels == [], f"matrix axis should carry no row labels, got {mat_labels}"


def test_smoke_heatmap_and_gamma(adata_basic, mini_results):
    genes = list(adata_basic.var_names[:3])
    out = scat.pl.active_genes_heatmap(adata_basic, genes=genes, groupby="condition", show=False)
    if out is not None:
        fig = out[0] if isinstance(out, tuple) else out
        if hasattr(fig, "savefig"):
            plt.close(fig)
    df = mini_results.copy()
    df["gamma_shrinkage_weight"] = np.linspace(0.1, 0.9, len(df))
    df["total_us_counts"] = np.linspace(10, 500, len(df))
    df["effective_gamma"] = np.linspace(0.2, 0.8, len(df))
    fig, ax = scat.pl.gamma_shrinkage_plot(df, show=False)
    assert fig is not None
    plt.close(fig)
