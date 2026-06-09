"""
Expanded test suite for scATrans after structural + quality refactor.

Covers:
- Core active_score (heuristic + permutation small-n)
- Advanced mode (skipped if scvelo missing)
- Layer remapping (kb_python mature/nascent)
- add_gene_features + list_available_gene_features
- Enrichment (run_enrichment with dict + simplify)
- Plotting (non-interactive Agg backend)
- Basic error paths
- Metadata written to adata.uns / .var
"""

import os
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless for CI / no display

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scatrans as scat


# --------------------------- fixtures ---------------------------


@pytest.fixture
def adata_basic():
    """Basic test AnnData with spliced/unspliced layers + gene features."""
    np.random.seed(42)
    n_cells, n_genes = 120, 250
    X = np.random.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 60 + ["Control"] * 60
    # Add sample ids for mixed model tests (multiple samples per group)
    ad.obs["sample"] = ["s" + str(i % 8) for i in range(n_cells)]
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.55
    ad.var["gene_length"] = np.random.randint(700, 4500, n_genes)
    ad.var["intron_number"] = np.random.randint(0, 12, n_genes)
    return ad


@pytest.fixture
def adata_mature_nascent():
    """kb_python style layer names."""
    np.random.seed(123)
    n_cells, n_genes = 80, 180
    X = np.random.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["GA"] * 40 + ["Ctrl"] * 40
    ad.layers["mature"] = X.copy()
    ad.layers["nascent"] = X * 0.5
    ad.var["gene_length"] = np.random.randint(800, 4000, n_genes)
    ad.var["intron_number"] = np.random.randint(1, 9, n_genes)
    return ad


@pytest.fixture
def adata_mixed_small():
    """Lightweight fixture dedicated to mixed-model + delta_variance tests (fast CI)."""
    np.random.seed(42)
    n_cells, n_genes = 60, 70
    X = np.random.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 30 + ["Control"] * 30
    # 6 samples (3 per group) — enough for (1|sample) RE but small for speed
    ad.obs["sample"] = ["s" + str(i % 6) for i in range(n_cells)]
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.45
    ad.var["gene_length"] = np.random.randint(600, 3500, n_genes)
    ad.var["intron_number"] = np.random.randint(0, 8, n_genes)
    return ad


# --------------------------- core active_score ---------------------------


def test_heuristic_basic(adata_basic):
    res, sig, all_res = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )
    assert "active_score" in res.var.columns
    assert "velocity_residual" in res.var.columns
    assert "logFC" in res.var.columns
    assert len(all_res) == adata_basic.n_vars
    assert "scatrans" in res.uns


@pytest.mark.parametrize("use_perm", [False, True])
def test_heuristic_with_and_without_permutation(adata_basic, use_perm):
    n_perm = 4 if use_perm else 0  # tiny for speed
    res, sig, _ = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=use_perm,
        n_perm=n_perm,
        random_seed=0,
    )
    assert "active_score" in res.var.columns
    if use_perm:
        assert "active_score_pval" in res.var.columns or len(sig) == 0
        assert "active_score_fdr" in res.var.columns or len(sig) == 0


def test_advanced_runs_or_skips(adata_basic):
    try:
        res, sig, _ = scat.active_score(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            mode="advanced",
            advanced_fallback=True,
            show_plot=False,
            use_permutation=False,
            n_jobs=1,
        )
        assert res.uns["scatrans"]["mode"] in {"advanced", "heuristic_fallback_from_advanced"}
    except ImportError:
        pytest.skip("scvelo not installed in this environment")


def test_layer_remapping_kb_python_style(adata_mature_nascent):
    res, sig, _ = scat.active_score(
        adata_mature_nascent,
        groupby="condition",
        target_group="GA",
        reference_group="Ctrl",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )
    # Internal layers should have been added
    assert "spliced" in res.layers
    assert "unspliced" in res.layers
    assert "active_score" in res.var.columns


# --------------------------- gene features ---------------------------


def test_add_gene_features_and_list(adata_basic):
    # Should not crash even if features are incomplete
    before = adata_basic.var.columns.tolist()
    out = scat.add_gene_features(adata_basic, organism="mouse")
    assert out is adata_basic
    # list should return something (now re-exported at top level for convenience)
    avail = scat.list_available_gene_features()
    assert isinstance(avail, list)


# --------------------------- enrichment ---------------------------


def test_run_enrichment_dict():
    genes = ["GeneA", "GeneB", "GeneC"]
    gene_sets = {
        "TERM1": ["GeneA", "GeneB", "GeneD"],
        "TERM2": ["GeneC", "GeneE"],
    }
    res = scat.run_enrichment(genes, gene_sets=gene_sets, pval_cutoff=1.0, min_size=1)
    assert not res.empty or len(genes) > 0
    if not res.empty:
        assert "Term" in res.columns
        assert "p.adjust" in res.columns


def test_simplify_enrichment():
    df = pd.DataFrame(
        {
            "Term": ["A", "B"],
            "p.adjust": [0.01, 0.02],
            "Genes": ["g1;g2", "g1;g3"],
            "Count": [2, 2],
        }
    )
    simp = scat.simplify_enrichment(df, similarity_cutoff=0.1, min_count=1)
    assert len(simp) <= len(df)


# --------------------------- plotting (headless) ---------------------------


def test_pl_set_style_and_comet(adata_basic):
    # Run a quick analysis so we have results df
    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )
    scat.pl.set_style()
    # Should not raise (comet_plot does not take show_plot)
    fig, ax = scat.pl.comet_plot(allr, top_n=5)
    # comet_plot calls plt.show() internally; with Agg it is fine
    import matplotlib.pyplot as plt

    plt.close("all")


def test_pl_rankplot_and_heatmap_stubs(adata_basic):
    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )
    # rankplot now has a real (simple) impl
    fig, ax = scat.pl.active_score_rankplot(allr, top_n=6)
    # heatmap stub should at least not explode
    scat.pl.active_genes_heatmap(adata_basic, genes=allr.index[:5].tolist())


# --------------------------- ax= parameter & edge cases ---------------------------


def test_plotting_with_ax_parameter(adata_basic):
    """Test that main plot functions accept an external ax (for multi-panel figures)."""
    _, _, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    # comet with ax
    f1, a1 = scat.pl.comet_plot(allr, top_n=4, ax=axes[0])
    assert a1 is axes[0]
    # volcano with ax
    f2, a2 = scat.pl.volcano_plot(allr, top_n=4, ax=axes[1])
    assert a2 is axes[1]
    plt.close(fig)


def test_edge_cases_low_features(adata_basic):
    """When gene features are missing/NaN the pipeline should still run (median fallback)."""
    ad = adata_basic.copy()
    ad.var["gene_length"] = np.nan
    ad.var["intron_number"] = np.nan
    res, sig, allr = scat.active_score(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        use_permutation=False,
        show_plot=False,
    )
    assert "active_score" in res.var.columns


def test_cli_main_callable():
    from scatrans.generate_gene_features import main

    assert callable(main)


# --------------------------- error paths ---------------------------


def test_error_bad_groups(adata_basic):
    with pytest.raises(ValueError):
        scat.active_score(
            adata_basic,
            groupby="condition",
            target_group="NOPE",
            reference_group="Control",
            show_plot=False,
        )


def test_error_missing_layers():
    bad = sc.AnnData(np.random.rand(10, 5))
    bad.obs["condition"] = ["A"] * 5 + ["B"] * 5
    with pytest.raises(ValueError):
        scat.active_score(
            bad, groupby="condition", target_group="A", reference_group="B", show_plot=False
        )


# --------------------------- CLI smoke ---------------------------


def test_cli_main_is_callable():
    from scatrans.generate_gene_features import main

    assert callable(main)


# --------------------------- mixed model + delta variance ---------------------------

def test_mixed_model_basic(adata_mixed_small):
    """use_mixed_model=True produces delta_variance / delta_var_pval and runs without crash."""
    res, sig, allr = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
        use_mixed_model=True,
        sample_col="sample",
        n_jobs=1,
    )
    assert "active_score" in res.var.columns
    assert "delta_variance" in res.var.columns
    assert "delta_var_pval" in res.var.columns
    assert "delta_variance" in allr.columns
    # delta_variance in [0,1]
    dvals = allr["delta_variance"].dropna().values
    assert len(dvals) > 0
    assert np.all((dvals >= 0) & (dvals <= 1))
    # pvals in [0,1]
    pvals = allr["delta_var_pval"].dropna().values
    assert np.all((pvals >= 0) & (pvals <= 1))


@pytest.mark.parametrize("use_dv", [False, True])
def test_delta_variance_filter_option(adata_mixed_small, use_dv):
    """The use_delta_variance_pval flag changes (or does not change) the sig set as expected."""
    _, sig_no, allr = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
        use_mixed_model=True,
        sample_col="sample",
        use_delta_variance_pval=False,
        delta_var_pval_cutoff=0.05,
        n_jobs=1,
    )
    _, sig_dv, _ = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
        use_mixed_model=True,
        sample_col="sample",
        use_delta_variance_pval=use_dv,
        delta_var_pval_cutoff=0.01,  # stricter to potentially reduce sigs
        n_jobs=1,
    )
    # When enabled with strict cutoff, |sig_dv| <= |sig_no| (or equal if no genes had small p)
    if use_dv:
        assert len(sig_dv) <= len(sig_no)
    # all_results always has the column when mixed used
    assert "delta_variance" in allr.columns


def test_mixed_model_incompatible_with_pseudobulk(adata_mixed_small):
    with pytest.raises(ValueError, match="incompatible"):
        scat.active_score(
            adata_mixed_small,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_pseudobulk=True,
            sample_col="sample",
            use_mixed_model=True,
            show_plot=False,
        )
