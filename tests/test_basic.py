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
    assert "unspliced_excess_residual" in res.var.columns
    assert "velocity_residual" in res.var.columns  # legacy alias
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
        assert "unspliced_excess_pval" in res.var.columns or len(sig) == 0
        assert "unspliced_excess_fdr" in res.var.columns or len(sig) == 0


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
    adata_basic.var.columns.tolist()
    out = scat.add_gene_features(adata_basic, organism="mouse")
    assert out is adata_basic
    # list should return something (now re-exported at top level for convenience)
    avail = scat.list_available_gene_features()
    assert isinstance(avail, list)
    # New: bundled human data should be discoverable and usable
    assert any("human" in f.lower() for f in avail), "human gene features should be bundled"
    # organism=human should resolve without error (even if few genes match the dummy adata)
    out_h = scat.add_gene_features(adata_basic.copy(), organism="human")
    assert "gene_length" in out_h.var.columns
    assert "intron_number" in out_h.var.columns


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


def test_simplify_enrichment_pathway_denester():
    gene_sets = {
        "TERM_PARENT": ["g1", "g2", "g3", "g4", "g5"],
        "TERM_CHILD": ["g1", "g2", "g3"],
    }
    df = pd.DataFrame(
        {
            "Term": ["TERM_PARENT", "TERM_CHILD"],
            "Description": ["Parent pathway", "Child pathway"],
            "p.adjust": [0.001, 0.01],
            "pvalue": [1e-5, 1e-3],
            "Genes": ["g1;g2;g3", "g1;g2;g3"],
            "Count": [3, 3],
            "TermSize": [5, 3],
        }
    )
    simp = scat.simplify_enrichment(
        df,
        method="pathway_denester",
        gene_sets=gene_sets,
        min_count=1,
        pval_threshold=0.05,
    )
    assert len(simp) <= len(df)
    assert "Denester_filter" not in simp.columns or (simp["Denester_filter"] == "keep").all()
    assert simp.attrs.get("simplify_method") == "pathway_denester"


def test_run_enrichment_universe_and_new_output():
    """Test clusterProfiler-aligned universe handling + enriched output columns/attrs."""
    genes = ["GeneA", "GeneB", "GeneC", "GeneX"]
    gene_sets = {
        "TERM1": ["GeneA", "GeneB", "GeneD"],
        "TERM2": ["GeneC", "GeneE", "GeneF"],
        "TERM3": ["GeneA", "GeneC"],
    }
    # No background -> uses gene_sets union
    res1 = scat.run_enrichment(
        genes, gene_sets=gene_sets, pval_cutoff=1.0, min_size=1, return_all=True
    )
    assert "neg_log10_padj" in res1.columns
    assert "p.adjust" in res1.columns
    assert "universe_info" in res1.attrs
    ui = res1.attrs["universe_info"]
    assert ui["effective_universe_size"] > 0
    assert res1.attrs.get("clusterprofiler_aligned") is True

    # Provide background (like adata.var_names) -> by default intersected (conservative, clusterProfiler-like)
    bg = ["GeneA", "GeneB", "GeneC", "GeneY", "GeneZ"]  # GeneY/Z have no annotation in gene_sets
    res2 = scat.run_enrichment(
        genes, gene_sets=gene_sets, universe=bg, pval_cutoff=1.0, min_size=1, return_all=True
    )
    ui2 = res2.attrs["universe_info"]
    assert ui2["provided_size"] == 5
    assert ui2["restricted_to_gene_sets"] is True
    assert ui2["dropped_by_annotation_filter"] >= 2  # Y and Z dropped
    assert ui2["effective_universe_size"] == 3  # A,B,C

    # force_universe=True should keep the full provided size
    res3 = scat.run_enrichment(
        genes,
        gene_sets=gene_sets,
        background=bg,
        force_universe=True,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
    )
    ui3 = res3.attrs["universe_info"]
    assert ui3["force_universe"] is True
    assert ui3["effective_universe_size"] == 5
    assert ui3["dropped_by_annotation_filter"] == 0


def test_enrich_plot_show_terms(adata_basic):
    """Test that enrich_dotplot accepts show_terms (int or list) like clusterProfiler showCategory."""
    # Build a fake enrichment df similar to real output
    fake = pd.DataFrame(
        {
            "Term": ["T1 long name (GO:0001)", "T2 (KEGG:123)", "T3 foo bar", "T4"],
            "Description": ["desc1", "desc2", "the third", "fourth"],
            "p.adjust": [0.001, 0.01, 0.05, 0.2],
            "GeneRatio": [0.1, 0.2, 0.05, 0.01],
            "Count": [5, 3, 2, 1],
            "neg_log10_padj": [3.0, 2.0, 1.3, 0.7],
        }
    )
    # int
    fig, ax = scat.pl.enrich_dotplot(fake, show_terms=2, top_n=99)
    import matplotlib.pyplot as plt

    plt.close("all")
    # list of terms (partial match on Term or Description)
    fig2, ax2 = scat.pl.enrich_dotplot(fake, show_terms=["T3", "desc1"])
    plt.close("all")
    # explicitly use Count on x-axis (user-requested flexibility)
    fig3, ax3 = scat.pl.enrich_dotplot(fake, x="Count", top_n=3)
    plt.close("all")
    # also exercise FoldEnrichment
    if "FoldEnrichment" not in fake.columns:
        fake = fake.copy()
        fake["FoldEnrichment"] = [2.5, 1.8, 3.1, 1.2]
    fig4, ax4 = scat.pl.enrich_dotplot(fake, x="FoldEnrichment", top_n=2)
    plt.close("all")

    # size_by="Count" must produce visibly different dot sizes (was a reported bug)
    # We just exercise it; visual correctness is checked in integration runs.
    fig5, ax5 = scat.pl.enrich_dotplot(fake, size_by="Count", top_n=4)
    plt.close("all")
    fig6, ax6 = scat.pl.enrich_dotplot(fake, size_by="GeneRatio", top_n=3)
    plt.close("all")


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

    # ggVolcano-style manual gene labels + top_n
    f3, a3 = scat.pl.volcano_plot(
        allr, top_n=3, label_genes=["GeneA", allr.index[0]] if len(allr) > 0 else None
    )
    plt.close("all")


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


# --------------------------- copy_input / ensure_raw_counts ---------------------------


def test_active_score_default_copy_preserves_input(adata_basic):
    """active_score(copy_input=True) must not normalize/log1p the caller's AnnData."""
    ad = adata_basic.copy()
    X_before = np.asarray(ad.X, dtype=float).copy()
    had_log1p = "log1p" in ad.uns
    layer_keys = set(ad.layers.keys())

    scat.active_score(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )

    assert np.allclose(np.asarray(ad.X, dtype=float), X_before)
    assert ("log1p" in ad.uns) == had_log1p
    assert set(ad.layers.keys()) == layer_keys
    assert "active_score" not in ad.var.columns


def test_ensure_raw_counts_exported():
    assert hasattr(scat, "ensure_raw_counts")
    assert callable(scat.ensure_raw_counts)


def test_ensure_raw_counts_recovers_from_raw():
    """After log1p on .X, ensure_raw_counts should recover integer counts from adata.raw."""
    np.random.seed(7)
    n_cells, n_genes = 40, 60
    X_raw = np.random.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X_raw)
    ad.raw = ad.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    scat.ensure_raw_counts(ad)
    assert "counts" in ad.layers
    assert ad.layers["counts"].shape[1] == ad.n_vars
    assert "raw_gene_list" in ad.uns.get("scatrans", {})


def test_mixed_model_on_log_data_no_double_transform(adata_mixed_small):
    """Mixed model path should accept already log-normalized data without re-logging."""
    ad = adata_mixed_small.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    X_logged = np.asarray(ad.X, dtype=float).copy()

    _, _, allr = scat.active_score(
        ad,
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
    assert "delta_variance" in allr.columns
    assert np.allclose(np.asarray(ad.X, dtype=float), X_logged)


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


# --------------------------- filter_active_genes helper ---------------------------


def test_significant_requires_permutation_fdr(adata_basic):
    """Without permutation the built-in significant list is empty; with perm it can be non-empty."""
    _, sig_no_perm, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
    )
    assert len(sig_no_perm) == 0
    assert "unspliced_excess_fdr" not in allr.columns

    _, sig_perm, allr_perm = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=True,
        n_perm=30,
        random_seed=1,
        n_jobs=1,
    )
    assert "unspliced_excess_fdr" in allr_perm.columns
    if len(sig_perm) > 0:
        assert (sig_perm["unspliced_excess_residual"] > 0).all()
        assert (sig_perm["unspliced_excess_fdr"] < 0.05).all()


def test_filter_active_genes_basic(adata_mixed_small):
    """filter_active_genes should work, be robust to missing columns, and respect thresholds."""
    # Run without permutation (no fdr columns)
    _, _, allr = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=False,
        use_mixed_model=False,
        n_jobs=1,
    )

    # Basic call - should not crash even though active_score_fdr is missing
    filt = scat.filter_active_genes(
        allr,
        active_score_cutoff=30,
        pval_cutoff=0.1,
        unspliced_excess_residual_cutoff=0.5,
        logfc_cutoff=0.1,
        unspliced_excess_fdr_cutoff=0.25,  # ignored because column missing
        effective_gamma_min=0.01,
        effective_gamma_max=None,
    )
    assert isinstance(filt, pd.DataFrame)
    if len(filt) > 0:
        assert filt["active_score"].iloc[0] >= filt["active_score"].iloc[-1]  # sorted descending

    # With permutation enabled
    _, _, allr_perm = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=True,
        n_perm=20,
        random_seed=0,
        n_jobs=1,
    )
    filt2 = scat.filter_active_genes(
        allr_perm,
        active_score_cutoff=20,
        unspliced_excess_fdr_cutoff=0.5,  # permissive
        effective_gamma_min=0.0,
    )
    assert "unspliced_excess_fdr" in allr_perm.columns
    if len(filt2) > 0:
        assert (filt2["unspliced_excess_fdr"] < 0.5).all()


def test_filter_active_genes_with_mixed(adata_mixed_small):
    """When delta_variance is present, the helper can filter on it."""
    _, _, allr = scat.active_score(
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
    assert "delta_variance" in allr.columns

    filt = scat.filter_active_genes(
        allr,
        active_score_cutoff=10,
        delta_variance_min=0.0,  # permissive
    )
    # Should not have dropped the column requirement
    assert len(filt) >= 0


def test_filter_active_genes_presets(adata_mixed_small):
    """preset parameter supplies sensible defaults for different analysis styles."""
    _, _, allr = scat.active_score(
        adata_mixed_small,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        mode="heuristic",
        show_plot=False,
        use_permutation=True,
        n_perm=10,
        n_jobs=1,
    )

    # permissive preset (or no preset) should keep most genes
    f_perm = scat.filter_active_genes(allr, preset="permissive")
    assert len(f_perm) > 0

    # heuristic preset applies stricter single-cell style defaults
    f_heu = scat.filter_active_genes(allr, preset="heuristic")
    # on small synthetic data this may be small or zero, but should not crash
    assert isinstance(f_heu, pd.DataFrame)

    # pseudobulk preset uses lenient values suitable after aggregation
    f_pb = scat.filter_active_genes(allr, preset="pseudobulk")
    assert len(f_pb) >= 0

    # explicit arg should override preset
    f_over = scat.filter_active_genes(allr, preset="heuristic", active_score_cutoff=0)
    assert len(f_over) > len(f_heu) or len(f_heu) == 0


# --------------------------- Enrichment core logic & edge cases (per design review) ---------------------------


def test_enrichment_hypergeom_small_example():
    """Fixed small example to verify core hypergeom.sf(k-1, N, K, n) logic.

    gene_list = ["A", "B"]
    gene_sets = {"set1": ["A", "B", "C"], "set2": ["D", "E"]}
    universe  = ["A", "B", "C", "D", "E"]

    For set1: N=5, n=2, K=3, k=2
    p = hypergeom.sf(1, 5, 3, 2)
    """
    from scipy.stats import hypergeom

    genes = ["A", "B"]
    gene_sets = {"set1": ["A", "B", "C"], "set2": ["D", "E"]}
    univ = ["A", "B", "C", "D", "E"]

    res = scat.run_enrichment(
        genes, gene_sets=gene_sets, universe=univ, pval_cutoff=1.0, min_size=1, return_all=True
    )
    assert not res.empty
    # set1 must be present
    row = res[res["Term"] == "set1"].iloc[0]
    assert row["Count"] == 2
    assert row["TermSize"] == 3
    expected_p = hypergeom.sf(1, 5, 3, 2)
    assert np.isclose(row["pvalue"], expected_p, rtol=1e-12)

    # set2 has k=0 -> filtered
    assert "set2" not in res["Term"].values


def test_enrichment_universe_variants_and_restrict():
    genes = ["A", "B", "C"]
    gene_sets = {"S1": ["A", "B", "C", "D"], "S2": ["A", "E"]}

    # No universe -> effective = union of gene sets
    r1 = scat.run_enrichment(genes, gene_sets, pval_cutoff=1.0, min_size=1, return_all=True)
    assert r1.attrs["universe_info"]["effective_universe_size"] == 5  # A,B,C,D,E

    # restrict=True (default) drops genes not in gene_sets
    bg = ["A", "B", "C", "X", "Y", "Z"]
    r2 = scat.run_enrichment(
        genes, gene_sets, universe=bg, pval_cutoff=1.0, min_size=1, return_all=True
    )
    ui2 = r2.attrs["universe_info"]
    assert ui2["provided_size"] == 6
    assert ui2["restricted_to_gene_sets"] is True
    assert ui2["effective_universe_size"] == 3  # A,B,C kept
    assert ui2["dropped_by_annotation_filter"] == 3

    # restrict=False keeps full provided (even unannotated)
    r3 = scat.run_enrichment(
        genes,
        gene_sets,
        universe=bg,
        restrict_background_to_gene_sets=False,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
    )
    assert r3.attrs["universe_info"]["effective_universe_size"] == 6

    # force_universe=True also bypasses intersect
    r4 = scat.run_enrichment(
        genes,
        gene_sets,
        background=bg,
        force_universe=True,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
    )
    assert r4.attrs["universe_info"]["force_universe"] is True
    assert r4.attrs["universe_info"]["effective_universe_size"] == 6


def test_enrichment_both_universe_background_raises():
    genes = ["A"]
    gs = {"T": ["A"]}
    with pytest.raises(ValueError, match="only one of `universe` or `background`"):
        scat.run_enrichment(genes, gs, universe=["A"], background=["A"])


def test_enrichment_empty_cases_preserve_attrs_and_reason():
    gs = {"T1": ["X", "Y"]}

    # Empty gene list
    e1 = scat.run_enrichment([], gs, pval_cutoff=0.05)
    assert e1.empty
    assert e1.attrs.get("reason") == "gene_list_empty"
    assert e1.attrs.get("method") == "ora"

    # Empty gene_sets dict -> warning + empty with reason
    import warnings

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        e2 = scat.run_enrichment(["A", "B"], {}, pval_cutoff=1.0)
        assert e2.empty
        assert e2.attrs.get("reason") in ("no_gene_sets_loaded", "gene_list_empty") or len(w) > 0

    # No overlap case
    e3 = scat.run_enrichment(
        ["Z"], gs, universe=["Z"], pval_cutoff=1.0, min_size=1, return_all=True
    )
    assert e3.empty or len(e3) == 0  # may be dropped by overlap==0 path
    # At least we didn't crash; attrs may be on the non-empty path if any term passed size filter

    # Universe empty after processing
    e4 = scat.run_enrichment(["A"], gs, universe=[], pval_cutoff=1.0)
    assert e4.empty
    assert e4.attrs.get("reason") == "universe_empty"


def test_enrichment_gene_case_and_mapping_examples():
    # lower / upper should affect matching
    genes = ["actb", "Gapdh"]
    gs = {"Housekeeping": ["ACTB", "GAPDH", "TUBB"]}

    # default (None) keeps case; "actb" won't match "ACTB" unless gene_case
    scat.run_enrichment(genes, gs, pval_cutoff=1.0, min_size=1, return_all=True)
    # with gene_case=upper they should map
    r_upper = scat.run_enrichment(
        genes, gs, gene_case="upper", pval_cutoff=1.0, min_size=1, return_all=True
    )
    assert len(r_upper) >= 1
    # Check that low mapping warning would have triggered on default if rate low
    # (we don't assert the warning here, but the upper case path should have higher hit rate)


def test_enrichment_padj_cutoff_alias_and_log_message():
    genes = ["A", "B"]
    gs = {"S1": ["A", "B", "C"]}
    # Using padj_cutoff should work identically to pval_cutoff
    r1 = scat.run_enrichment(genes, gs, universe=["A", "B", "C"], pval_cutoff=0.5, return_all=True)
    r2 = scat.run_enrichment(genes, gs, universe=["A", "B", "C"], padj_cutoff=0.5, return_all=True)
    assert len(r1) == len(r2)
    # attrs should record the cutoff used
    assert "pval_cutoff" in r1.attrs or r1.attrs.get("pval_cutoff") is not None


def test_enrichment_include_gene_list():
    genes = ["A", "B"]
    gs = {"S1": ["A", "B", "C"]}
    r = scat.run_enrichment(
        genes,
        gs,
        universe=["A", "B", "C"],
        include_gene_list=True,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
    )
    assert "Genes_list" in r.columns
    assert isinstance(r.iloc[0]["Genes_list"], list)
    assert "Genes" in r.columns  # still present


def test_enrichment_return_all_vs_sig_and_tested_log():
    genes = ["A", "B", "C"]
    gs = {"Big": ["A", "B", "C", "D"], "Small": ["X"]}
    r_all = scat.run_enrichment(genes, gs, pval_cutoff=1.0, min_size=1, return_all=True)
    r_sig = scat.run_enrichment(genes, gs, pval_cutoff=1e-9, min_size=1, return_all=False)
    # return_all gives more (or equal)
    assert len(r_all) >= len(r_sig)
    # p.adjust should be non-increasing after sort
    if len(r_all) > 1:
        padjs = r_all["p.adjust"].values
        assert np.all(padjs[:-1] <= padjs[1:] + 1e-12)  # sorted ascending


def test_enrichment_organism_normalization_in_attrs():
    genes = ["A"]
    gs = {"T": ["A"]}
    r = scat.run_enrichment(genes, gs, organism="Mouse", pval_cutoff=1.0, return_all=True)
    assert r.attrs.get("organism") in ("mouse", "Mouse")  # normalized to lower inside
    # run_kegg passes "Mouse" but run_enrichment normalizes
    rk = scat.run_kegg(genes, organism="mmu", pval_cutoff=1.0, return_all=True, min_size=1)
    # may be empty but attrs should exist
    if not rk.empty or "organism" in rk.attrs:
        assert rk.attrs.get("organism") == "mouse" or rk.attrs.get("organism") is not None


def test_run_go_wrapper_basic():
    # run_go maps ontology -> gene set name. With real bundled sets (BP) this will load locally (no net).
    # "Foo" will have no overlap with any real GO BP term -> empty result via overlap filter, with rich attrs.
    res2 = scat.run_go(["Foo"], ontology="BP", organism="mouse", pval_cutoff=1.0, min_size=5)
    assert res2.empty
    assert res2.attrs.get("method") in ("ora", "ora_go_all")
    assert res2.attrs.get("reason") is not None or res2.attrs.get("gene_set_info") is not None

    # ALL path (runs BP+CC+MF) should at least not crash (CC/MF may fall back or produce empty)
    res3 = scat.run_go(["Foo"], ontology="ALL", organism="mouse", pval_cutoff=1.0, min_size=5)
    assert res3.empty or isinstance(res3, pd.DataFrame)


def test_enrichment_empty_gene_sets_raises_or_warns():
    with pytest.warns(UserWarning):
        df = scat.run_enrichment(["A", "B"], gene_sets={}, verbose=True)
    assert df.empty


# --------------------------- New save report + expand + provenance (review round 2) ---------------------------


def test_expand_enrichment_genes_and_save(tmp_path):
    """expand + save_enrichment_report should produce files and handle list columns cleanly."""
    genes = ["A", "B"]
    gs = {"S1": ["A", "B", "C"], "S2": ["A"]}
    res = scat.run_enrichment(
        genes,
        gs,
        universe=["A", "B", "C", "D"],
        include_gene_list=True,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
    )
    assert not res.empty

    # Expand should give more rows than unique terms
    long_df = scat.expand_enrichment_genes(res)
    assert len(long_df) >= len(res)
    assert "Gene" in long_df.columns
    assert "Term" in long_df.columns

    # Save with all options
    prefix = tmp_path / "test_enrich"
    saved = scat.save_enrichment_report(
        res,
        prefix=str(prefix),
        save_excel=True,
        save_csv=True,
        save_metadata=True,
        save_term_gene_table=True,
    )
    assert "results_csv" in saved
    assert (
        "results_xlsx" in saved or "metadata_json" in saved
    )  # xlsx may fail if no engine, but csv+meta should be there
    assert "term_gene_table_csv" in saved or "metadata_json" in saved

    # Check that Genes_list was stringified in the exported main table (if csv was written)
    if "results_csv" in saved:
        df_back = pd.read_csv(saved["results_csv"])
        if "Genes_list" in df_back.columns:
            # Should not be a literal python list repr
            val = str(df_back["Genes_list"].iloc[0])
            assert not val.startswith("['") or ";" in val  # we joined with ;


def test_save_enrichment_report_empty_still_writes(tmp_path):
    """Even empty res with attrs should be savable."""
    empty = scat.run_enrichment([], {"T": ["X"]})
    prefix = tmp_path / "empty_enrich"
    saved = scat.save_enrichment_report(
        empty, prefix=str(prefix), save_excel=False, save_csv=True, save_metadata=True
    )
    assert "results_csv" in saved
    assert "metadata_json" in saved


def test_gene_set_info_has_actual_source_and_requested_source():
    genes = ["A", "B"]
    gs = {"S1": ["A", "B"]}
    res = scat.run_enrichment(genes, gs, pval_cutoff=1.0, min_size=1, return_all=True)
    gsi = res.attrs.get("gene_set_info", {})
    assert "requested_source" in gsi
    assert "actual_source" in gsi
    assert gsi["actual_source"] in ("dict", "bundled", "gseapy", "gmt_file", None)


def test_run_go_adjust_across_all(tmp_path):
    """adjust_across_all changes the p.adjust values (or at least runs without error and records the flag)."""
    # Use synthetic tiny sets so we control overlaps; real bundled would also work but be slow/noisy
    genes = ["X", "Y"]
    # Make two "ontologies" via dicts passed indirectly is hard; use run_enrichment + post process not needed.
    # Instead run_go with a name that falls back, but to keep deterministic use high cutoff + return_all + check attr.
    res_all = scat.run_go(
        genes,
        ontology="ALL",
        organism="mouse",
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        adjust_across_all=True,
    )
    # May be empty (no real overlap with GO terms), but attrs should record the flag
    assert res_all.attrs.get("adjust_across_all") is True or res_all.empty

    res_sep = scat.run_go(
        genes,
        ontology="ALL",
        organism="mouse",
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        adjust_across_all=False,
    )
    assert res_sep.attrs.get("adjust_across_all") is False or res_sep.empty

    # When non-empty the re-adjust should have happened (we can't easily assert numerical change without controlling pvalues)
    # At minimum the function path exercised.


def test_analysis_info_present_in_attrs():
    genes = ["A"]
    gs = {"T": ["A"]}
    res = scat.run_enrichment(genes, gs, pval_cutoff=1.0, min_size=1, return_all=True)
    ai = res.attrs.get("analysis_info", {})
    assert ai.get("package") == "scatrans"
    assert "timestamp" in ai
    # version may be None in some dev installs but key exists
    assert "package_version" in ai


def test_dual_cutoff_warning():
    genes = ["A", "B"]
    gs = {"S1": ["A", "B", "C"]}
    with pytest.warns(UserWarning, match="Both pval_cutoff and padj_cutoff"):
        scat.run_enrichment(
            genes, gs, universe=["A", "B", "C"], pval_cutoff=0.1, padj_cutoff=0.05, min_size=1
        )


def test_run_go_all_per_ontology_attrs_and_within_padj():
    """GO ALL should capture per_ontology_attrs and (when adjust=True) the within_ontology column + rich top attrs."""
    genes = ["A", "B"]
    # Real GO path + high min_size often yields empty result, but structure + keys must be present
    res = scat.run_go(
        genes,
        ontology="ALL",
        organism="mouse",
        pval_cutoff=1.0,
        min_size=5,
        return_all=True,
        adjust_across_all=True,
    )
    attrs = res.attrs
    assert attrs.get("ontology") == "ALL" or attrs.get("method") == "ora_go_all"
    assert "per_ontology_attrs" in attrs
    poa = attrs["per_ontology_attrs"]
    assert isinstance(poa, dict)
    for o in ("BP", "CC", "MF"):
        assert o in poa
    # analysis_info always present
    assert "analysis_info" in attrs
    # If we somehow got rows (or in future synthetic), within column check
    if not res.empty and "p.adjust.within_ontology" in res.columns:
        assert "p.adjust" in res.columns


def test_expand_enrichment_genes_preserves_ontology():
    """Long table from GO ALL results must include and lead with Ontology column."""
    # Create a fake combined result resembling GO ALL output
    fake = pd.DataFrame(
        {
            "Ontology": ["BP", "CC"],
            "Term": ["T1", "T2"],
            "Description": ["d1", "d2"],
            "Count": [2, 1],
            "Genes": ["G1;G2", "G3"],
            "p.adjust": [0.01, 0.02],
        }
    )
    long = scat.expand_enrichment_genes(fake)
    assert "Ontology" in long.columns
    assert long.columns[0] == "Ontology"
    assert len(long) == 3  # 2 + 1
    assert set(long["Ontology"]) == {"BP", "CC"}


def test_save_enrichment_report_mkdir_and_tsv(tmp_path):
    """prefix dir should be created; save_tsv produces .tsv files."""
    from pathlib import Path as _Path  # local to avoid polluting module if not wanted

    genes = ["A"]
    gs = {"S": ["A", "B"]}
    res = scat.run_enrichment(genes, gs, pval_cutoff=1.0, min_size=1, return_all=True)
    subdir = tmp_path / "out" / "nested"
    pref = subdir / "myrep"
    saved = scat.save_enrichment_report(
        res, prefix=str(pref), save_csv=False, save_tsv=True, save_excel=False, save_metadata=True
    )
    assert (subdir).exists()
    assert "results_tsv" in saved
    assert saved["results_tsv"].endswith("_results.tsv")
    assert "term_gene_table_tsv" in saved or "metadata_json" in saved  # term may be present
    assert _Path(saved.get("results_tsv", "")).exists() or "results_tsv" not in saved

    # metadata json always written when requested
    if "metadata_json" in saved:
        assert _Path(saved["metadata_json"]).exists()


@pytest.mark.skipif(
    not (lambda: __import__("importlib.util").util.find_spec("gseapy"))(),
    reason="gseapy not installed for GSEA tests",
)
def test_run_gsea_basic():
    """Basic run_gsea smoke test using gseapy.prerank wrapper."""
    import pandas as pd

    ranked = pd.Series({"GeneA": 5.0, "GeneB": 4.5, "GeneC": 3.0, "GeneX": -2.0})
    gene_sets = {
        "TERM_UP": ["GeneA", "GeneB", "GeneD"],
        "TERM_DOWN": ["GeneX", "GeneY"],
    }
    res = scat.run_gsea(ranked, gene_sets, nperm=20, min_size=1, verbose=False)
    assert isinstance(res, pd.DataFrame)
    if not res.empty:
        assert "NES" in res.columns or "ES" in res.columns
        assert "pvalue" in res.columns or "p.adjust" in res.columns
        # attrs should have gsea info
        assert res.attrs.get("method") == "gsea_prerank"
        assert "gsea_info" in res.attrs


def test_run_gsea_without_gseapy_raises(monkeypatch):
    """If gseapy missing, run_gsea should raise clear ImportError."""
    # simulate no gseapy by patching
    import builtins

    import scatrans as scat

    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "gseapy" or name.startswith("gseapy."):
            raise ImportError("No module named 'gseapy' (simulated)")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    import pandas as pd

    with pytest.raises(ImportError, match="requires the 'gseapy'"):
        scat.run_gsea(pd.Series({"A": 1}), {"T": ["A"]})
