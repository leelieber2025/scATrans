"""Tests for enrichment public API: compare, extract, bundled sets, kegg."""

import warnings

import pandas as pd
import pytest

import scatrans as scat


@pytest.fixture
def tiny_gene_sets():
    return {
        "TERM_A": ["GeneA", "GeneB", "GeneC", "GeneD"],
        "TERM_B": ["GeneB", "GeneC", "GeneE"],
        "TERM_C": ["GeneX", "GeneY"],
    }


def test_resolve_gene_set_name_historical_versions_passthrough():
    """Explicit Enrichr year versions must not be silently remapped to bundled 2026."""
    from scatrans.enrich._data import _resolve_gene_set_name

    assert _resolve_gene_set_name("KEGG_2021", "scatrans", "mouse") == "KEGG_2021"
    assert _resolve_gene_set_name("GO_Biological_Process_2023", "scatrans", "human") == (
        "GO_Biological_Process_2023"
    )
    # Base names still map to bundled defaults
    assert _resolve_gene_set_name("KEGG", "scatrans", "mouse") == "Mm_KEGG_2026"
    assert _resolve_gene_set_name("GO_BP", "scatrans", "human") == "Hs_GO_Biological_Process_2026"


def test_list_bundled_gene_sets():
    names = scat.list_bundled_gene_sets(verbose=False)
    assert isinstance(names, list)
    assert len(names) > 0
    assert any("GO" in n or "KEGG" in n for n in names)
    # Guard against stale fallback lists referencing non-shipped "ghost" files
    assert "GO_Biological_Process_scATrans.gmt" not in names
    assert "KEGG_scATrans.gmt" not in names
    # At least the 2026 organism sets should be present
    assert any(n.endswith("_2026.txt") for n in names)


def test_extract_gene_lists_single_df():
    df = pd.DataFrame(
        {"logFC": [1.2, -0.9, 0.1], "p_adj": [0.01, 0.02, 0.5]},
        index=["G1", "G2", "G3"],
    )
    out = scat.extract_gene_lists(df, logfc_cutoff=0.5, pval_cutoff=0.05)
    assert isinstance(out, dict)
    assert len(out) >= 1


def test_run_enrichment_warns_on_mechanism_class_table(tiny_gene_sets):
    """P0: ORA on mechanism_class-partitioned tables must warn (induction trap)."""
    df = pd.DataFrame(
        {
            "logFC": [2.0, 2.0, 2.0, 2.0],
            "mechanism_class": [
                "stabilization-driven",
                "stabilization-driven",
                "stabilization-driven",
                "stabilization-driven",
            ],
        },
        index=["GeneA", "GeneB", "GeneC", "GeneD"],
    )
    with pytest.warns(UserWarning, match="mechanism_class|induction"):
        scat.run_enrichment(
            df,
            gene_sets=tiny_gene_sets,
            universe=["GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneX"],
            verbose=False,
            min_size=1,
            max_size=50,
        )


def test_run_enrichment_allow_mechanism_class_ora_silences(tiny_gene_sets):
    df = pd.DataFrame(
        {
            "logFC": [2.0, 2.0, 2.0],
            "mechanism_class": ["transcription-driven"] * 3,
        },
        index=["GeneA", "GeneB", "GeneC"],
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        scat.run_enrichment(
            df,
            gene_sets=tiny_gene_sets,
            universe=["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"],
            verbose=False,
            min_size=1,
            max_size=50,
            allow_mechanism_class_ora=True,
        )


def test_extract_gene_lists_prefers_gene_column_over_range_index():
    """Scanpy-style DE tables use a 'gene'/'names' column with a default RangeIndex."""
    df = pd.DataFrame(
        {
            "gene": ["G_up", "G_down", "G_ns"],
            "logFC": [1.2, -0.9, 0.1],
            "p_adj": [0.01, 0.02, 0.5],
        },
    )
    out = scat.extract_gene_lists(df, logfc_cutoff=0.5, pval_cutoff=0.05, logfc_direction="up")
    assert out["contrast"] == ["G_up"]


def test_extract_gene_lists_separate_directions_single_df_matches_dict():
    df = pd.DataFrame(
        {"logFC": [1.0, -1.0, 2.0, -2.0], "p_adj": [0.01] * 4},
        index=["G1", "G2", "G3", "G4"],
    )
    single = scat.extract_gene_lists(
        df, logfc_direction="up", separate_directions=True, name_prefix="X"
    )
    multi = scat.extract_gene_lists({"X": df}, logfc_direction="up", separate_directions=True)
    assert single == multi
    assert single == {"X_up": ["G1", "G3"], "X_down": ["G2", "G4"]}


def test_extract_gene_lists_separate_directions_single_df_aligned():
    df = pd.DataFrame(
        {"logFC": [1.0, -1.0], "p_adj": [0.01, 0.01]},
        index=["GeneA", "GeneB"],
    )
    out = scat.extract_gene_lists(
        df,
        logfc_cutoff=0.5,
        pval_cutoff=0.05,
        logfc_direction="both",
        separate_directions=True,
    )
    assert out["up"] == ["GeneA"]
    assert out["down"] == ["GeneB"]


def test_extract_gene_lists_multi_and_separate_directions():
    df1 = pd.DataFrame(
        {"logFC": [1.0, -1.0], "p_adj": [0.01, 0.01]},
        index=["Up1", "Down1"],
    )
    df2 = pd.DataFrame(
        {"logFC": [0.8, -0.7], "p_adj": [0.02, 0.03]},
        index=["Up2", "Down2"],
    )
    out = scat.extract_gene_lists(
        {"A": df1, "B": df2},
        logfc_cutoff=0.5,
        pval_cutoff=0.05,
        logfc_direction="both",
        separate_directions=True,
    )
    assert any("_up" in k or "_down" in k for k in out)


def test_compare_enrichment_clusters_attrs_only_successful(tiny_gene_sets):
    """attrs['clusters'] must list only clusters that contributed rows, not skipped/failed/empty."""
    clusters = {
        "ClusterA": ["GeneA", "GeneB"],
        "ClusterB_empty": [],
        "ClusterC_badgenes": ["NOTAREALGENE123"],
    }
    res = scat.compare_enrichment(
        clusters,
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    if not res.empty:
        assert sorted(res["Cluster"].unique()) == ["ClusterA"]
    assert res.attrs["n_clusters"] == 1
    assert res.attrs["clusters"] == ["ClusterA"]
    assert len(res.attrs["clusters"]) == res.attrs["n_clusters"]
    per_cluster = res.attrs.get("per_cluster", {})
    assert "ClusterB_empty" in per_cluster
    assert "ClusterC_badgenes" in per_cluster


def test_compare_enrichment(tiny_gene_sets):
    clusters = {"C1": ["GeneA", "GeneB"], "C2": ["GeneC", "GeneE"]}
    res = scat.compare_enrichment(
        clusters,
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    assert isinstance(res, pd.DataFrame)
    if not res.empty:
        assert "Cluster" in res.columns
    assert res.attrs.get("method") == "compare_enrichment"


def test_compare_enrichment_adjust_across_clusters(tiny_gene_sets):
    clusters = {"G1": ["GeneA", "GeneB"], "G2": ["GeneB", "GeneC"]}
    res = scat.compare_enrichment(
        clusters,
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        adjust_across_clusters=True,
        return_all=True,
        verbose=False,
    )
    sc_meta = res.attrs.get("scatrans", {})
    assert sc_meta.get("adjust_across_clusters") is True
    # Default method is fdr_bh; must be recorded (not a silent hardcoded label only)
    mt = sc_meta.get("multiple_testing", {})
    assert mt.get("scope") == "all_clusters"
    assert mt.get("method") == "fdr_bh"


def test_compare_enrichment_adjust_across_clusters_honors_bonferroni(tiny_gene_sets):
    """Global re-adjust must use p_adjust_method, not always BH."""
    import numpy as np

    from scatrans.enrich._data import _apply_p_adjust

    clusters = {"G1": ["GeneA", "GeneB"], "G2": ["GeneB", "GeneC"]}
    res = scat.compare_enrichment(
        clusters,
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        adjust_across_clusters=True,
        return_all=True,
        p_adjust_method="bonferroni",
        verbose=False,
    )
    assert res.attrs.get("p_adjust_method") == "bonferroni"
    mt = res.attrs.get("scatrans", {}).get("multiple_testing", {})
    assert mt.get("method") == "bonferroni"
    if not res.empty and "pvalue" in res.columns:
        expected = _apply_p_adjust(res["pvalue"].to_numpy(dtype=float), method="bonferroni")
        np.testing.assert_allclose(
            res["p.adjust"].to_numpy(dtype=float), expected, rtol=1e-10, atol=1e-12
        )
        assert "p.adjust.within_cluster" in res.columns


def test_concat_compare_results_skips_none_and_empty_clusters():
    from scatrans.enrich import concat_compare_results

    df1 = pd.DataFrame({"Term": ["t1"], "p.adjust": [0.01], "Count": [5]})
    results = {
        "GroupA": df1,
        "GroupB_empty": pd.DataFrame(columns=["Term", "p.adjust", "Count"]),
        "GroupC_none": None,
    }
    combined = concat_compare_results(results)
    assert sorted(combined["Cluster"].unique()) == ["GroupA"]
    assert combined.attrs["n_clusters"] == 1
    assert combined.attrs["clusters"] == ["GroupA"]
    assert len(combined.attrs["clusters"]) == combined.attrs["n_clusters"]


def test_resolve_gseapy_weight_broad_semantics():
    from scatrans.enrich._data import _resolve_gseapy_weight

    assert _resolve_gseapy_weight() == 1.0
    assert _resolve_gseapy_weight(weighted_score_type="classic") == 0.0
    assert _resolve_gseapy_weight(weighted_score_type="weighted") == 1.0
    assert _resolve_gseapy_weight(weighted_score_type="unweighted") == 0.0


def test_compare_enrichment_attrs_clusters_match_returned_rows(tiny_gene_sets):
    """After global padj filter, attrs['clusters'] must match rows actually returned."""
    clusters = {"clusterA": ["GeneA", "GeneB"], "clusterB": ["GeneC"]}
    res = scat.compare_enrichment(
        clusters,
        gene_sets=tiny_gene_sets,
        padj_cutoff=1e-12,
        min_size=1,
        adjust_across_clusters=True,
        return_all=False,
        verbose=False,
    )
    if res.empty:
        assert res.attrs.get("clusters") == []
        assert res.attrs.get("n_clusters") == 0
    else:
        assert set(res.attrs["clusters"]) == set(res["Cluster"].astype(str).unique())
        assert res.attrs["n_clusters"] == len(res.attrs["clusters"])


def test_concat_compare_results_duplicate_cluster_names():
    from scatrans.enrich import concat_compare_results

    df1 = pd.DataFrame({"Term": ["t1"], "p.adjust": [0.01], "Count": [5]})
    df2 = pd.DataFrame({"Term": ["t2"], "p.adjust": [0.02], "Count": [3]})
    # dict keys cannot repeat; use list-of-tuples API for duplicate display names.
    combined = concat_compare_results([("dup", df1), ("dup", df2)])
    assert combined.attrs["n_clusters"] == 2
    assert combined.attrs["clusters"] == ["dup", "dup"]
    assert len(combined.attrs["clusters"]) == combined.attrs["n_clusters"]


def test_concat_compare_results(tiny_gene_sets):
    r1 = scat.run_enrichment(
        ["GeneA", "GeneB"],
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    r2 = scat.run_enrichment(
        ["GeneC"],
        gene_sets=tiny_gene_sets,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    combined = scat.concat_compare_results({"A": r1, "B": r2})
    assert "Cluster" in combined.columns or combined.empty
    assert combined.attrs.get("method") == "compare_concat"


@pytest.mark.slow
def test_run_kegg_bundled_smoke():
    res = scat.run_kegg(
        ["GeneA"],
        organism="mouse",
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    assert isinstance(res, pd.DataFrame)
    gsi = res.attrs.get("gene_set_info", {})
    assert gsi.get("requested_source") is not None or res.empty


def test_run_enrichment_with_adata_universe(adata_basic, tiny_gene_sets):
    ad = adata_basic.copy()
    scat.store_raw_counts(ad, layer="counts")
    res = scat.run_enrichment(
        ad.var_names[:5].tolist(),
        gene_sets=tiny_gene_sets,
        adata=ad,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    ui = res.attrs.get("universe_info", {})
    assert ui.get("effective_universe_size", 0) >= 0
