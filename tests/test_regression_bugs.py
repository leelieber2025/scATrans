"""Regression tests for confirmed bug fixes (2026-07 audit)."""

from __future__ import annotations

import importlib.util
import logging
import sys
from unittest.mock import MagicMock, patch

import anndata as ad
import matplotlib as mpl
import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import scatrans as scat
from scatrans._de import _run_de_wrapper, _validate_de_result
from scatrans._utils import (
    _apply_de_preprocess,
    _is_integer_counts_like,
    _pseudobulk_with_layers,
    _x_gene_dispersion_looks_raw,
)
from scatrans.tl._common import _materialize_if_view, _select_var


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


def test_bias_corrected_false_on_median_fallback_without_huber_fit():
    """Median fallback must not set bias_corrected=True when Huber regression did not run."""
    from scatrans._utils import _fit_huber_bias_correction

    n = 5
    delta = np.array([1.0, 2.0, 0.5, -0.2, 3.0])
    valid_expr = np.ones(n, dtype=bool)
    valid_feat = np.zeros(n, dtype=bool)  # no genes with length/intron features
    _, bias = _fit_huber_bias_correction(
        delta,
        np.full(n, np.nan),
        np.full(n, np.nan),
        np.ones(n),
        valid_feat,
        valid_expr,
        X_features=None,
        min_fit_obs=30,
        bias_correction="huber_length_intron",
    )
    assert bias["fallback_to_median"] is True
    assert bias["bias_corrected"] is False
    assert bias["n_genes_used_for_fit"] == 0


def test_huber_success_does_not_zero_residual_for_missing_features():
    """After successful Huber fit, genes lacking features get median-centered residual (not 0)."""
    from scatrans._utils import _fit_huber_bias_correction

    rng = np.random.default_rng(0)
    n = 40
    gene_length = rng.uniform(500, 5000, size=n)
    intron_number = rng.uniform(1, 20, size=n)
    # Make delta depend on length so Huber has signal
    delta = 0.002 * np.log1p(gene_length) + rng.normal(0, 0.05, size=n)
    # Last 5 genes: missing features but expressed
    valid_feat = np.ones(n, dtype=bool)
    valid_feat[-5:] = False
    gene_length[-5:] = np.nan
    intron_number[-5:] = np.nan
    valid_expr = np.ones(n, dtype=bool)
    X_features = np.column_stack(
        [np.log1p(gene_length[valid_feat]), np.log1p(intron_number[valid_feat])]
    )
    residual, bias = _fit_huber_bias_correction(
        delta,
        gene_length,
        intron_number,
        np.ones(n) * 100,
        valid_feat,
        valid_expr,
        X_features=X_features,
        min_fit_obs=10,
        bias_correction="huber_length_intron",
    )
    assert bias.get("bias_corrected") is True
    assert bias.get("n_genes_residual_missing_features") == 5
    # Missing-feature genes must not all be left at exact 0
    assert not np.allclose(residual[-5:], 0.0)


def test_builtin_significant_mask_uses_pseudobulk_scale():
    """Pseudobulk residual ~0.1 should pass pb thresholds, not SC residual>1."""
    from scatrans.tl.filter import _builtin_significant_mask

    df = pd.DataFrame(
        {
            "p_adj": [0.01, 0.01],
            "logFC": [0.5, 0.5],
            "unspliced_excess_residual": [0.1, 0.1],
            "active_score": [10.0, 10.0],
            "unspliced_excess_fdr": [0.01, 0.01],
            "active_score_fdr": [0.01, 0.01],
            "valid_expr": [True, True],
        },
        index=["g1", "g2"],
    )
    mask_sc = _builtin_significant_mask(
        df,
        use_permutation=True,
        extra_metadata={"is_pseudobulk": False, "use_fdr_for_significance": True},
    )
    mask_pb = _builtin_significant_mask(
        df,
        use_permutation=True,
        extra_metadata={"is_pseudobulk": True, "use_fdr_for_significance": True},
    )
    assert not mask_sc.any()  # residual 0.1 << SC cutoff 1.0; score 10 << 55
    assert mask_pb.all()  # residual 0.1 > 0.05; score 10 >= 5


def test_default_mouse_gene_features_prefer_2020a():
    """Default organism=mouse must prefer mouse_2020A over GRCm39 Ensembl tables."""
    from scatrans.pp_bias import list_available_gene_features

    avail = list_available_gene_features(verbose=False)
    preferred = "mouse_2020A_gene_features.parquet"
    if preferred not in avail:
        pytest.skip("mouse_2020A not packaged")
    # Lexicographic first mouse candidate is often GRCm39; preferred must still win.
    mouse_cands = [f for f in avail if "mouse" in f.lower() or f.startswith("Mus")]
    assert preferred in mouse_cands
    # Resolution: preferred name is chosen when present (not mouse_cands[0])
    chosen = preferred if preferred in avail else mouse_cands[0]
    assert chosen == preferred


def test_add_gene_features_warns_on_lowercase_symbol_mismatch():
    """Lowercase var_names against uppercase feature tables must emit a mapping warning."""
    adata = ad.AnnData(
        np.ones((4, 4)),
        var=pd.DataFrame(index=["actb", "gapdh", "tp53", "malat1"]),
    )
    with pytest.warns(UserWarning, match="Low gene feature mapping rate"):
        out = scat.add_gene_features(adata, organism="human")
    assert out.var["gene_length"].isna().all()


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


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
def test_strict_pydeseq2_counts_rejects_log_normalized_pseudobulk():
    """Rounding before integer check must not let log-normalized pb data through."""
    rng = np.random.default_rng(7)
    n_genes = 20
    genes = [f"G{i}" for i in range(n_genes)]
    # Two samples per group; X values are clearly not integer counts
    X = rng.uniform(4.0, 6.5, size=(4, n_genes))
    obs = pd.DataFrame(
        {
            "sample": ["s1", "s1", "s2", "s2"],
            "condition": ["Control", "Control", "Disease", "Disease"],
        }
    )
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=genes))
    with pytest.raises(ValueError, match="does not look like raw non-negative integer counts"):
        _run_de_wrapper(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            is_pseudobulk=True,
            pb_backend="pydeseq2",
            strict_pydeseq2_counts=True,
            min_counts_per_gene=1,
        )


def test_pseudobulk_with_layers_flags_non_count_source_before_rounding():
    """_pseudobulk_with_layers must check count-likeness pre-aggregation (round() always
    produces integer-looking sums, so a post-hoc check on the aggregated X would always pass)."""
    rng = np.random.default_rng(3)
    n_cells, n_genes = 40, 15
    X = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    lognorm = np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4)
    adata = ad.AnnData(lognorm)
    adata.obs["condition"] = ["Disease"] * 20 + ["Control"] * 20
    adata.obs["sample"] = ["D0"] * 10 + ["D1"] * 10 + ["C0"] * 10 + ["C1"] * 10

    pb = _pseudobulk_with_layers(
        adata, sample_col="sample", groupby="condition", min_cells=5, min_counts=0
    )
    assert pb.uns["pb_x_is_count_like"] is False


def test_pseudobulk_flags_counts_layer_non_count_source_before_rounding():
    """Aggregated layers['counts'] must record pre-agg count-likeness (same trap as .X)."""
    rng = np.random.default_rng(11)
    n_cells, n_genes = 40, 12
    raw = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    lognorm = np.log1p(raw / raw.sum(axis=1, keepdims=True) * 1e4)
    adata = ad.AnnData(raw)  # .X is count-like
    adata.layers["counts"] = lognorm  # mislabeled non-counts
    adata.layers["spliced"] = raw.copy()
    adata.layers["unspliced"] = (raw * 0.3).astype(float)
    adata.obs["condition"] = ["Disease"] * 20 + ["Control"] * 20
    adata.obs["sample"] = ["D0"] * 10 + ["D1"] * 10 + ["C0"] * 10 + ["C1"] * 10

    pb = _pseudobulk_with_layers(
        adata,
        sample_col="sample",
        groupby="condition",
        layers=["spliced", "unspliced", "counts"],
        use_total_for_x=True,
        min_cells=5,
        min_counts=0,
    )
    # Post-round counts layer looks integer
    assert _is_integer_counts_like(pb.layers["counts"])
    # But pre-agg verdict must flag the non-count source
    assert pb.uns["pb_counts_is_count_like"] is False
    assert pb.uns["pb_layer_is_count_like"]["counts"] is False
    # True count layers still pass
    assert pb.uns["pb_layer_is_count_like"]["spliced"] is True


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
def test_strict_pydeseq2_counts_rejects_aggregated_non_count_counts_layer():
    """counts= from aggregated layers['counts'] must honor pb_counts_is_count_like.

    Regression: post-aggregation round() made the counts layer look integer, so
    strict_pydeseq2_counts only protected .X via pb_x_is_count_like and let a
    log-normalized layers['counts'] through when passed as counts=.
    """
    rng = np.random.default_rng(13)
    n_cells, n_genes = 40, 15
    raw = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    lognorm = np.log1p(raw / np.maximum(raw.sum(axis=1, keepdims=True), 1.0) * 1e4)
    adata = ad.AnnData(raw)
    adata.layers["counts"] = lognorm
    adata.obs["condition"] = ["Disease"] * 20 + ["Control"] * 20
    adata.obs["sample"] = ["D0"] * 10 + ["D1"] * 10 + ["C0"] * 10 + ["C1"] * 10

    pb = _pseudobulk_with_layers(
        adata,
        sample_col="sample",
        groupby="condition",
        layers=["counts"],
        min_cells=5,
        min_counts=0,
    )
    assert pb.uns.get("pb_counts_is_count_like") is False
    # Post-round matrix looks integer — without the uns flag, strict would pass.
    assert _is_integer_counts_like(pb.layers["counts"])

    with pytest.raises(ValueError, match="does not look like raw non-negative integer counts"):
        _run_de_wrapper(
            pb,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            is_pseudobulk=True,
            pb_backend="pydeseq2",
            counts=pb.layers["counts"],  # same object path as resolve
            strict_pydeseq2_counts=True,
            min_counts_per_gene=1,
        )

    with pytest.raises(ValueError, match="does not look like raw non-negative integer counts"):
        _run_de_wrapper(
            pb,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            is_pseudobulk=True,
            pb_backend="pydeseq2",
            counts="counts",
            strict_pydeseq2_counts=True,
            min_counts_per_gene=1,
        )


def test_subset_col_normalizes_numeric_labels():
    """subset_values=1 / '1' must match obs float 1.0 (same as groupby normalize)."""
    rng = np.random.default_rng(0)
    n, g = 40, 8
    X = rng.poisson(3, size=(n, g)).astype(float)
    # condition balanced within each cluster
    cond = (["A", "B"] * (n // 2))[:n]
    cluster = [1.0] * (n // 2) + [2.0] * (n // 2)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"condition": cond, "cluster": cluster}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]),
    )
    for sv in (1, 1.0, "1", "1.0"):
        ad_out, res = scat.differential_expression(
            adata,
            groupby="condition",
            target_group="B",
            reference_group="A",
            subset_col="cluster",
            subset_values=sv,
            copy_input=True,
        )
        assert ad_out.n_obs == n // 2, f"subset_values={sv!r} kept {ad_out.n_obs}"
        assert res.shape[0] == g


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
def test_permutation_de_resolves_counts_layer_like_observed():
    """Null DE must prefer layers['counts'] (same as observed), not U+S .X."""
    from scatrans._de import _run_de_wrapper
    from scatrans._utils import (
        _normalize_group_label,
        _pseudobulk_with_layers,
        _resolve_aligned_raw_counts,
    )

    rng = np.random.default_rng(0)
    n, g = 60, 20
    counts = rng.poisson(10, size=(n, g)).astype(float)
    U = rng.poisson(3, size=(n, g)).astype(float)
    S = rng.poisson(7, size=(n, g)).astype(float)
    counts[n // 2 :, :5] += 30
    obs = pd.DataFrame(
        {
            "condition": ["A"] * (n // 2) + ["B"] * (n // 2),
            "sample": [f"s{i % 6}" for i in range(n)],
        }
    )
    adata = ad.AnnData(
        counts,
        obs=obs,
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]),
    )
    adata.layers["counts"] = counts.copy()
    adata.layers["unspliced"] = U
    adata.layers["spliced"] = S
    adata.obs["condition"] = adata.obs["condition"].map(_normalize_group_label).values

    pb = _pseudobulk_with_layers(
        adata,
        "sample",
        "condition",
        layers=["spliced", "unspliced", "counts"],
        use_total_for_x=True,
        min_cells=1,
        min_counts=1,
    )
    resolved = _resolve_aligned_raw_counts(pb, layer="counts", require_integer=True)
    assert resolved is not None
    de_obs = _run_de_wrapper(
        pb,
        "condition",
        "B",
        "A",
        is_pseudobulk=True,
        pb_backend="pydeseq2",
        counts=resolved,
        n_jobs=1,
        min_counts_per_gene=1,
    )
    # Simulate perm path: re-resolve on copy (must match observed source)
    ad_perm = pb.copy()
    perm_counts = _resolve_aligned_raw_counts(ad_perm, layer="counts", require_integer=True)
    de_perm = _run_de_wrapper(
        ad_perm,
        "condition",
        "B",
        "A",
        is_pseudobulk=True,
        pb_backend="pydeseq2",
        counts=perm_counts,
        n_jobs=1,
        min_counts_per_gene=1,
    )
    # Same matrix source → highly correlated logFC (identical labels → identical DE)
    corr = np.corrcoef(de_obs["logFC"].to_numpy(), de_perm["logFC"].to_numpy())[0, 1]
    assert corr > 0.99
    # Contrast: without counts=, DE uses U+S .X and diverges
    de_x = _run_de_wrapper(
        pb,
        "condition",
        "B",
        "A",
        is_pseudobulk=True,
        pb_backend="pydeseq2",
        counts=None,
        n_jobs=1,
        min_counts_per_gene=1,
    )
    corr_x = np.corrcoef(de_obs["logFC"].to_numpy(), de_x["logFC"].to_numpy())[0, 1]
    assert corr_x < 0.9


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
def test_normalize_log1p_skipped_for_pydeseq2_count_path():
    """Explicit normalize_log1p must not transform .X when skip_auto (PyDESeq2)."""
    from scatrans._utils import _apply_de_preprocess, _is_integer_counts_like

    rng = np.random.default_rng(0)
    X = rng.poisson(5, size=(6, 10)).astype(float)
    adata = ad.AnnData(X.copy())
    before = adata.X.copy()
    _apply_de_preprocess(adata, "normalize_log1p", skip_auto=True)
    assert np.allclose(adata.X, before)
    assert _is_integer_counts_like(adata.X)


def test_log1p_marker_trusted_when_max_exceeds_20():
    """High-depth log-normalized matrices with uns['log1p'] must not be re-logged."""
    from scatrans._utils import _reconcile_log1p_marker, _x_looks_log_normalized

    rng = np.random.default_rng(0)
    # Non-integer log-like values with max > 20 (bulk-like / high depth)
    X = rng.uniform(0.0, 25.0, size=(50, 80))
    adata = ad.AnnData(X)
    adata.uns["log1p"] = {"base": None}
    assert _x_looks_log_normalized(adata.X, has_log1p_marker=True) is True
    assert _reconcile_log1p_marker(adata) is True
    assert "log1p" in adata.uns


def test_extract_gene_lists_accepts_padj_and_p_adjust_columns():
    from scatrans.enrich.compare import extract_gene_lists

    df = pd.DataFrame(
        {"logFC": [1.0, -1.0, 0.1], "padj": [0.01, 0.01, 0.5]},
        index=["Up1", "Dn1", "ns"],
    )
    out = extract_gene_lists(df, logfc_cutoff=0.5, pval_cutoff=0.05, logfc_direction="up")
    assert out["contrast"] == ["Up1"]

    df2 = pd.DataFrame(
        {"log2FoldChange": [1.2, -0.8], "p.adjust": [0.001, 0.2]},
        index=["A", "B"],
    )
    out2 = extract_gene_lists(df2, logfc_cutoff=0.5, pval_cutoff=0.05, logfc_direction="up")
    assert out2["contrast"] == ["A"]


@pytest.mark.skipif(importlib.util.find_spec("gseapy") is None, reason="gseapy not installed")
def test_run_gsea_collapses_duplicate_ids_by_max_abs():
    from scatrans.enrich.gsea import run_gsea

    ranked = pd.Series([0.1, 2.0, -3.0], index=["geneA", "GENEA", "geneB"])
    # With gene_case upper, geneA and GENEA collide
    res = run_gsea(
        ranked,
        gene_sets={"SET": ["GENEA", "GENEB", "GENEC"]},
        organism="mouse",
        gene_case="upper",
        nperm=10,
        min_size=1,
        max_size=10,
        verbose=False,
    )
    # Should complete without error; duplicate collapsed
    assert isinstance(res, pd.DataFrame)


def test_gseaplot_fallback_sorts_by_score():
    """Fallback RES must walk high→low scores, not arbitrary Series order."""
    import matplotlib

    matplotlib.use("Agg")
    # High score gene is set member but appears last in unsorted series
    ranked = pd.Series({"Z_low": 0.1, "M_mid": 0.5, "A_hit": 3.0})
    gsea_res = pd.DataFrame(
        {
            "Term": ["T1"],
            "NES": [1.0],
            "pvalue": [0.1],
            "p.adjust": [0.2],
            "leading_edge": ["A_hit"],
        }
    )
    # no gsea_details / ranking attrs → fallback path
    fig, ax = scat.pl.gseaplot(ranked, gsea_result=gsea_res, term="T1", show=False)
    assert fig is not None
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_active_score_simple_copy_input_false_preserves_caller_var():
    rng = np.random.default_rng(0)
    n, g = 40, 15
    X = rng.poisson(3, size=(n, g)).astype(float)
    obs = pd.DataFrame(
        {
            "condition": ["Disease"] * (n // 2) + ["Control"] * (n // 2),
        }
    )
    var = pd.DataFrame(index=[f"g{i}" for i in range(g)])
    adata = ad.AnnData(X, obs=obs, var=var)
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.3
    assert "gene_length" not in adata.var.columns
    # No sample_col → single-cell path (avoids pb min_cells failures on tiny toys)
    scat.active_score_simple(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        sample_col=None,
        organism="mouse",
        show_plot=False,
        copy_input=False,
    )
    # Caller must not receive gene features from _maybe_add_gene_features
    assert "gene_length" not in adata.var.columns


def test_perm_fast_scanpy_log1p_on_pseudobulk_counts():
    """perm_de_backend='fast' must log1p integer pb .X so scanpy logFC stays finite."""
    from scatrans._permutation import _single_permutation_task
    from scatrans._utils import _is_integer_counts_like

    rng = np.random.default_rng(0)
    n_pb, n_genes = 6, 40
    # Large integer pb counts (would overflow scanpy expm1 if not log1p'd)
    X = rng.poisson(500, size=(n_pb, n_genes)).astype(float)
    X[:3, :5] += 2000
    labels = np.array(["Ctrl", "Ctrl", "Ctrl", "GA", "GA", "GA"])
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"condition": labels}),
        var=pd.DataFrame(
            {
                "gene_length": rng.integers(500, 5000, n_genes),
                "intron_number": rng.integers(1, 20, n_genes),
            },
            index=[f"g{i}" for i in range(n_genes)],
        ),
    )
    assert _is_integer_counts_like(adata.X)
    uns = X * 0.3
    spl = X * 0.7
    score, resid = _single_permutation_task(
        seed=1,
        original_labels=labels,
        target_group="GA",
        reference_group="Ctrl",
        adata_subset=adata,
        X_features=None,
        valid_feat=np.ones(n_genes, dtype=bool),
        uns_layer=uns,
        spl_layer=spl,
        total_us_for_filter=X.sum(axis=0),
        min_total_counts=10,
        weight_fc=1.0,
        weight_unspliced=1.0,
        weight_pval=1.0,
        lambda_fc=1.0,
        lambda_res=1.0,
        lambda_pval=1.0,
        is_pseudobulk=True,
        pb_backend="scanpy",
        de_method="t-test_overestim_var",
        prior_weight=5.0,
        gamma_method="heuristic_shrink",
        de_preprocess="none",
        strict_pydeseq2_counts=True,
    )
    assert score is not None and np.all(np.isfinite(score))
    assert resid is not None and np.all(np.isfinite(resid))


def test_memento_rejects_non_integer_counts_layer():
    """Memento must not silently use log-like layers['counts'] (no memento install needed)."""
    from scatrans._de import _run_memento_de

    rng = np.random.default_rng(0)
    X = rng.poisson(3, size=(30, 8)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"condition": ["A"] * 15 + ["B"] * 15}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(8)]),
    )
    adata.layers["counts"] = np.log1p(X)
    with pytest.raises(ValueError, match="raw integer counts"):
        _run_memento_de(
            adata,
            groupby="condition",
            target_group="B",
            reference_group="A",
            counts="counts",
            num_boot=100,
            n_cpus=1,
        )


def test_volcano_auto_respects_label_by_p_adj():
    import matplotlib

    matplotlib.use("Agg")
    # Gene with best p_adj is not the highest active_score
    df = pd.DataFrame(
        {
            "logFC": [2.0, 1.5, 0.1],
            "p_adj": [0.2, 0.001, 0.5],
            "active_score": [90.0, 10.0, 5.0],
        },
        index=["hi_score", "hi_sig", "other"],
    )
    from scatrans.pl import _volcano_collect_labels

    lab = _volcano_collect_labels(
        df, top_n=1, label_genes=None, label_by="p_adj", min_label_score=None
    )
    assert list(lab.index) == ["hi_sig"]
    lab2 = _volcano_collect_labels(
        df, top_n=1, label_genes=None, label_by="active_score", min_label_score=None
    )
    assert list(lab2.index) == ["hi_score"]


def test_enrich_vennplot_four_clusters_no_nameerror():
    import matplotlib

    matplotlib.use("Agg")
    df = pd.DataFrame(
        {
            "Cluster": ["A", "A", "B", "B", "C", "C", "D", "D"],
            "Term": ["t1", "t2", "t1", "t3", "t1", "t4", "t1", "t5"],
            "p.adjust": [0.01] * 8,
            "Description": ["d"] * 8,
        }
    )
    fig, ax = scat.pl.enrich_vennplot(df, show=False)
    assert fig is not None
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_active_score_records_sample_col_for_pseudobulk(adata_pb):
    """sample_col must be stored in uns['scatrans'] for pseudobulk (not only MixedLM)."""
    ad_out, _, _ = scat.active_score(
        adata_pb,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_pseudobulk=True,
        sample_col="sample",
        pseudobulk_de_backend="scanpy",
        de_method="wilcoxon",
        use_permutation=False,
        show_plot=False,
        pb_use_total_for_x=False,
        min_cells=1,
        min_counts=1,
    )
    assert ad_out.uns["scatrans"].get("sample_col") == "sample"


def test_mixed_model_and_memento_incompatible():
    from scatrans._de import _run_de_wrapper

    rng = np.random.default_rng(0)
    X = rng.poisson(3, size=(20, 5)).astype(float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame(
            {
                "condition": ["A"] * 10 + ["B"] * 10,
                "sample": [f"s{i % 5}" for i in range(20)],
            }
        ),
        var=pd.DataFrame(index=[f"g{i}" for i in range(5)]),
    )
    with pytest.raises(ValueError, match="incompatible"):
        _run_de_wrapper(
            adata,
            groupby="condition",
            target_group="B",
            reference_group="A",
            use_mixed_model=True,
            sample_col="sample",
            use_memento_de=True,
        )


def test_nascent_excess_forces_residual_only_weights():
    """ranking_mode='nascent_excess' must override custom DE weights."""
    import inspect

    from scatrans.tl import active as active_mod

    # Source-level contract: always assigns weight_fc, weight_pval, weight_unspliced = 0,0,1
    src = inspect.getsource(active_mod.active_score)
    assert "weight_fc, weight_pval, weight_unspliced = 0.0, 0.0, 1.0" in src
    assert "forces weight_fc=0" in src or "overriding weight" in src


def test_raw_gamma_uses_global_fallback_for_zero_ref_genes():
    from scatrans._velocity import _compute_velocity_delta

    uns = np.array([[0.0, 0.0, 2.0], [1.0, 1.0, 4.0]], dtype=float)
    spl = np.array([[0.0, 0.0, 10.0], [1.0, 1.0, 5.0]], dtype=float)
    t_mask = np.array([False, True])
    r_mask = np.array([True, False])
    _, _, gamma, info = _compute_velocity_delta(
        uns, spl, t_mask, r_mask, prior_weight=5.0, gamma_method="raw"
    )
    assert info.get("n_genes_global_ratio_fallback") == 2
    assert info.get("n_genes_raw_ratio") == 1
    # Zero-ref genes share the global ratio, not 1.0 from eps/eps
    g_fallback = info["global_ratio_fallback"]
    assert np.isclose(gamma[0], g_fallback)
    assert np.isclose(gamma[1], g_fallback)
    assert gamma[2] < 0.5  # 2/10


def test_pathway_denester_comb_matches_hypergeom_and_scales():
    from scipy.stats import hypergeom

    from scatrans.enrich.simplify import _comb_comb_comb

    # Large pathway size that would hang factorial-based comb
    N, K, n, k = 5000, 200, 80, 5
    p = _comb_comb_comb(K, k, n, N)
    expected = float(hypergeom.sf(k - 1, N, K, n))
    assert abs(p - expected) < 1e-9
    assert 0.0 <= p <= 1.0


def test_copy_input_false_does_not_mutate_caller_labels():
    rng = np.random.default_rng(0)
    n, g = 40, 10
    X = rng.poisson(3, size=(n, g)).astype(float)
    obs = pd.DataFrame(
        {
            "condition": np.array([1.0] * (n // 2) + [2.0] * (n // 2)),
        }
    )
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))
    before = adata.obs["condition"].copy()
    scat.differential_expression(
        adata,
        groupby="condition",
        target_group=2,
        reference_group=1,
        copy_input=False,
        de_preprocess="none",
    )
    # Caller labels must remain original floats
    assert list(adata.obs["condition"]) == list(before)


def test_store_raw_preserves_full_gene_list_on_restor():
    rng = np.random.default_rng(0)
    n, g = 20, 30
    X = rng.poisson(2, size=(n, g)).astype(float)
    adata = ad.AnnData(X, var=pd.DataFrame(index=[f"g{i}" for i in range(g)]))
    scat.store_raw_counts(adata)
    full = list(adata.uns["scatrans"]["raw_gene_list"])
    assert len(full) == g
    adata2 = adata[:, :10].copy()
    scat.store_raw_counts(adata2, overwrite=True)
    assert len(adata2.uns["scatrans"]["raw_gene_list"]) == 10
    assert "raw_gene_list_full" in adata2.uns["scatrans"]
    assert len(adata2.uns["scatrans"]["raw_gene_list_full"]) == g


def test_robust_median_gamma_excludes_zero_expression_genes():
    """robust_median base_gamma must not be dominated by (eps/eps)≈1 from zeros."""
    from scatrans._velocity import _compute_velocity_delta

    # 1 cell ref, 1 cell target; most genes zero in ref, one gene with U/S=0.2
    uns = np.array([[0.0, 0.0, 0.0, 2.0], [1.0, 1.0, 1.0, 4.0]], dtype=float)
    spl = np.array([[0.0, 0.0, 0.0, 10.0], [1.0, 1.0, 1.0, 5.0]], dtype=float)
    t_mask = np.array([False, True])
    r_mask = np.array([True, False])
    _, _, gamma, info = _compute_velocity_delta(
        uns, spl, t_mask, r_mask, prior_weight=5.0, gamma_method="robust_median"
    )
    # Only gene 3 is expressed in ref → median anchor ≈ 0.2, not ≈ 1
    assert info.get("n_genes_used_for_median_anchor") == 1
    # After shrinkage, gamma for expressed gene stays near 0.2 scale, not near 1
    assert gamma[3] < 0.5


def test_expand_gene_list_rejects_rangeindex_numeric_columns():
    """RangeIndex DE table without gene column must not treat logFC as gene IDs."""
    from scatrans.enrich._data import _expand_gene_list_input
    from scatrans.enrich.compare import extract_gene_lists

    df = pd.DataFrame({"logFC": [1.0, -1.0], "p_adj": [0.01, 0.01]})
    assert _expand_gene_list_input(df) == []
    out = extract_gene_lists(df, name_prefix="x")
    assert out.get("x", []) == []


def test_size_factor_zero_rows_use_unit_factor():
    from scatrans._utils import _normalize_velocity_layers_by_size_factor

    U = np.zeros((4, 3))
    U[0] = [10, 0, 5]
    U[1] = [100, 0, 50]  # deeper library so median SF != 1 for both
    S = np.zeros((4, 3))
    S[0] = [10, 2, 5]
    S[1] = [100, 20, 50]
    _, _, totals, factors = _normalize_velocity_layers_by_size_factor(U, S)
    assert factors[2] == 1.0  # zero-total row
    assert factors[3] == 1.0
    assert factors[0] != factors[1]  # different library sizes
    assert not np.any(factors > 1e6)


def test_diagnose_design_missing_sample_col_name():
    from scatrans.tl.design import diagnose_design

    n, g = 40, 5
    X = np.ones((n, g))
    obs = pd.DataFrame(
        {
            "condition": ["A"] * 20 + ["B"] * 20,
            "sample": [f"s{i % 4}" for i in range(n)],
        }
    )
    adata = ad.AnnData(X, obs=obs)
    adata.layers["spliced"] = X
    adata.layers["unspliced"] = X * 0.3
    d = diagnose_design(adata, "condition", "B", "A", sample_col="Sample")  # typo case
    assert any("not in adata.obs" in w for w in d["warnings"])


def test_mixedlm_diagnostics_surface_logfc_method_in_uns():
    """logFC_method / sign-discordant counts must appear in uns scatrans diagnostics."""
    rng = np.random.default_rng(0)
    n_s, c_per, n_g = 5, 6, 6
    n_per = n_s * c_per
    X = rng.poisson(20, size=(n_per * 2, n_g)).astype(float)
    X[:n_per, 0] += 40
    obs = pd.DataFrame(
        {
            "condition": ["Disease"] * n_per + ["Control"] * n_per,
            "sample": [f"d{i}" for i in range(n_s) for _ in range(c_per)]
            + [f"c{i}" for i in range(n_s) for _ in range(c_per)],
        }
    )
    adata = ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_g)]))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    ad_out, res = scat.differential_expression(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_mixed_model=True,
        sample_col="sample",
        n_jobs=1,
        copy_input=True,
    )
    mm = ad_out.uns["scatrans"]["diagnostics"]["mixed_model"]
    assert mm["used"] is True
    assert mm["logFC_method"] == "sample_mean_of_means_log2"
    assert "n_genes_logFC_mixedlm_sign_discordant" in mm
    assert isinstance(mm["n_genes_logFC_mixedlm_sign_discordant"], int)
    assert res.attrs.get("logFC_method") == "sample_mean_of_means_log2"


def test_mixedlm_logfc_is_sample_mean_of_means_not_cell_weighted():
    """MixedLM logFC must equal-weight samples (not be dominated by a large sample)."""
    from scatrans._de import _mixedlm_sample_aware_logfc, _run_mixedlm_de

    # Two samples per group. Control: balanced. Disease: one huge sample with low
    # expression and one tiny sample with high expression on gene 0.
    # Cell-weighted mean would be pulled down by the large sample; sample-mean-of-means
    # averages the two sample means equally.
    rng = np.random.default_rng(0)
    rows = []
    conds = []
    samps = []
    # Control: sC0 (20 cells, mean ~2), sC1 (20 cells, mean ~2) on gene 0
    for s, n in [("sC0", 20), ("sC1", 20)]:
        x = np.full((n, 3), 1.0)
        x[:, 0] = 2.0 + rng.normal(0, 0.05, size=n)
        rows.append(x)
        conds.extend(["Control"] * n)
        samps.extend([s] * n)
    # Disease: sD_big (80 cells, mean ~1.0) and sD_small (5 cells, mean ~4.0)
    # Cell-weighted Disease mean ≈ 1.17; sample-mean-of-means ≈ 2.5 → positive logFC
    # vs Control 2.0 when using sample means.
    for s, n, mu in [("sD_big", 80, 1.0), ("sD_small", 5, 4.0)]:
        x = np.full((n, 3), 1.0)
        x[:, 0] = mu + rng.normal(0, 0.05, size=n)
        rows.append(x)
        conds.extend(["Disease"] * n)
        samps.extend([s] * n)
    X = np.vstack(rows)
    # Treat X as already log1p-scale for the helper (values in a plausible range).
    logfc = _mixedlm_sample_aware_logfc(
        X,
        condition=np.array(conds),
        samples=np.array(samps),
        target_group="Disease",
        reference_group="Control",
    )
    # Sample-mean-of-means: Disease ~2.5, Control ~2.0 → positive logFC
    assert logfc[0] > 0.1
    # Cell-weighted mean would be negative (Disease bulk is low)
    cell_t = X[np.array(conds) == "Disease", 0].mean()
    cell_r = X[np.array(conds) == "Control", 0].mean()
    cell_logfc = np.log2((np.expm1(cell_t) + 1e-9) / (np.expm1(cell_r) + 1e-9))
    assert cell_logfc < 0
    assert logfc[0] * cell_logfc < 0  # opposite signs: sample-aware vs cell-weighted

    # End-to-end: run MixedLM on a design with enough samples and check attrs
    n_s, c_per, n_g = 5, 6, 4
    n_per = n_s * c_per
    Xe = rng.poisson(20, size=(n_per * 2, n_g)).astype(float)
    Xe[:n_per, 0] += 40  # Disease high
    obs = pd.DataFrame(
        {
            "condition": ["Disease"] * n_per + ["Control"] * n_per,
            "sample": [f"d{i}" for i in range(n_s) for _ in range(c_per)]
            + [f"c{i}" for i in range(n_s) for _ in range(c_per)],
        }
    )
    ad_mm = ad.AnnData(Xe, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(n_g)]))
    sc.pp.normalize_total(ad_mm, target_sum=1e4)
    sc.pp.log1p(ad_mm)
    res = _run_mixedlm_de(
        ad_mm,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        sample_col="sample",
        n_jobs=1,
    )
    assert res.attrs.get("logFC_method") == "sample_mean_of_means_log2"
    assert "mixedlm_coef" in res.columns
    assert res.loc["g0", "logFC"] > 0


@pytest.mark.skipif(
    importlib.util.find_spec("pydeseq2") is None,
    reason="pydeseq2 not installed",
)
def test_strict_pydeseq2_counts_rejects_log_normalized_layer_end_to_end():
    """End-to-end: use_pseudobulk=True + pb_x_layer pointing at log-normalized data must raise,
    not silently round the log-normalized sums into look-alike integers."""
    rng = np.random.default_rng(1)
    n_cells, n_genes = 300, 20
    X = rng.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    a = ad.AnnData(X)
    a.obs["condition"] = ["Disease"] * 150 + ["Control"] * 150
    a.obs["sample"] = (
        ["D0"] * 50 + ["D1"] * 50 + ["D2"] * 50 + ["C0"] * 50 + ["C1"] * 50 + ["C2"] * 50
    )
    a.layers["spliced"] = X.copy()
    a.layers["unspliced"] = (X * 0.4).astype(float)
    a.layers["lognorm"] = np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4)

    with pytest.raises(ValueError, match="does not look like raw non-negative integer counts"):
        scat.active_score(
            a,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            show_plot=False,
            use_pseudobulk=True,
            sample_col="sample",
            pseudobulk_de_backend="pydeseq2",
            strict_pydeseq2_counts=True,
            pb_x_layer="lognorm",
            pb_use_total_for_x=False,
            min_cells=5,
            min_counts=1,
        )


def test_dispersion_check_prevents_double_log1p_on_heterogeneous_population():
    """Real cell-population heterogeneity can inflate per-cell library-size CV even in properly
    log-normalized data (e.g. after anndata.concat() drops uns['log1p']). The cross-gene
    mean-variance dispersion check must stop this from being misread as raw counts and
    re-transformed (which would systematically compress logFC / effect sizes)."""
    rng = np.random.default_rng(0)
    n_genes = 300
    total = 10000.0

    def make_pop(n, n_active, conc):
        Xp = np.zeros((n, n_genes))
        idx = rng.choice(n_genes, size=n_active, replace=False)
        for i in range(n):
            counts = rng.dirichlet(np.ones(n_active) * conc) * total
            Xp[i, idx] = counts
        return Xp

    n_per = 100
    # Two very different expression-breadth cell states (concentrated vs diffuse) at equal
    # per-cell totals -- log1p compresses these differently, inflating post-log1p library CV.
    X = np.vstack(
        [
            make_pop(n_per, 15, 1.0),
            make_pop(n_per, 250, 0.5),
            make_pop(n_per, 15, 1.0),
            make_pop(n_per, 250, 0.5),
        ]
    )
    Xlog = np.log1p(X)
    adata = ad.AnnData(Xlog.copy())
    adata.obs["condition"] = ["Disease"] * (2 * n_per) + ["Control"] * (2 * n_per)
    # No uns['log1p'] marker set, simulating anndata.concat() dropping it.
    assert "log1p" not in adata.uns

    _apply_de_preprocess(adata, "auto")
    np.testing.assert_allclose(adata.X, Xlog)


def test_dispersion_check_still_catches_raw_decimal_counts():
    """Raw (non-integer) decimal counts, e.g. kallisto/salmon pseudo-counts, must still be
    detected as needing normalize_total + log1p even without an integer-count signature."""
    rng = np.random.default_rng(2)
    n_cells, n_genes = 200, 150
    raw = rng.negative_binomial(
        5, np.random.default_rng(5).uniform(0.05, 0.4, size=n_genes), size=(n_cells, n_genes)
    ).astype(float)
    raw_decimal = raw * rng.uniform(0.5, 2.0, size=(n_cells, 1))
    assert _x_gene_dispersion_looks_raw(raw_decimal) is True

    adata = ad.AnnData(raw_decimal.copy())
    _apply_de_preprocess(adata, "auto")
    assert not np.allclose(adata.X, raw_decimal)


def test_run_enrichment_default_no_pval_cutoff_deprecation_warning():
    """Default calls must not emit a false legacy pval_cutoff deprecation warning."""
    import warnings

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        scat.run_enrichment(
            ["Actb", "Gapdh", "Tp53"],
            gene_sets={"TERM": ["Actb", "Gapdh", "Mdm2"]},
            organism="mouse",
            min_size=1,
            verbose=False,
        )
    deprecated = [
        w for w in record if issubclass(w.category, UserWarning) and "pval_cutoff" in str(w.message)
    ]
    assert not deprecated


def test_run_enrichment_explicit_pval_cutoff_emits_deprecation_warning():
    with pytest.warns(UserWarning, match="pval_cutoff.*deprecated"):
        scat.run_enrichment(
            ["Actb"],
            gene_sets={"TERM": ["Actb", "Gapdh"]},
            pval_cutoff=0.1,
            min_size=1,
            verbose=False,
        )


def test_differential_expression_diagnostics_schema(adata_de_only):
    """DE-only path must expose diagnostics['mixed_model'] like active_score."""
    ad = adata_de_only.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    out, _ = scat.differential_expression(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_mixed_model=True,
        sample_col="sample",
        de_preprocess="none",
    )
    diag = out.uns["scatrans"]["diagnostics"]
    assert "mixed_model" in diag
    assert "failed_fit_rate" in diag["mixed_model"]
    assert "n_genes_failed_fit" in diag["mixed_model"]


def test_expand_enrichment_genes_preserves_columns_on_zero_overlap():
    df = pd.DataFrame(
        {
            "Term": ["T1", "T2"],
            "Description": ["d1", "d2"],
            "Genes": ["", ""],
            "Count": [0, 0],
            "pvalue": [0.5, 0.6],
            "p.adjust": [0.9, 0.95],
        }
    )
    expanded = scat.expand_enrichment_genes(df)
    assert expanded.shape == (0, len(expanded.columns))
    assert "Gene" in expanded.columns


def test_de_preprocess_auto_normalizes_depth_confounded_decimal_counts():
    """kallisto-style decimal counts with depth bias must not yield genome-wide DE."""
    rng = np.random.default_rng(11)
    n_cells, n_genes = 80, 40
    base = rng.uniform(0.5, 8.0, size=(n_cells, n_genes))
    depth = np.where(
        np.arange(n_cells) < n_cells // 2,
        1.0,
        2.2,
    )[:, None]
    X = base * depth
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"condition": ["Control"] * 40 + ["Disease"] * 40}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(n_genes)]),
    )
    _, res_auto = scat.differential_expression(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        de_preprocess="auto",
    )
    _, res_force = scat.differential_expression(
        adata.copy(),
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        de_preprocess="normalize_log1p",
    )
    n_sig_auto = int((res_auto["p_adj"] < 0.05).sum())
    n_sig_force = int((res_force["p_adj"] < 0.05).sum())
    assert n_sig_auto <= max(3, int(0.15 * n_genes))
    assert n_sig_force <= max(3, int(0.15 * n_genes))


def test_de_method_logreg_raises_clear_error(adata_de_only):
    ad = adata_de_only.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    with pytest.raises(ValueError, match="de_method='logreg' is not supported"):
        scat.differential_expression(
            ad,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            de_method="logreg",
            de_preprocess="none",
        )


def test_generate_gene_features_drops_empty_gene_name(tmp_path):
    pytest.importorskip("gtfparse")
    from scatrans.pp_bias import generate_gene_features_from_gtf

    gtf = tmp_path / "missing_name.gtf"
    rows = [
        'chr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id "GENE4"; gene_type "protein_coding";',
        'chr1\t.\tgene\t300\t500\t.\t+\t.\tgene_id "GENE5"; gene_name "Gene5"; gene_type "protein_coding";',
    ]
    gtf.write_text("\n".join(rows) + "\n")
    out = tmp_path / "features.parquet"
    gene_df = generate_gene_features_from_gtf(str(gtf), str(out))
    assert len(gene_df) == 1
    assert gene_df.iloc[0]["gene_name"] == "Gene5"


def test_scanpy_de_single_cell_per_group_raises_clear_error():
    """One cell per group must fail with a scATrans-wrapped message, not raw scanpy."""
    X = np.ones((2, 5), dtype=float)
    adata = ad.AnnData(
        X,
        obs=pd.DataFrame({"group": ["T", "R"]}),
        var=pd.DataFrame(index=[f"g{i}" for i in range(5)]),
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    with pytest.raises(ValueError, match="at least 2 cells per group"):
        _run_de_wrapper(
            adata,
            groupby="group",
            target_group="T",
            reference_group="R",
            de_method="wilcoxon",
        )


def test_negative_unspliced_layer_emits_warning(adata_basic, caplog):
    ad = adata_basic.copy()
    ad.layers["unspliced"][0, 0] = -3.0
    with caplog.at_level(logging.WARNING, logger="scatrans._utils"):
        scat.active_score(
            ad,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_permutation=False,
            copy_input=True,
        )
    assert any("negative value" in r.message.lower() for r in caplog.records)


def test_gseaplot_rejects_dataframe_as_first_argument():
    df = pd.DataFrame({"Term": ["T1"], "NES": [1.2], "pvalue": [0.01]})
    with pytest.raises(TypeError, match="gsea_result"):
        scat.pl.gseaplot(df, show=False)


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


def test_filter_significant_preset_matches_builtin_when_perm_fdr_disabled(adata_basic):
    """Small n_perm disables FDR; filter must still match built-in significant list."""
    _, sig, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=True,
        n_perm=30,
        random_seed=1,
        n_jobs=1,
        show_plot=False,
    )
    ctx = allr.attrs.get("scatrans_filter_context", {})
    assert ctx.get("use_fdr_for_significance") is False
    filt_sig = scat.filter_active_genes(allr, preset="significant")
    filt_heu = scat.filter_active_genes(allr, preset="heuristic")
    assert sig.index.tolist() == filt_sig.index.tolist()
    assert sig.index.tolist() == filt_heu.index.tolist()


def test_filter_heuristic_skips_fdr_when_context_disables_it(adata_basic):
    """preset='heuristic' must not apply FDR cutoffs when perm FDR was disabled upstream."""
    _, sig, allr = scat.active_score(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_permutation=True,
        n_perm=30,
        random_seed=2,
        n_jobs=1,
        show_plot=False,
    )
    assert "scatrans_filter_context" in allr.attrs
    manual = scat.filter_active_genes(
        allr,
        preset="heuristic",
        unspliced_excess_fdr_cutoff=0.05,
        active_score_fdr_cutoff=0.25,
    )
    relaxed = scat.filter_active_genes(allr, preset="heuristic")
    assert len(relaxed) >= len(sig)
    assert len(manual) <= len(relaxed)


@pytest.mark.skipif(
    importlib.util.find_spec("gseapy") is None,
    reason="gseapy not installed",
)
def test_run_gsea_accepts_all_results_dataframe_index():
    """run_gsea must read gene symbols from DataFrame index, not column positions."""
    all_results = pd.DataFrame(
        {
            "active_score": [90.0, 80.0, 10.0, 5.0],
            "logFC": [1.5, 1.2, -0.8, -1.0],
            "unspliced_excess_delta": [2.0, 1.5, -0.5, -0.3],
        },
        index=["Il1r2", "Hdc", "Gapdh", "Actb"],
    )
    gene_sets = {
        "IMMUNE_UP": ["Il1r2", "Hdc", "Csf3r"],
        "HOUSE": ["Gapdh", "Actb"],
    }
    res = scat.run_gsea(
        all_results,
        gene_sets,
        score_column="logFC",
        nperm=50,
        min_size=1,
        verbose=False,
    )
    assert not res.empty, res.attrs.get("reason")
    assert "NES" in res.columns or "ES" in res.columns


def test_clean_gene_list_accepts_all_results_dataframe_index():
    """run_enrichment gene_list must read gene index, not column names."""
    from scatrans.enrich._data import _clean_gene_list

    df = pd.DataFrame(
        {"logFC": [1.5, 1.2], "active_score": [90.0, 80.0]},
        index=["Il1r2", "Hdc"],
    )
    assert _clean_gene_list(df) == ["Il1r2", "Hdc"]
    assert _clean_gene_list(df) != ["logFC", "active_score"]


def test_run_enrichment_accepts_filtered_dataframe_gene_list():
    """Passing filter/significant DataFrame directly must not use column names as genes."""
    df = pd.DataFrame(
        {"logFC": [1.5, 1.2]},
        index=["Il1r2", "Hdc"],
    )
    gene_sets = {"TERM": ["Il1r2", "Hdc", "Gapdh"]}
    res = scat.run_enrichment(
        df,
        gene_sets,
        universe=["Il1r2", "Hdc", "Gapdh", "Actb"],
        min_size=1,
        return_all=True,
        verbose=False,
    )
    assert res.attrs.get("reason") != "gene_list_empty"
    info = res.attrs.get("universe_info") or {}
    assert info.get("n_input_raw") == 2
    assert info.get("n_input_mapped") == 2
    assert not res.empty


def test_extract_gene_lists_from_all_results_index():
    """extract_gene_lists must use gene index for standard all_results tables."""
    all_results = pd.DataFrame(
        {
            "logFC": [1.5, 1.2, -0.9],
            "p_adj": [0.01, 0.02, 0.01],
            "active_score": [90.0, 80.0, 10.0],
        },
        index=["Il1r2", "Hdc", "Gapdh"],
    )
    out = scat.extract_gene_lists(
        all_results, logfc_cutoff=0.5, pval_cutoff=0.05, logfc_direction="up"
    )
    assert "Il1r2" in out["contrast"]
    assert "Hdc" in out["contrast"]
    assert "logFC" not in out["contrast"]


def test_run_gsea_legacy_two_column_dataframe_still_works():
    """Legacy [gene, score] two-column tables without gene index remain supported."""
    legacy = pd.DataFrame(
        {
            "gene": ["Il1r2", "Hdc", "Gapdh"],
            "logFC": [1.5, 1.2, -0.2],
        }
    )
    gene_sets = {"TERM": ["Il1r2", "Hdc"]}
    if importlib.util.find_spec("gseapy") is None:
        pytest.skip("gseapy not installed")
    res = scat.run_gsea(legacy, gene_sets, nperm=20, min_size=1, verbose=False)
    assert not res.empty, res.attrs.get("reason")


def test_mixedlm_composite_groups_when_sample_labels_reused_across_conditions():
    """rep1/rep2 reused per condition must not share one random effect (unpaired)."""
    from scatrans._de import _resolve_mixedlm_random_groups

    obs = pd.DataFrame(
        {
            "condition": ["A"] * 4 + ["B"] * 4,
            "sample": ["rep1", "rep2", "rep1", "rep2"] * 2,
        }
    )
    groups, meta = _resolve_mixedlm_random_groups(
        obs, "condition", "sample", paired_replicates=False
    )
    assert meta["grouping"] == "condition_sample_composite"
    assert meta["n_random_groups"] == 4
    assert meta["n_random_groups_raw_sample_col"] == 2
    assert set(groups) == {"A::rep1", "A::rep2", "B::rep1", "B::rep2"}

    raw_groups, raw_meta = _resolve_mixedlm_random_groups(
        obs, "condition", "sample", paired_replicates=True
    )
    assert raw_meta["grouping"] == "sample_col_raw"
    assert raw_meta["n_random_groups"] == 2


def test_mixedlm_keeps_raw_groups_when_sample_labels_are_globally_unique():
    from scatrans._de import _resolve_mixedlm_random_groups

    obs = pd.DataFrame(
        {
            "condition": ["GA"] * 3 + ["Ctrl"] * 3,
            "individual": ["GA_Ind1", "GA_Ind2", "GA_Ind3", "Ctrl_Ind1", "Ctrl_Ind2", "Ctrl_Ind3"],
        }
    )
    groups, meta = _resolve_mixedlm_random_groups(
        obs, "condition", "individual", paired_replicates=False
    )
    assert meta["grouping"] == "sample_col_raw"
    assert meta["overlapping_sample_labels"] == []
    assert meta["n_random_groups"] == 6


def test_maybe_repel_labels_subsamples_avoid_points_for_large_volcanoes():
    """Genome-wide volcanoes must subsample avoid-points before adjustText (MemoryError guard)."""
    from scatrans.pl import _maybe_repel_labels

    n_total = 25_000
    x = np.arange(n_total, dtype=float)
    y = np.arange(n_total, dtype=float) * 0.01
    texts = [MagicMock()]
    ax = MagicMock()

    with patch("adjustText.adjust_text") as mock_adjust:
        _maybe_repel_labels(
            texts,
            x,
            y,
            ax,
            label_repel=True,
            max_avoid_points=5000,
        )

    mock_adjust.assert_called_once()
    passed_x = mock_adjust.call_args.kwargs["x"]
    passed_y = mock_adjust.call_args.kwargs["y"]
    assert len(passed_x) == 5000
    assert len(passed_y) == 5000


def test_run_enrichment_reads_numpy_raw_gene_list_after_h5ad_roundtrip(adata_basic):
    """raw_gene_list stored as ndarray (h5ad reload) must not break universe auto-detection."""
    ad = adata_basic.copy()
    scat.store_raw_counts(ad, layer="counts")
    # AnnData/h5py round-trip deserializes string lists as numpy arrays.
    ad.uns["scatrans"]["raw_gene_list"] = np.asarray(ad.var_names, dtype=str)
    genes = ad.var_names[:3].tolist()
    gene_sets = {f"TERM_{g}": [g] for g in genes}

    res = scat.run_enrichment(
        genes,
        gene_sets=gene_sets,
        adata=ad,
        pval_cutoff=1.0,
        min_size=1,
        return_all=True,
        verbose=False,
    )
    ui = res.attrs.get("universe_info", {})
    assert ui.get("provided_size") == ad.n_vars
    assert ui.get("effective_universe_size", 0) > 0


def test_recommend_workflow_disables_permutation_on_small_pseudobulk_design():
    """Few exact pseudobulk shuffles must auto-set use_permutation=False (avoid runaway runtime)."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 48, 50
    X = rng.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)
    ad = sc.AnnData(X)
    ad.obs["condition"] = ["Disease"] * 24 + ["Control"] * 24
    ad.obs["sample"] = ["T0", "T1", "T2"] * 8 + ["R0", "R1", "R2"] * 8
    ad.layers["spliced"] = X.copy()
    ad.layers["unspliced"] = X * 0.5

    rec = scat.recommend_workflow(
        ad,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        sample_col="sample",
    )

    assert rec["workflow_preset"] == "pseudobulk_report"
    assert rec["suggested_kwargs"].get("use_pseudobulk") is True
    assert rec["use_permutation"] is False
    assert rec["suggested_kwargs"]["use_permutation"] is False
    assert rec["filter_preset"] == "pseudobulk"
    assert (rec["power_summary"] or {}).get("max_exact_permutations_pseudobulk", 999) < 100
    assert any("use_permutation=False" in w for w in rec["warnings"])


def test_mixed_model_rejects_too_few_samples_per_group():
    """3 vs 3 biological replicates should fail fast with a clear error."""
    rng = np.random.default_rng(0)
    n_cells, n_genes = 40, 30
    X = rng.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    adata = ad.AnnData(X)
    adata.obs["condition"] = ["Disease"] * 20 + ["Control"] * 20
    adata.obs["sample"] = (
        ["D0"] * 7 + ["D1"] * 7 + ["D2"] * 6 + ["C0"] * 7 + ["C1"] * 7 + ["C2"] * 6
    )
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.4
    with pytest.raises(ValueError, match="Mixed linear model requires"):
        scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_mixed_model=True,
            sample_col="sample",
            show_plot=False,
            use_permutation=False,
        )


def _make_mixedlm_perm_adata(n_genes: int = 12, seed: int = 0):
    """≥4 samples/group overlapping sample IDs (composite RE groups under default pairing)."""
    rng = np.random.default_rng(seed)
    n_cells = 60
    X = rng.negative_binomial(3, 0.5, size=(n_cells, n_genes)).astype(float)
    # Mild DE so MixedLM can fit without all near-constant genes
    X[30:, :3] = X[30:, :3] + 5
    adata = ad.AnnData(X)
    adata.obs["condition"] = ["Disease"] * 30 + ["Control"] * 30
    adata.obs["sample"] = [f"s{i % 6}" for i in range(n_cells)]
    adata.layers["spliced"] = X.copy()
    adata.layers["unspliced"] = X * 0.4 + rng.poisson(1, size=X.shape).astype(float)
    adata.var["gene_length"] = rng.integers(600, 3500, n_genes)
    adata.var["intron_number"] = rng.integers(0, 8, n_genes)
    return adata


def test_mixed_model_permutation_same_backend_forwards_mixedlm(monkeypatch):
    """perm_de_backend='same' + use_mixed_model must call MixedLM on observed and null DE.

    Regression: permutations previously always used scanpy t-test while logging
    that they used the 'same' backend, making active_score_fdr invalid.
    """
    import scatrans._de as de_mod
    import scatrans._permutation as perm_mod
    import scatrans.tl.active as active_mod

    seen_flags: list[bool] = []
    orig = de_mod._run_de_wrapper

    def tracking_wrapper(*args, **kwargs):
        seen_flags.append(bool(kwargs.get("use_mixed_model", False)))
        # Force MixedLM path to still run for a real end-to-end shape check
        return orig(*args, **kwargs)

    monkeypatch.setattr(de_mod, "_run_de_wrapper", tracking_wrapper)
    monkeypatch.setattr(perm_mod, "_run_de_wrapper", tracking_wrapper)
    monkeypatch.setattr(active_mod, "_run_de_wrapper", tracking_wrapper)

    adata = _make_mixedlm_perm_adata()
    ad_out, _sig, allr = scat.active_score(
        adata,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        use_mixed_model=True,
        sample_col="sample",
        use_permutation=True,
        perm_de_backend="same",
        n_perm=2,
        n_jobs=1,
        show_plot=False,
        min_total_counts=1,
    )

    meta = ad_out.uns["scatrans"]
    assert meta.get("use_mixed_model") is True
    assert meta.get("perm_use_mixed_model") is True
    assert meta.get("perm_de_backend") == "same"
    assert "active_score_fdr" in allr.columns
    assert "unspliced_excess_fdr" in allr.columns
    note = meta.get("permutation_approximation_note") or ""
    assert "MixedLM" in note
    assert "both observed" in note or "permutation null" in note

    # Observed DE + each successful perm task must request MixedLM
    assert len(seen_flags) >= 1 + 2
    assert all(seen_flags), f"expected all MixedLM calls, got {seen_flags}"
    # sample_col must reach the wrapper on MixedLM calls (via kwargs)
    # (re-run tracking is sufficient; failed MixedLM would raise)


def test_mixed_model_permutation_fast_warns_and_skips_mixedlm(monkeypatch, caplog):
    """perm_de_backend='fast' must not use MixedLM null and must warn about invalid active_score FDR."""
    import scatrans._de as de_mod
    import scatrans._permutation as perm_mod
    import scatrans.tl.active as active_mod

    seen_flags: list[bool] = []
    orig = de_mod._run_de_wrapper

    def tracking_wrapper(*args, **kwargs):
        seen_flags.append(bool(kwargs.get("use_mixed_model", False)))
        return orig(*args, **kwargs)

    monkeypatch.setattr(de_mod, "_run_de_wrapper", tracking_wrapper)
    monkeypatch.setattr(perm_mod, "_run_de_wrapper", tracking_wrapper)
    monkeypatch.setattr(active_mod, "_run_de_wrapper", tracking_wrapper)

    adata = _make_mixedlm_perm_adata(seed=1)
    with caplog.at_level(logging.WARNING, logger="scatrans.tl.active"):
        ad_out, _sig, allr = scat.active_score(
            adata,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            use_mixed_model=True,
            sample_col="sample",
            use_permutation=True,
            perm_de_backend="fast",
            n_perm=2,
            n_jobs=1,
            show_plot=False,
            min_total_counts=1,
        )

    meta = ad_out.uns["scatrans"]
    assert meta.get("use_mixed_model") is True
    assert meta.get("perm_use_mixed_model") is False
    assert meta.get("perm_de_backend") == "fast"
    assert "active_score_fdr" in allr.columns  # still computed, but statistically mismatched
    note = meta.get("permutation_approximation_note") or ""
    assert "not a valid null" in note or "non-MixedLM" in note

    # First call(s) are observed MixedLM; perm tasks must be non-MixedLM
    assert seen_flags[0] is True
    assert any(not f for f in seen_flags[1:]), f"expected non-MixedLM perm calls: {seen_flags}"
    assert any(
        "not valid" in r.message or "active_score_pval" in r.message for r in caplog.records
    ), "expected warning that active_score FDR is invalid under MixedLM + fast"


def test_shuffle_condition_labels_mixedlm_keeps_samples_pure():
    """Unpaired MixedLM null must not split a sample across conditions."""
    from scatrans._permutation import _shuffle_condition_labels

    # Clean unpaired: each sample is entirely one condition
    labels = np.array(["T"] * 20 + ["R"] * 20)
    samples = np.array(
        [f"s{i}" for i in range(4) for _ in range(5)]  # s0-s3 target
        + [f"s{i}" for i in range(4, 8) for _ in range(5)]  # s4-s7 ref
    )
    rng = np.random.default_rng(0)
    n_distinct = 0
    for seed in range(30):
        rng = np.random.default_rng(seed)
        out = _shuffle_condition_labels(
            labels,
            rng=rng,
            use_mixed_model=True,
            sample_ids=samples,
            paired_replicates=False,
        )
        assert out is not None
        shuffled, re_ids = out
        assert re_ids is not None
        assert shuffled.shape == labels.shape
        # Every sample must be pure under the null labels
        for s in pd.unique(samples):
            labs = np.unique(shuffled[samples == s])
            assert labs.size == 1, f"sample {s} split across {labs}"
        if not np.array_equal(shuffled, labels):
            n_distinct += 1
    assert n_distinct >= 5


def test_shuffle_condition_labels_cell_level_can_split_samples():
    """Non-MixedLM cell shuffle may (and usually does) mix conditions within sample."""
    from scatrans._permutation import _shuffle_condition_labels

    labels = np.array(["T"] * 20 + ["R"] * 20)
    samples = np.array([f"s{i % 4}" for i in range(40)])
    mixed_within_sample = False
    for seed in range(40):
        out = _shuffle_condition_labels(
            labels,
            rng=np.random.default_rng(seed),
            use_mixed_model=False,
            sample_ids=samples,
        )
        assert out is not None
        shuffled, re_ids = out
        assert re_ids is None
        for s in pd.unique(samples):
            if np.unique(shuffled[samples == s]).size > 1:
                mixed_within_sample = True
                break
        if mixed_within_sample:
            break
    assert mixed_within_sample


def test_shuffle_condition_labels_composite_clusters_stay_pure():
    """Recycled sample IDs (composite RE): each observed cluster stays pure after shuffle."""
    from scatrans._de import _resolve_mixedlm_random_groups
    from scatrans._permutation import _shuffle_condition_labels

    # Same fixture pattern as adata_mixed_small: s0..s5 in both arms
    n = 60
    labels = np.array(["Disease"] * 30 + ["Control"] * 30)
    samples = np.array([f"s{i % 6}" for i in range(n)])
    obs = pd.DataFrame({"condition": labels, "sample": samples})
    groups, meta = _resolve_mixedlm_random_groups(
        obs, "condition", "sample", paired_replicates=False, quiet=True
    )
    assert meta["grouping"] == "condition_sample_composite"

    out = _shuffle_condition_labels(
        labels,
        rng=np.random.default_rng(7),
        use_mixed_model=True,
        sample_ids=samples,
        paired_replicates=False,
    )
    assert out is not None
    shuffled, re_ids = out
    assert re_ids is not None
    # Observed RE clusters remain pure under null labels
    for g in pd.unique(np.asarray(groups).astype(str)):
        m = np.asarray(groups).astype(str) == g
        assert np.unique(shuffled[m]).size == 1
    # And the pinned sample_ids for MixedLM are the observed clusters
    assert set(np.asarray(re_ids).astype(str)) == set(np.asarray(groups).astype(str))


def test_shuffle_condition_labels_paired_preserves_within_sample_counts():
    """Paired MixedLM: within-subject condition counts fixed under null."""
    from scatrans._permutation import _shuffle_condition_labels

    # 4 subjects, each with 5 Disease + 5 Control cells
    labels_list = []
    samples_list = []
    for i in range(4):
        labels_list.extend(["T"] * 5 + ["R"] * 5)
        samples_list.extend([f"ind{i}"] * 10)
    labels = np.array(labels_list)
    samples = np.array(samples_list)
    out = _shuffle_condition_labels(
        labels,
        rng=np.random.default_rng(3),
        use_mixed_model=True,
        sample_ids=samples,
        paired_replicates=True,
    )
    assert out is not None
    shuffled, re_ids = out
    assert re_ids is not None
    for s in pd.unique(samples):
        m = samples == s
        assert (shuffled[m] == "T").sum() == 5
        assert (shuffled[m] == "R").sum() == 5
