"""Regression tests for confirmed bug fixes (2026-07 audit)."""

from __future__ import annotations

import logging
import sys
from unittest.mock import patch

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
from scatrans.tl import _materialize_if_view, _select_var


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
