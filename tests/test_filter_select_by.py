"""Tests for filter_active_genes(select_by="de") — "DE SELECTS, proxy ANNOTATES".

select_by="de" must decide gene-list membership from the DE gates only
(p_adj / logFC + direction), skipping the nascent-proxy gates (active_score,
unspliced_excess_residual, the FDR columns, effective_gamma, delta_variance),
while keeping those proxy columns on the output as annotations. Default
select_by="composite" must reproduce the pre-existing behavior exactly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scatrans._utils import UNSPLICED_EXCESS_FDR_COL, UNSPLICED_EXCESS_RESIDUAL_COL
from scatrans.tl.filter import filter_active_genes


def _make_table(seed: int = 0, n: int = 200) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(
        {
            "logFC": rng.normal(0, 2, n),
            "p_adj": rng.uniform(0, 1, n),
            "active_score": rng.uniform(0, 100, n),
            UNSPLICED_EXCESS_RESIDUAL_COL: rng.normal(0, 1, n),
            UNSPLICED_EXCESS_FDR_COL: rng.uniform(0, 1, n),
        },
        index=[f"g{i}" for i in range(n)],
    )
    # g0: strong DE-up but BAD proxy (low residual, high FDR, low active_score)
    df.loc["g0", ["logFC", "p_adj"]] = [3.0, 1e-6]
    df.loc["g0", [UNSPLICED_EXCESS_RESIDUAL_COL, UNSPLICED_EXCESS_FDR_COL, "active_score"]] = [
        -2.0,
        0.99,
        1.0,
    ]
    # g1: strong DE-up but proxy columns are NaN (not computable)
    df.loc["g1", ["logFC", "p_adj"]] = [2.5, 1e-5]
    df.loc["g1", [UNSPLICED_EXCESS_RESIDUAL_COL, UNSPLICED_EXCESS_FDR_COL, "active_score"]] = [
        np.nan,
        np.nan,
        np.nan,
    ]
    return df


def test_de_selects_gene_with_bad_proxy():
    df = _make_table()
    de = filter_active_genes(df, select_by="de")
    comp = filter_active_genes(df, preset="pseudobulk")
    # DE-up gene with a bad proxy is selected on the DE axis but excluded by composite
    assert "g0" in de.index
    assert "g0" not in comp.index


def test_de_keeps_gene_with_nan_proxy():
    df = _make_table()
    de = filter_active_genes(df, select_by="de")
    # a strong-DE gene whose proxy columns are NaN must remain selectable
    assert "g1" in de.index
    assert not de.loc["g1", ["logFC", "p_adj"]].isna().any()


def test_de_gates_are_the_agreed_defaults():
    df = _make_table()
    de = filter_active_genes(df, select_by="de")
    # default DE standard: padj < 0.05 AND log2FC > 1 (up)
    assert (de["p_adj"] < 0.05).all()
    assert (de["logFC"] > 1).all()


def test_de_retains_proxy_columns_as_annotation():
    df = _make_table()
    de = filter_active_genes(df, select_by="de")
    for col in ("active_score", UNSPLICED_EXCESS_RESIDUAL_COL, UNSPLICED_EXCESS_FDR_COL):
        assert col in de.columns


def test_de_sorted_by_padj_not_active_score():
    df = _make_table()
    de = filter_active_genes(df, select_by="de")
    padj = de["p_adj"].to_numpy()
    assert np.all(np.diff(padj) >= 0)  # ascending p_adj (DE ranking), not active_score


def test_composite_default_unchanged():
    df = _make_table()
    a = filter_active_genes(df, preset="pseudobulk")
    b = filter_active_genes(df, preset="pseudobulk", select_by="composite")
    assert a.index.equals(b.index)


def test_de_direction_down():
    df = _make_table()
    df.loc["g2", ["logFC", "p_adj"]] = [-3.0, 1e-6]  # strong down gene
    de_down = filter_active_genes(df, select_by="de", logfc_direction="down")
    assert "g2" in de_down.index
    assert (de_down["logFC"] < -1).all()
    assert (de_down["p_adj"] < 0.05).all()


def test_de_incompatible_with_significant_preset():
    df = _make_table()
    with pytest.raises(ValueError, match="incompatible with preset='significant'"):
        filter_active_genes(df, preset="significant", select_by="de")


def test_invalid_select_by_raises():
    df = _make_table()
    with pytest.raises(ValueError, match="select_by must be"):
        filter_active_genes(df, select_by="nonsense")


def test_de_explicit_cutoffs_override_defaults():
    df = _make_table()
    # relax to keep more genes; explicit args beat the 0.05/1.0 defaults
    lenient = filter_active_genes(df, select_by="de", padj_cutoff=0.5, logfc_cutoff=0.0)
    strict = filter_active_genes(df, select_by="de")
    assert len(lenient) > len(strict)


def test_pipeline_select_by_de_end_to_end(adata_basic):
    import scatrans as scat

    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
        select_by="de",
    )
    assert res.meta["select_by"] == "de"
    cand = res.candidates
    # candidates chosen by DE gates; proxy columns still present as annotations
    if len(cand) > 0:
        assert (cand["p_adj"] < 0.05).all()
        assert "unspliced_excess_residual" in cand.columns
    # composite default still works and records its mode
    res_c = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
    )
    assert res_c.meta["select_by"] == "composite"
