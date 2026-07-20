"""Tests for tl.mechanism: annotate_mechanism_class, threshold_sensitivity,
program_mechanism (transcription-vs-stabilization annotation — never gates)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import scatrans as scat
from scatrans._utils import UNSPLICED_EXCESS_RESIDUAL_COL as RC
from scatrans.tl.bias import RESID_COL as ABNORM
from scatrans.tl.mechanism import CLASS_COL, CONF_COL, SUPPORT_COL


def _table(seed: int = 1, n: int = 300) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(
        {
            "logFC": rng.normal(0.5, 1.5, n),
            "p_adj": rng.uniform(0, 1, n),
            RC: rng.normal(0, 1, n),
        },
        index=[f"g{i}" for i in range(n)],
    )
    df.loc["g0", ["logFC", "p_adj", RC]] = [2.0, 1e-4, 3.0]  # up + high excess
    df.loc["g1", ["logFC", "p_adj", RC]] = [2.0, 1e-4, -3.0]  # up + low excess
    df.loc["g2", ["logFC", "p_adj", RC]] = [-2.0, 1e-4, 1.0]  # down
    return df


# ---------------- annotate_mechanism_class (2) ----------------
def test_mechanism_class_labels():
    out, diag = scat.annotate_mechanism_class(_table())
    assert out.loc["g0", CLASS_COL] == "transcription-driven"
    assert out.loc["g1", CLASS_COL] == "stabilization-driven"
    assert out.loc["g2", CLASS_COL] == "unclassified_down"
    assert diag["n_classified"] > 0


def test_mechanism_adds_columns_and_preserves_membership():
    df = _table()
    out, _ = scat.annotate_mechanism_class(df)
    for c in (SUPPORT_COL, CLASS_COL, CONF_COL):
        assert c in out.columns
    # annotation only: same rows, same order
    assert out.index.equals(df.index)


def test_mechanism_reliability_scales_confidence():
    out0, _ = scat.annotate_mechanism_class(_table(), reliability=0.0)
    out1, _ = scat.annotate_mechanism_class(_table(), reliability=1.0)
    up = out1[CONF_COL].notna()
    assert np.allclose(out0.loc[up, CONF_COL].fillna(0), 0.0)
    assert (out1.loc[up, CONF_COL] >= out0.loc[up, CONF_COL].fillna(0)).all()


def test_mechanism_prefers_bias_corrected_residual():
    df = _table()
    df[ABNORM] = -df[RC]  # opposite sign; if preferred, g0 flips to stabilization
    out, diag = scat.annotate_mechanism_class(df)
    assert diag["residual_col"] == ABNORM
    assert out.loc["g0", CLASS_COL] == "stabilization-driven"


def test_mechanism_missing_residual_raises():
    df = _table().drop(columns=[RC])
    with pytest.raises(KeyError, match="residual"):
        scat.annotate_mechanism_class(df)


def test_mechanism_invalid_params_raise():
    with pytest.raises(ValueError):
        scat.annotate_mechanism_class(_table(), reliability=2.0)
    with pytest.raises(ValueError):
        scat.annotate_mechanism_class(_table(), class_threshold=-1.0)


# ---------------- threshold_sensitivity (3) ----------------
def test_threshold_sensitivity_grid_and_reference():
    ts = scat.threshold_sensitivity(
        _table(),
        padj_grid=(0.01, 0.05, 0.1),
        logfc_grid=(0.58, 1.0, 1.5),
        reference=(0.05, 1.0),
    )
    assert len(ts) == 9
    assert set(ts.columns) == {
        "padj_cutoff",
        "logfc_cutoff",
        "n_selected",
        "jaccard_vs_reference",
        "is_reference",
    }
    ref = ts[ts["is_reference"]]
    assert len(ref) == 1
    assert np.isclose(ref["jaccard_vs_reference"].iloc[0], 1.0)


def test_threshold_sensitivity_looser_padj_selects_more():
    ts = scat.threshold_sensitivity(_table(), padj_grid=(0.01, 0.1), logfc_grid=(1.0,))
    n_strict = ts[np.isclose(ts.padj_cutoff, 0.01)]["n_selected"].iloc[0]
    n_loose = ts[np.isclose(ts.padj_cutoff, 0.1)]["n_selected"].iloc[0]
    assert n_loose >= n_strict


# ---------------- program_mechanism (5) ----------------
def _table_with_program(seed: int = 2):
    rng = np.random.default_rng(seed)
    df = _table(seed)
    trans = [f"t{i}" for i in range(30)]
    extra = pd.DataFrame({"logFC": 2.0, "p_adj": 1e-3, RC: rng.normal(1.5, 0.8, 30)}, index=trans)
    return pd.concat([df, extra]), trans


def test_program_mechanism_detects_transcription_program():
    df, trans = _table_with_program()
    gene_sets = {"TRANS": trans, "OTHER": [f"g{i}" for i in range(30)]}
    pm = scat.program_mechanism(df, gene_sets, min_genes=5)
    row = pm[pm["program"] == "TRANS"].iloc[0]
    assert row["direction"] == "transcription-driven"
    assert row["significant"]
    assert row["mean_support"] > row["bg_mean_support"]
    assert {"fdr", "p_value", "n_genes"} <= set(pm.columns)


def test_program_mechanism_min_genes_filters():
    df, trans = _table_with_program()
    gene_sets = {"TRANS": trans, "TINY": trans[:3]}
    pm = scat.program_mechanism(df, gene_sets, min_genes=5)
    assert "TINY" not in set(pm["program"])
    assert "TRANS" in set(pm["program"])


def test_program_mechanism_null_program_not_inflated_by_strong_one():
    # A strong transcription program must NOT make a genuinely-null program read
    # "stabilization-driven" via background inflation (exclude_other_programs=True).
    df, trans = _table_with_program()
    null_genes = [f"g{i}" for i in range(30)]  # ~mean-0 support, truly null
    gene_sets = {"TRANS": trans, "NULL": null_genes}
    pm = scat.program_mechanism(df, gene_sets, min_genes=5)  # default excludes
    assert pm[pm.program == "TRANS"].iloc[0]["direction"] == "transcription-driven"
    assert pm[pm.program == "NULL"].iloc[0]["direction"] == "ns"


def test_program_mechanism_sorted_by_pvalue():
    df, trans = _table_with_program()
    gene_sets = {
        "TRANS": trans,
        "A": [f"g{i}" for i in range(30)],
        "B": [f"g{i}" for i in range(30, 60)],
    }
    pm = scat.program_mechanism(df, gene_sets, min_genes=5)
    assert list(pm["p_value"]) == sorted(pm["p_value"])


# ---------------- pipeline integration ----------------
def test_pipeline_annotate_mechanism_flag(adata_basic):
    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
        annotate_mechanism=True,
    )
    assert CLASS_COL in res.all_results.columns
    assert SUPPORT_COL in res.all_results.columns
    assert "mechanism" in res.meta
    # default (flag off) does not add the column
    res_off = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        run_go_enrichment=False,
        show_plot=False,
    )
    assert CLASS_COL not in res_off.all_results.columns
