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


def _flag_table(seed: int = 0) -> pd.DataFrame:
    """Large, all-significant-induced table so the induction-confound flag engages
    (>= _MIN_INDUCED_FOR_FLAG). Mostly stabilization-driven (negative residual)."""
    rng = np.random.default_rng(seed)
    n = 400
    df = pd.DataFrame(
        {
            "logFC": rng.uniform(1.0, 6.0, n),
            "p_adj": np.full(n, 1e-4),
            RC: rng.normal(-2.0, 0.5, n),
        },
        index=[f"g{i}" for i in range(n)],
    )
    df.loc["hi_stab", ["logFC", "p_adj", RC]] = [10.0, 1e-6, -3.0]  # extreme induction -> FLAG
    df.loc["lo_stab", ["logFC", "p_adj", RC]] = [1.1, 1e-6, -3.0]  # mild induction -> not flagged
    df.loc["hi_txn", ["logFC", "p_adj", RC]] = [10.0, 1e-6, 3.0]  # extreme but transcription
    return df


# ---------------- induction-confound flag (6) ----------------
def test_induction_flag_default_flags_extreme_stabilization():
    from scatrans.tl.mechanism import INDUCTION_CONFOUND_COL

    out, diag = scat.annotate_mechanism_class(_flag_table())
    assert INDUCTION_CONFOUND_COL in out.columns
    assert diag["n_induction_confounded"] > 0
    # extreme-induction stabilization call is flagged; mild one and the txn call are not
    assert bool(out.loc["hi_stab", INDUCTION_CONFOUND_COL]) is True
    assert bool(out.loc["lo_stab", INDUCTION_CONFOUND_COL]) is False
    assert bool(out.loc["hi_txn", INDUCTION_CONFOUND_COL]) is False


def test_induction_flag_discounts_only_flagged_confidence():
    from scatrans.tl.mechanism import INDUCTION_CONFOUND_COL

    on, _ = scat.annotate_mechanism_class(_flag_table())
    off, _ = scat.annotate_mechanism_class(_flag_table(), flag_induction_confound=False)
    assert not off[INDUCTION_CONFOUND_COL].any()
    # flagged gene: confidence strictly reduced; unflagged mild stab: unchanged
    assert on.loc["hi_stab", CONF_COL] < off.loc["hi_stab", CONF_COL]
    assert np.isclose(on.loc["lo_stab", CONF_COL], off.loc["lo_stab", CONF_COL])
    # graded floor: never below 0.3x the base confidence
    assert on.loc["hi_stab", CONF_COL] >= 0.3 * off.loc["hi_stab", CONF_COL] - 1e-9


def test_induction_flag_is_program_invariant():
    """The flag must not touch transcription_support (program-level calls stay put)."""
    on, _ = scat.annotate_mechanism_class(_flag_table())
    off, _ = scat.annotate_mechanism_class(_flag_table(), flag_induction_confound=False)
    assert np.allclose(on[SUPPORT_COL], off[SUPPORT_COL], equal_nan=True)
    assert on[CLASS_COL].equals(off[CLASS_COL])


def test_induction_flag_smooth_penalty():
    on, diag = scat.annotate_mechanism_class(_flag_table(), induction_confound_penalty="smooth")
    off, _ = scat.annotate_mechanism_class(_flag_table(), flag_induction_confound=False)
    assert diag["induction_confound_penalty"] == "smooth"
    assert on.loc["hi_stab", CONF_COL] < off.loc["hi_stab", CONF_COL]


def test_induction_flag_skipped_on_small_table():
    # the default fixture has too few significant induced genes -> flag is a no-op
    _, diag = scat.annotate_mechanism_class(_table())
    assert diag["n_induction_confounded"] == 0


def test_induction_flag_validation():
    with pytest.raises(ValueError):
        scat.annotate_mechanism_class(_flag_table(), induction_confound_penalty="nope")
    with pytest.raises(ValueError):
        scat.annotate_mechanism_class(_flag_table(), induction_confound_quantile=1.5)


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
    with pytest.raises(ValueError, match="preset"):
        scat.annotate_mechanism_class(_table(), preset="nope")


def test_high_precision_preset_sets_threshold():
    _, diag = scat.annotate_mechanism_class(_table(), preset="high_precision")
    assert diag["class_threshold"] == 1.0
    assert diag["preset"] == "high_precision"


def test_high_precision_preset_fewer_positive_calls():
    # Raising the threshold 0.5 -> 1.0 can only shrink the set of hard calls
    # (borderline genes fall back to "ambiguous"); it never creates new ones.
    base, _ = scat.annotate_mechanism_class(_table())
    hp, _ = scat.annotate_mechanism_class(_table(), preset="high_precision")
    called = {"transcription-driven", "stabilization-driven"}
    n_base = base[CLASS_COL].isin(called).sum()
    n_hp = hp[CLASS_COL].isin(called).sum()
    assert n_hp <= n_base
    # every high-precision call is also a base call (threshold only tightened)
    hp_called = hp.index[hp[CLASS_COL].isin(called)]
    assert (base.loc[hp_called, CLASS_COL] == hp.loc[hp_called, CLASS_COL]).all()


def test_explicit_args_override_preset():
    # An explicit class_threshold beats the preset's 1.0.
    _, diag = scat.annotate_mechanism_class(_table(), preset="high_precision", class_threshold=0.3)
    assert diag["class_threshold"] == 0.3


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


def test_hard_labels_suppressed_when_reliability_near_zero():
    out, diag = scat.annotate_mechanism_class(_table(), reliability=0.0)
    hard = {"transcription-driven", "stabilization-driven"}
    assert not out[CLASS_COL].isin(hard).any()
    assert diag["hard_labels_suppressed"] is True
    assert diag["n_hard_labels_suppressed"] >= 1
    # support still computed; confidence scaled to ~0
    assert SUPPORT_COL in out.columns
    up = out[CONF_COL].notna()
    assert np.allclose(out.loc[up, CONF_COL].fillna(0), 0.0)


def test_hard_labels_kept_when_suppression_disabled():
    out, diag = scat.annotate_mechanism_class(
        _table(), reliability=0.0, suppress_hard_labels_when_unreliable=False
    )
    assert diag["hard_labels_suppressed"] is False
    assert out.loc["g0", CLASS_COL] == "transcription-driven"
    assert out.loc["g1", CLASS_COL] == "stabilization-driven"


def test_hard_labels_not_suppressed_at_full_reliability():
    out, diag = scat.annotate_mechanism_class(_table(), reliability=1.0)
    assert diag["hard_labels_suppressed"] is False
    assert out.loc["g0", CLASS_COL] == "transcription-driven"


# ---------------- program_mechanism_induction_matched ----------------
def _induction_matched_table(seed: int = 7):
    """Induced genes: program genes have lower support at matched logFC."""
    rng = np.random.default_rng(seed)
    n_bg, n_prog = 200, 25
    bg_idx = [f"bg{i}" for i in range(n_bg)]
    pg_idx = [f"pg{i}" for i in range(n_prog)]
    # logFC ~ U(1, 4) for both; program residual lower by ~2
    lfc_bg = rng.uniform(1.0, 4.0, n_bg)
    lfc_pg = rng.uniform(1.0, 4.0, n_prog)
    # residual anticorrelated with logFC mildly + program offset
    rc_bg = rng.normal(0.5, 0.4, n_bg) - 0.1 * (lfc_bg - 2.5)
    rc_pg = rng.normal(-1.5, 0.4, n_prog) - 0.1 * (lfc_pg - 2.5)
    df = pd.DataFrame(
        {
            "logFC": np.concatenate([lfc_bg, lfc_pg]),
            "p_adj": 1e-4,
            RC: np.concatenate([rc_bg, rc_pg]),
        },
        index=bg_idx + pg_idx,
    )
    return df, pg_idx


def test_program_mechanism_induction_matched_detects_stabilization():
    df, prog = _induction_matched_table()
    # annotate support first (function can also derive from residual)
    df, _ = scat.annotate_mechanism_class(df, flag_induction_confound=False)
    res = scat.program_mechanism_induction_matched(
        df, {"STAB_PROG": prog}, min_genes=10, methods=("regression", "nearest")
    )
    assert len(res) == 1
    row = res.iloc[0]
    assert row["program"] == "STAB_PROG"
    assert row["regression_beta"] < 0
    assert row["regression_p"] < 0.05
    assert row["direction"] == "stabilization-driven"
    assert row["nearest_median_delta"] < 0


def test_program_mechanism_induction_matched_null_program_ns():
    df, _ = _induction_matched_table()
    df, _ = scat.annotate_mechanism_class(df, flag_induction_confound=False)
    # random background genes as "program" — should be ns
    null = [f"bg{i}" for i in range(25)]
    res = scat.program_mechanism_induction_matched(
        df, {"NULL": null}, min_genes=10, methods=("regression",)
    )
    assert len(res) == 1
    assert res.iloc[0]["direction"] == "ns" or res.iloc[0]["regression_p"] > 0.01


def test_program_mechanism_induction_matched_methods_validation():
    df, prog = _induction_matched_table()
    with pytest.raises(ValueError, match="methods"):
        scat.program_mechanism_induction_matched(df, {"P": prog}, methods=("nope",))


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
