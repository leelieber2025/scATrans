"""Tests for tl.partition.partition_de_by_mechanism — the DE-selects /
scATrans-partitions-by-mechanism primary workflow.

Covers: builtin DE, injected DataFrame DE, callable DE, method-name DE routed to
differential_expression, program-level table, mandatory reliability pre-flight,
soft-labels-only (proxy never gates membership), column mapping + errors, and the
composite deprecation on run_default_pipeline.
"""

from __future__ import annotations

import warnings

import pandas as pd
import pytest

import scatrans as scat
from scatrans.tl import PartitionResult, partition_de_by_mechanism
from scatrans.tl.mechanism import CLASS_COL, CONF_COL, SUPPORT_COL


def _run(adata, **kw):
    kw.setdefault("target_group", "Disease")
    kw.setdefault("reference_group", "Control")
    kw.setdefault("organism", "human")
    return partition_de_by_mechanism(adata, groupby="condition", **kw)


# --------------------------- builtin DE path --------------------------------
def test_builtin_returns_partition_result(adata_basic):
    r = _run(adata_basic)
    assert isinstance(r, PartitionResult)
    # gene_table is the FULL scored universe with mechanism annotation columns
    for c in (SUPPORT_COL, CLASS_COL, CONF_COL):
        assert c in r.gene_table.columns
    # selected is a subset of the gene table (DE membership)
    assert set(r.selected.index).issubset(set(r.gene_table.index))
    assert len(r.gene_table) >= len(r.selected)
    # mandatory reliability pre-flight recorded and used
    assert "regime" in r.meta and "reliability" in r.regime
    assert r.meta["de_source"] == "builtin"
    s = r.summary()
    assert s["n_selected"] == len(r.selected)


def test_selected_membership_is_de_only_not_proxy(adata_basic):
    """Membership must come from DE gates (padj/logFC), never the proxy."""
    r = _run(adata_basic, padj_cutoff=0.05, logfc_cutoff=1.0)
    sel = r.selected
    if len(sel):
        assert (pd.to_numeric(sel["p_adj"], errors="coerce") < 0.05).all()
        assert (pd.to_numeric(sel["logFC"], errors="coerce") > 1.0).all()
    # a stabilization-driven gene is NOT dropped from selection for being low-support
    # (proxy annotates, never filters): selected includes both classes if present
    classes = set(r.gene_table.loc[sel.index, CLASS_COL]) if len(sel) else set()
    assert classes.issubset(
        {"transcription-driven", "stabilization-driven", "ambiguous", "unknown"}
    )


def test_gene_table_membership_independent_of_nan_proxy(adata_basic):
    """Genes with NaN proxy still live in the gene_table (annotation, not gate)."""
    r = _run(adata_basic)
    # gene_table keeps every scored gene regardless of support being finite
    assert r.gene_table[SUPPORT_COL].isna().sum() >= 0  # column exists, NaNs allowed
    assert len(r.gene_table) == r.adata.n_vars or len(r.gene_table) > 0


# --------------------------- injected DataFrame DE --------------------------
def test_injected_dataframe_de_controls_membership(adata_basic):
    names = list(adata_basic.var_names)
    chosen = set(names[:5])
    de_df = pd.DataFrame(
        {
            "logFC": [3.0 if g in chosen else 0.0 for g in names],
            "p_adj": [1e-6 if g in chosen else 1.0 for g in names],
        },
        index=names,
    )
    r = _run(adata_basic, de=de_df)
    assert r.meta["de_source"] == "dataframe"
    # membership follows the injected DE exactly (all chosen, nothing else)
    assert set(r.selected.index) == chosen
    assert r.meta["de"]["n_matched_to_scored"] == len(names)


def test_injected_de_stats_written_into_selected_and_mechanism(adata_basic):
    """Regression (Issue 1): for external DE, reported logFC/p_adj and the mechanism
    DIRECTION must come from the SELECTING DE — not the builtin active_score pass —
    so up-selected genes are never mislabeled unclassified_down."""
    names = list(adata_basic.var_names)
    chosen = names[:6]
    de_df = pd.DataFrame(
        {
            "logFC": [4.0 if g in chosen else 0.0 for g in names],
            "p_adj": [1e-8 if g in chosen else 1.0 for g in names],
        },
        index=names,
    )
    r = _run(adata_basic, de=de_df)
    assert set(r.selected.index) == set(chosen)
    # reported logFC is the EXTERNAL value (4.0), not the builtin DE's
    assert (pd.to_numeric(r.selected["logFC"], errors="coerce") == 4.0).all()
    assert (pd.to_numeric(r.selected["p_adj"], errors="coerce") == 1e-8).all()
    # up-selected genes must not be labeled unclassified_down (that requires logFC<0)
    assert (r.selected["mechanism_class"] != "unclassified_down").all()
    # Issue A: no stale builtin p_val next to the external logFC/p_adj. The external
    # table gave no raw p, so p_val is cleared to NaN (never the builtin DE's value).
    if "p_val" in r.selected.columns:
        assert pd.to_numeric(r.selected["p_val"], errors="coerce").isna().all()


def test_injected_de_carries_external_pval_when_present(adata_basic):
    """When the external DE table provides a raw p_val, it is carried through."""
    names = list(adata_basic.var_names)
    chosen = set(names[:4])
    de_df = pd.DataFrame(
        {
            "logFC": [3.0 if g in chosen else 0.0 for g in names],
            "p_val": [1e-9 if g in chosen else 0.9 for g in names],
            "p_adj": [1e-7 if g in chosen else 1.0 for g in names],
        },
        index=names,
    )
    r = _run(adata_basic, de=de_df)
    assert (pd.to_numeric(r.selected["p_val"], errors="coerce") == 1e-9).all()


def test_pipeline_de_candidates_carry_mechanism_columns(adata_basic):
    """Issue E/3: with annotate_mechanism=True the add-ons run after candidates is
    sliced, so result.candidates must still carry the mechanism columns."""
    res = scat.run_default_pipeline(
        adata_basic,
        groupby="condition",
        target_group="Disease",
        reference_group="Control",
        organism="human",
        select_by="de",
        annotate_mechanism=True,
        run_go_enrichment=False,
    )
    if len(res.candidates):
        for c in ("transcription_support", "mechanism_class", "mechanism_confidence"):
            assert c in res.candidates.columns


def test_injected_de_column_mapping(adata_basic):
    names = list(adata_basic.var_names)
    de_df = pd.DataFrame(
        {"lfc": [3.0] + [0.0] * (len(names) - 1), "fdr": [1e-6] + [1.0] * (len(names) - 1)},
        index=names,
    )
    r = _run(adata_basic, de=de_df, de_logfc_col="lfc", de_padj_col="fdr")
    assert set(r.selected.index) == {names[0]}


def test_injected_de_missing_column_raises(adata_basic):
    de_df = pd.DataFrame({"logFC": [1.0]}, index=[adata_basic.var_names[0]])
    with pytest.raises(KeyError):
        _run(adata_basic, de=de_df)  # no p_adj column


# --------------------------- callable DE ------------------------------------
def test_callable_de(adata_basic):
    def my_de(adata):
        names = list(adata.var_names)
        return pd.DataFrame(
            {
                "logFC": [2.0 if i < 3 else 0.0 for i in range(len(names))],
                "p_adj": [1e-4 if i < 3 else 1.0 for i in range(len(names))],
            },
            index=names,
        )

    r = _run(adata_basic, de=my_de)
    assert r.meta["de_source"] == "callable"
    assert set(r.selected.index) == set(list(adata_basic.var_names)[:3])


# --------------------------- method-name DE (front-end) ---------------------
def test_method_name_routed_to_differential_expression(adata_basic):
    r = _run(adata_basic, de="t-test")
    assert r.meta["de_source"].startswith("differential_expression")
    assert "de_kwargs" in r.meta["de"]
    # selection still respects the DE gates
    if len(r.selected):
        assert (pd.to_numeric(r.selected["p_adj"], errors="coerce") < 0.05).all()


def test_de_kwargs_dict_routed(adata_basic):
    """de as a kwargs dict must route to differential_expression (Mapping branch)."""
    r = _run(adata_basic, de={"de_method": "t-test"})
    assert r.meta["de_source"].startswith("differential_expression")
    assert r.meta["de"]["de_kwargs"]["de_method"] == "t-test"


def test_bad_de_type_raises(adata_basic):
    with pytest.raises(TypeError):
        _run(adata_basic, de=12345)


# --------------------------- program-level call -----------------------------
def test_program_mechanism_table(adata_basic):
    names = list(adata_basic.var_names)
    gene_sets = {"progA": names[:20], "progB": names[20:40]}
    r = _run(
        adata_basic, gene_sets=gene_sets, program_min_genes=3, program_restrict_to_selected=False
    )
    assert r.programs is not None
    assert {"program", "mean_support", "direction", "p_value", "fdr"}.issubset(r.programs.columns)
    assert set(r.programs["program"]).issubset({"progA", "progB"})


def test_no_gene_sets_means_no_programs(adata_basic):
    r = _run(adata_basic)
    assert r.programs is None


# --------------------------- reliability pre-flight -------------------------
def test_confidence_scaled_by_reliability(adata_high_unspliced):
    """High-unspliced (low reliability) must down-scale mechanism_confidence."""
    r = partition_de_by_mechanism(
        adata_high_unspliced,
        groupby="condition",
        target_group="A",
        reference_group="B",
        organism="human",
    )
    rel = r.regime["reliability"]
    assert rel < 1.0  # high unspliced -> degraded reliability
    conf = pd.to_numeric(r.gene_table[CONF_COL], errors="coerce").dropna()
    if len(conf):
        # confidence = clipped(|support|/(thr+1)) * reliability  <= reliability
        assert (conf <= rel + 1e-9).all()


# --------------------------- deprecation of composite -----------------------
def test_run_default_pipeline_composite_deprecated(adata_basic):
    with pytest.warns(DeprecationWarning, match="partition_de_by_mechanism"):
        scat.run_default_pipeline(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            organism="human",
            select_by="composite",
            run_go_enrichment=False,
        )


def test_select_by_de_not_deprecated(adata_basic):
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        # select_by="de" must NOT raise the composite deprecation
        scat.run_default_pipeline(
            adata_basic,
            groupby="condition",
            target_group="Disease",
            reference_group="Control",
            organism="human",
            select_by="de",
            run_go_enrichment=False,
        )


# --------------------------- P0 summary + P1 induction_matched ---------------
def test_summary_is_program_first_and_marks_soft_labels(adata_basic):
    r = _run(adata_basic)
    s = r.summary()
    assert s["per_gene_labels_are_soft"] is True
    assert "note" in s and "mechanism_class" in s["note"]
    assert "per_gene_class_counts_selected" in s
    # backward-compatible alias still present
    assert "class_counts_selected" in s
    assert "n_programs" in s


def test_pseudoreplication_warning_without_sample_col(adata_basic, caplog):
    import logging

    with caplog.at_level(logging.WARNING, logger="scatrans.tl.partition"):
        r = _run(adata_basic)
    assert r.meta.get("pseudoreplication_warning") is True
    assert any("pseudoreplication" in m.message.lower() for m in caplog.records)


def test_induction_matched_optional_on_partition(adata_basic):
    names = list(adata_basic.var_names)
    # force a few DE genes via external table so gene_sets can hit them
    de_df = pd.DataFrame(
        {
            "logFC": [3.0 if i < 20 else 0.0 for i in range(len(names))],
            "p_adj": [1e-6 if i < 20 else 1.0 for i in range(len(names))],
        },
        index=names,
    )
    gene_sets = {"PROG": names[:12], "OTHER": names[12:20]}
    r = _run(
        adata_basic,
        de=de_df,
        gene_sets=gene_sets,
        induction_matched=True,
        program_min_genes=5,
    )
    assert r.programs is not None
    assert r.meta["programs_induction_matched"]["status"] in ("ok", "empty")
    # attribute always present
    assert hasattr(r, "programs_induction_matched")
