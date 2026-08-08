"""Tests for the gene-level UpSet feature.

Covers the three-piece workflow added to ``scatrans.pl``:

    build_gene_membership  ->  gene_upsetplot  ->  common_genes

Deterministic gene sets with known overlaps make the intersection assertions
exact. All figures use the Agg backend.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba

import scatrans as scat

GENES = [f"G{i}" for i in range(10)]


def _make_de(up: list[str], down: list[str]) -> pd.DataFrame:
    """A minimal DE table: listed genes get a strong signed logFC, the rest 0.

    All p_adj are tiny so the direction/logFC cutoff alone decides membership.
    """
    logfc = dict.fromkeys(GENES, 0.0)
    for g in up:
        logfc[g] = 1.0
    for g in down:
        logfc[g] = -1.0
    return pd.DataFrame(
        {"logFC": [logfc[g] for g in GENES], "p_adj": [1e-3] * len(GENES)},
        index=GENES,
    )


@pytest.fixture
def three_de() -> dict[str, pd.DataFrame]:
    # G0,G1 up in every model; G8,G9 down in every model; plus model-specific ones.
    return {
        "A": _make_de(up=["G0", "G1", "G2"], down=["G7", "G8", "G9"]),
        "B": _make_de(up=["G0", "G1", "G3"], down=["G6", "G8", "G9"]),
        "C": _make_de(up=["G0", "G1", "G4"], down=["G5", "G8", "G9"]),
    }


# --------------------------------------------------------------------------- #
# build_gene_membership
# --------------------------------------------------------------------------- #
def test_membership_separate_splits_up_down(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    assert list(mem.columns) == [
        "A::up",
        "A::down",
        "B::up",
        "B::down",
        "C::up",
        "C::down",
    ]
    # values are strictly 0/1
    assert set(np.unique(mem.values)) <= {0, 1}
    # A::up must contain exactly G0,G1,G2
    assert set(mem.index[mem["A::up"] == 1]) == {"G0", "G1", "G2"}
    assert set(mem.index[mem["A::down"] == 1]) == {"G7", "G8", "G9"}
    # per-set gene lists are stashed in attrs
    assert set(mem.attrs["gene_sets"]["A::up"]) == {"G0", "G1", "G2"}


def test_membership_single_direction_one_set_per_model(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="up", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    assert list(mem.columns) == ["A", "B", "C"]
    assert set(mem.index[mem["A"] == 1]) == {"G0", "G1", "G2"}


def test_membership_accepts_plain_gene_lists():
    mem = scat.pl.build_gene_membership({"X": ["G1", "G2", "G3"], "Y": ["G2", "G3", "G4"]})
    assert list(mem.columns) == ["X", "Y"]
    assert set(mem.index) == {"G1", "G2", "G3", "G4"}
    assert set(mem.index[mem["X"] == 1]) == {"G1", "G2", "G3"}


def test_membership_bad_input_raises():
    with pytest.raises(ValueError):
        scat.pl.build_gene_membership(["not", "a", "mapping"])  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        scat.pl.build_gene_membership({"A": _make_de(["G0"], [])}, direction="sideways")


# --------------------------------------------------------------------------- #
# common_genes
# --------------------------------------------------------------------------- #
def test_common_genes_up_and_down(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    assert scat.pl.common_genes(mem, direction="up") == ["G0", "G1"]
    assert scat.pl.common_genes(mem, direction="down") == ["G8", "G9"]


def test_common_genes_min_sets_relaxes_intersection(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    # G2 is up only in A, G3 only in B, G4 only in C -> none shared by >=2 beyond G0,G1
    relaxed = scat.pl.common_genes(mem, direction="up", min_sets=2)
    assert set(relaxed) == {"G0", "G1"}
    # min_sets=1 = union of all up genes
    union_up = scat.pl.common_genes(mem, direction="up", min_sets=1)
    assert set(union_up) == {"G0", "G1", "G2", "G3", "G4"}


def test_common_genes_explicit_sets(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    # genes up in A and down in ... none; intersect A::up & B::up -> G0,G1
    assert scat.pl.common_genes(mem, sets=["A::up", "B::up"]) == ["G0", "G1"]


# --------------------------------------------------------------------------- #
# gene_upsetplot
# --------------------------------------------------------------------------- #
def test_gene_upsetplot_from_de_results(three_de):
    fig, ax = scat.pl.gene_upsetplot(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5, show=False
    )
    assert fig is not None and ax is not None
    plt.close(fig)


def test_gene_upsetplot_from_prebuilt_membership(three_de):
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    fig, ax = scat.pl.gene_upsetplot(membership=mem, show=False)
    plt.close(fig)


def test_gene_upsetplot_requires_some_input():
    with pytest.raises(ValueError):
        scat.pl.gene_upsetplot(show=False)


def test_gene_upsetplot_empty_membership_is_graceful():
    empty = pd.DataFrame()
    fig, ax = scat.pl.gene_upsetplot(membership=empty, show=False)
    assert fig is not None
    plt.close(fig)


def test_gene_upsetplot_custom_colors_applied(three_de):
    """set_color hits the set-size bars; a per-column intersection_color list
    colors individual intersection bars (exercises the highlight path)."""
    mem = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    inter_colors = ["#EE7733", "#009988"] + ["#BBBBBB"] * 30
    fig, _ = scat.pl.gene_upsetplot(
        membership=mem, set_color="#8888CC", intersection_color=inter_colors, show=False
    )
    # axes created in order: set-size bars, intersection bars, matrix
    ax_set, ax_inter = fig.axes[0], fig.axes[1]
    assert to_rgba(ax_set.patches[0].get_facecolor()) == to_rgba("#8888CC")
    assert to_rgba(ax_inter.patches[0].get_facecolor()) == to_rgba("#EE7733")
    assert to_rgba(ax_inter.patches[1].get_facecolor()) == to_rgba("#009988")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# normalize_external_de / external tables
# --------------------------------------------------------------------------- #
def test_normalize_external_de_deseq2_style():
    raw = pd.DataFrame(
        {
            "gene": ["A", "B", "C"],
            "log2FoldChange": [1.2, -0.8, 0.1],
            "padj": [1e-4, 1e-3, 0.2],
            "pvalue": [1e-5, 1e-4, 0.1],
        }
    )
    de = scat.pl.normalize_external_de(raw)
    assert list(de.index) == ["A", "B", "C"]
    assert "logFC" in de.columns and "p_adj" in de.columns
    assert de.loc["A", "logFC"] == pytest.approx(1.2)
    assert de.loc["B", "p_adj"] == pytest.approx(1e-3)


def test_normalize_external_de_seurat_style_index_gene():
    raw = pd.DataFrame(
        {"avg_log2FC": [1.0, -1.0], "p_val_adj": [0.01, 0.02]},
        index=["X", "Y"],
    )
    de = scat.pl.normalize_external_de(raw)
    assert de.loc["X", "logFC"] == pytest.approx(1.0)
    assert de.loc["Y", "p_adj"] == pytest.approx(0.02)


def test_build_membership_from_external_deseq2():
    a = pd.DataFrame(
        {
            "Gene": ["G0", "G1", "G2", "G8"],
            "log2FoldChange": [1.0, 1.0, 1.0, -1.0],
            "padj": [1e-3, 1e-3, 1e-3, 1e-3],
        }
    )
    b = pd.DataFrame(
        {
            "gene": ["G0", "G1", "G3", "G9"],
            "log2FoldChange": [1.0, 1.0, 1.0, -1.0],
            "p.adjust": [1e-3, 1e-3, 1e-3, 1e-3],
        }
    )
    mem = scat.pl.build_gene_membership(
        {"fracture": a, "osteomyelitis": b},
        direction="separate",
        pval_cutoff=0.05,
        logfc_cutoff=0.5,
    )
    assert set(mem.index[mem["fracture::up"] == 1]) == {"G0", "G1", "G2"}
    assert set(mem.index[mem["osteomyelitis::up"] == 1]) == {"G0", "G1", "G3"}
    assert scat.pl.common_genes(mem, direction="up") == ["G0", "G1"]
    assert scat.pl.common_genes(mem, direction="down") == []


def test_already_filtered_skips_pvalue_gate():
    # only logFC, no padj — pre-filtered DEG list
    a = pd.DataFrame({"gene": ["U1", "D1"], "log2FoldChange": [2.0, -2.0]})
    b = pd.DataFrame({"gene": ["U1", "D1"], "log2FoldChange": [1.5, -1.5]})
    mem = scat.pl.build_gene_membership(
        {"A": a, "B": b},
        direction="separate",
        logfc_cutoff=0.5,
        already_filtered=True,
    )
    assert set(mem.index[mem["A::up"] == 1]) == {"U1"}
    assert set(mem.index[mem["A::down"] == 1]) == {"D1"}
    assert scat.pl.common_genes(mem, direction="up") == ["U1"]
    assert scat.pl.common_genes(mem, direction="down") == ["D1"]


def test_p_type_nominal_vs_adj():
    """p_type chooses p_val vs p_adj when both columns exist."""
    df = pd.DataFrame(
        {
            "logFC": [1.0, 1.0, -1.0],
            "p_val": [1e-4, 0.04, 1e-4],
            "p_adj": [0.2, 0.04, 0.2],  # only middle gene is adj-significant
        },
        index=["U_nom", "U_both", "D_nom"],
    )
    mem_adj = scat.pl.build_gene_membership(
        {"A": df},
        direction="separate",
        pval_cutoff=0.05,
        logfc_cutoff=0.5,
        p_type="adj",
    )
    assert set(mem_adj.index[mem_adj["A::up"] == 1]) == {"U_both"}
    assert set(mem_adj.index[mem_adj["A::down"] == 1]) == set()

    mem_nom = scat.pl.build_gene_membership(
        {"A": df},
        direction="separate",
        pval_cutoff=0.05,
        logfc_cutoff=0.5,
        p_type="nominal",
    )
    assert set(mem_nom.index[mem_nom["A::up"] == 1]) == {"U_nom", "U_both"}
    assert set(mem_nom.index[mem_nom["A::down"] == 1]) == {"D_nom"}


def test_collect_gene_sets_two_and_three(three_de):
    mem2 = scat.pl.build_gene_membership(
        {"fracture": three_de["A"], "osteomyelitis": three_de["B"]},
        direction="separate",
        pval_cutoff=0.05,
        logfc_cutoff=0.5,
    )
    gs2 = scat.pl.collect_gene_sets(mem2)
    assert "fracture_up" in gs2 and "osteomyelitis_down" in gs2
    assert "common_up" in gs2 and "common_down" in gs2
    assert set(gs2["common_up"]) == {"G0", "G1"}
    assert set(gs2["common_down"]) == {"G8", "G9"}

    mem3 = scat.pl.build_gene_membership(
        three_de, direction="separate", pval_cutoff=0.05, logfc_cutoff=0.5
    )
    gs3 = scat.pl.collect_gene_sets(mem3)
    assert set(gs3["common_up"]) == {"G0", "G1"}
    assert set(gs3["common_down"]) == {"G8", "G9"}
    # three conditions → 3 up + 3 down + 2 common = 8 keys
    assert len(gs3) == 8
