#!/usr/bin/env python
"""Example: the DE→mechanism PRIMARY workflow (partition_de_by_mechanism).

scATrans does NOT compete with differential expression for gene discovery. The
recommended entry runs a standard DE to SELECT the changed genes, then partitions
those genes into transcription-driven vs stabilization-driven by the nascent
unspliced-excess signal — decisively at the program level.

Two demos:
  (A) synthetic known-truth  — induced-by-transcription vs induced-by-stabilization
  (B) pluggable DE front-end — same call, different DE sources (builtin / method
      name / precomputed DataFrame / callable)

Run:  PYTHONPATH=src python examples/partition_de_by_mechanism_example.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc

import scatrans as scat


def _synthetic(seed: int = 0):
    """Two groups; some genes induced via transcription (high unspliced), some
    via stabilization (high spliced, low unspliced), some unchanged."""
    rng = np.random.default_rng(seed)
    n = 400
    n_txn, n_stab, n_bg = 30, 30, 190
    genes = (
        [f"TXN{i}" for i in range(n_txn)]
        + [f"STAB{i}" for i in range(n_stab)]
        + [f"BG{i}" for i in range(n_bg)]
    )
    g = len(genes)

    def block(mu_s_ctrl, mu_s_dis, mu_u_ctrl, mu_u_dis, k):
        s_c = rng.poisson(mu_s_ctrl, (n, k)); s_d = rng.poisson(mu_s_dis, (n, k))
        u_c = rng.poisson(mu_u_ctrl, (n, k)); u_d = rng.poisson(mu_u_dis, (n, k))
        return s_c, s_d, u_c, u_d

    # transcription-driven: both spliced AND unspliced rise
    sTc, sTd, uTc, uTd = block(4, 16, 3, 14, n_txn)
    # stabilization-driven: spliced rises, unspliced barely moves
    sSc, sSd, uSc, uSd = block(4, 16, 3, 4, n_stab)
    # background: unchanged
    sBc, sBd, uBc, uBd = block(5, 5, 4, 4, n_bg)

    S_ctrl = np.hstack([sTc, sSc, sBc]); S_dis = np.hstack([sTd, sSd, sBd])
    U_ctrl = np.hstack([uTc, uSc, uBc]); U_dis = np.hstack([uTd, uSd, uBd])
    S = np.vstack([S_dis, S_ctrl]).astype(float)
    U = np.vstack([U_dis, U_ctrl]).astype(float)

    ad = sc.AnnData(S)
    ad.var_names = genes
    ad.obs["condition"] = ["Disease"] * n + ["Control"] * n
    ad.layers["spliced"] = S
    ad.layers["unspliced"] = U
    ad.var["gene_length"] = rng.integers(700, 4500, g)
    ad.var["intron_number"] = rng.integers(1, 12, g)
    return ad


def main():
    ad = _synthetic()
    txn = [g for g in ad.var_names if g.startswith("TXN")]
    stab = [g for g in ad.var_names if g.startswith("STAB")]

    print("=" * 70)
    print("(A) partition_de_by_mechanism — builtin DE + program-level call")
    print("=" * 70)
    r = scat.partition_de_by_mechanism(
        ad, groupby="condition", target_group="Disease", reference_group="Control",
        organism="human", gene_sets={"transcription_program": txn,
                                     "stabilization_program": stab},
        program_min_genes=5,
        # toy caveat: here the two programs ARE the whole selected list, so pool
        # the program test over all tested genes (else the background is empty).
        program_restrict_to_selected=False,
    )
    print("regime:", r.regime["regime"], "| reliability:", r.regime["reliability"])
    print("selected (DE):", len(r.selected), "genes")
    print("\nper-gene mechanism_class counts among selected:")
    print(r.selected["mechanism_class"].value_counts().to_string())
    print("\nPROGRAM-LEVEL call (the decisive layer):")
    cols = ["program", "n_genes", "mean_support", "direction", "p_value", "fdr"]
    if r.programs is not None and len(r.programs):
        print(r.programs[[c for c in cols if c in r.programs.columns]].to_string(index=False))
    else:
        print("(no programs met min_genes / background — see program_restrict_to_selected)")
    print(
        "\nexpected: transcription_program -> transcription-driven; "
        "stabilization_program -> stabilization-driven"
    )

    print("\n" + "=" * 70)
    print("(B) pluggable DE front-end — same call, different DE sources")
    print("=" * 70)
    # builtin
    n_builtin = len(scat.partition_de_by_mechanism(
        ad, target_group="Disease", reference_group="Control", organism="human").selected)
    # a scanpy DE method routed through differential_expression
    n_ttest = len(scat.partition_de_by_mechanism(
        ad, target_group="Disease", reference_group="Control", organism="human",
        de="t-test").selected)
    # a precomputed DataFrame (pretend it came from edgeR/DESeq2)
    de_df = pd.DataFrame(
        {"log2FC": [3.0 if (g.startswith("TXN") or g.startswith("STAB")) else 0.0
                    for g in ad.var_names],
         "FDR": [1e-5 if (g.startswith("TXN") or g.startswith("STAB")) else 1.0
                 for g in ad.var_names]},
        index=list(ad.var_names),
    )
    n_ext = len(scat.partition_de_by_mechanism(
        ad, target_group="Disease", reference_group="Control", organism="human",
        de=de_df, de_logfc_col="log2FC", de_padj_col="FDR").selected)
    print(f"selected genes — builtin: {n_builtin} | t-test: {n_ttest} | precomputed DF: {n_ext}")
    print("(same partition/annotation applied regardless of which DE selected the list)")


if __name__ == "__main__":
    main()
