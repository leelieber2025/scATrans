#!/usr/bin/env python3
"""
Tutorial: the "DE selects, proxy annotates" workflow (scATrans >= this release).

Demonstrates the decoupled workflow and every function it added:

  1. scat.qc.regime_diagnosis(adata)         -- proxy-reliability pre-flight
  2. run_default_pipeline(select_by="de",     -- DE SELECTS the gene list,
                          annotate_mechanism=True)    proxy ANNOTATES mechanism
  3. scat.threshold_sensitivity(...)          -- how list size moves with cutoffs
  4. scat.program_mechanism(..., gene_sets)   -- threshold-free, program-level
                                                 transcription-vs-stabilization call

It runs TWO demos:
  A) a small SYNTHETIC object with known ground truth (a transcription-driven and
     a stabilization-driven program) -> shows the functions return the RIGHT calls;
  B) the real scNT-seq KCl-neuron object if present
     (~/cellranger/scNT_work/scNT_KCl_velocity_plus_labeled.h5ad) -> shows it runs
     on real data AND how the regime pre-flight flags a low-unspliced-capture
     dataset where the per-gene proxy direction should NOT be trusted.

Run:
    python examples/select_annotate_workflow_example.py
"""

from __future__ import annotations

import os

import numpy as np
import anndata as ad
import scipy.sparse as sp

import scatrans as scat

RANDOM = np.random.default_rng(0)
SCNT_PATH = os.path.expanduser("~/cellranger/scNT_work/scNT_KCl_velocity_plus_labeled.h5ad")

# Classic KCl / activity-induced immediate-early genes (mouse); their up-regulation
# is transcription-driven (a burst of new transcription).
IEG_PROGRAM = [
    "Fos", "Fosb", "Junb", "Egr1", "Egr2", "Arc", "Npas4",
    "Bdnf", "Nr4a1", "Dusp1", "Nr4a2", "Nr4a3", "Egr4", "Gadd45b",
]


def make_synthetic(n_cells: int = 800, n_genes: int = 400):
    """Disease vs Control with KNOWN ground truth:
      g0..g19  transcription-driven (up in spliced AND unspliced),
      g20..g29 stabilization-like    (up in spliced only),
      rest     unchanged.
    """
    cond = np.array(["Disease"] * (n_cells // 2) + ["Control"] * (n_cells // 2))
    spl = RANDOM.poisson(3, (n_cells, n_genes)).astype(float)
    uns = RANDOM.poisson(1, (n_cells, n_genes)).astype(float)
    dis = cond == "Disease"
    for j in range(20):
        spl[dis, j] += RANDOM.poisson(8, dis.sum())
        uns[dis, j] += RANDOM.poisson(6, dis.sum())
    for j in range(20, 30):
        spl[dis, j] += RANDOM.poisson(8, dis.sum())
    a = ad.AnnData(
        X=sp.csr_matrix(spl),
        layers={"spliced": sp.csr_matrix(spl), "unspliced": sp.csr_matrix(uns)},
    )
    a.obs["condition"] = cond
    a.var_names = [f"g{i}" for i in range(n_genes)]
    gene_sets = {
        "induced_transcription": [f"g{i}" for i in range(20)],
        "stabilization_set": [f"g{i}" for i in range(20, 30)],
        "random_control": [f"g{i}" for i in range(100, 130)],
    }
    return a, "Disease", "Control", gene_sets


def load_real_scnt(max_cells: int = 3000, min_counts: int = 40):
    """scNT KCl neurons, 0min (reference) vs 15min (target); subset for speed."""
    adata = ad.read_h5ad(SCNT_PATH)
    adata = adata[adata.obs["condition"].isin(["Neu-Kcl-0min", "Neu-Kcl-15min"])].copy()
    tot = np.asarray(adata.layers["spliced"].sum(0) + adata.layers["unspliced"].sum(0)).ravel()
    gkeep = tot >= min_counts
    for g in IEG_PROGRAM:
        if g in adata.var_names:
            gkeep[adata.var_names.get_loc(g)] = True
    adata = adata[:, gkeep].copy()
    if adata.n_obs > max_cells:
        adata = adata[np.sort(RANDOM.choice(adata.n_obs, max_cells, replace=False))].copy()
    gene_sets = {
        "KCl_immediate_early": [g for g in IEG_PROGRAM if g in adata.var_names],
        "random_control": list(RANDOM.choice(adata.var_names, 30, replace=False)),
    }
    return adata, "Neu-Kcl-15min", "Neu-Kcl-0min", gene_sets


def run_demo(adata, target, reference, gene_sets, title, caveat=None):
    print("\n" + "#" * 72)
    print("#", title)
    print("#" * 72)
    print(f"{adata.n_obs} cells x {adata.n_vars} genes; contrast {target} vs {reference}")

    # 1) pre-flight regime / proxy reliability
    regime = scat.qc.regime_diagnosis(adata)
    print(f"\n1) regime_diagnosis: unspliced={regime['unspliced_fraction']:.1%}  "
          f"regime={regime['regime']}  reliability={regime['reliability']:.2f}")
    print(f"   {regime['message']}")

    # 2) DE selects, proxy annotates
    res = scat.run_default_pipeline(
        adata, groupby="condition", target_group=target, reference_group=reference,
        organism="mouse", run_go_enrichment=False,
        select_by="de", annotate_mechanism=True,
    )
    cand = res.candidates
    print(f"\n2) select_by='de': {len(cand)} DE candidates (padj<0.05 & |log2FC|>1, up)")
    if "mechanism_class" in res.all_results.columns and len(cand):
        counts = res.all_results.loc[cand.index, "mechanism_class"].value_counts().to_dict()
        print(f"   mechanism_class among candidates: {counts}")
        cols = [c for c in ("logFC", "p_adj", "transcription_support",
                            "mechanism_class", "mechanism_confidence") if c in cand.columns]
        print(res.all_results.loc[cand.index, cols].head(8).to_string())

    # 3) threshold sensitivity
    print("\n3) threshold_sensitivity:")
    ts = scat.threshold_sensitivity(res.all_results, padj_grid=(0.01, 0.05, 0.1),
                                    logfc_grid=(0.58, 1.0, 1.5))
    print(ts.to_string(index=False))

    # 4) program-level, threshold-free mechanism
    print("\n4) program_mechanism (threshold-free pooling):")
    gene_sets = {k: v for k, v in gene_sets.items() if len(v) >= 5}
    pm = scat.program_mechanism(res.all_results, gene_sets, min_genes=5)
    print("   (no program met min_genes)" if pm.empty else
          pm[["program", "n_genes", "mean_support", "bg_mean_support",
              "p_value", "fdr", "direction"]].to_string(index=False))
    if caveat:
        print(f"\n   CAVEAT: {caveat}")


def main():
    # A) synthetic with known ground truth -> the functions should return the
    #    correct transcription- vs stabilization-driven calls.
    a, tgt, ref, gs = make_synthetic()
    run_demo(a, tgt, ref, gs, "A) SYNTHETIC (known ground truth)")

    # B) real data (optional): shows it runs, and the regime pre-flight catching a
    #    low-capture dataset where the per-gene proxy direction is unreliable.
    if os.path.exists(SCNT_PATH):
        a2, tgt2, ref2, gs2 = load_real_scnt()
        run_demo(
            a2, tgt2, ref2, gs2,
            "B) REAL scNT-seq KCl neurons (15min vs 0min)",
            caveat=(
                "scNT is Drop-seq -> low unspliced capture (regime=low_unspliced). "
                "The proxy is noise-dominated here, so the mechanism DIRECTION for "
                "the IEGs is unreliable (they are truly transcription-driven). This "
                "is exactly why regime_diagnosis down-weights the annotation. A "
                "high-capture protocol (10x with introns / plate-based) is needed to "
                "trust the per-gene call — the synthetic demo (A) shows the correct "
                "behavior when capture is adequate."
            ),
        )
    else:
        print(f"\n(real scNT object not found at {SCNT_PATH}; ran synthetic only)")


if __name__ == "__main__":
    main()
