#!/usr/bin/env python
"""
Example demonstrating Memento integration in scATrans.

Primary focus (main line):
- Standard active transcription analysis with active_score on data that HAS spliced/unspliced layers.
- How to choose between traditional DE (t-test / wilcoxon via de_method) and Memento inside the velocity workflow.

Additional capability (shown at the end):
- Using the package for pure differential expression when you have NO velocity data at all
  (via the new differential_expression() entry point, which also supports Memento).

This keeps the original scATrans active transcription story as the main narrative,
while documenting the new Memento / general DE options as add-ons.

Requires (for full demo):
  pip install "scatrans[memento]"   # pulls memento-de
  (the rest are core dependencies)
"""

import numpy as np
import pandas as pd
import anndata as ad
import scatrans as scat

print("=== scATrans + Memento Example ===")
print(f"scATrans version: {scat.__version__}")

# =============================================================================
# PART 1: Main line — data WITH spliced/unspliced layers (recommended original workflow)
# =============================================================================

print("\n" + "="*70)
print("PART 1: Original main workflow (data has velocity layers)")
print("Focus: active_score + choosing DE backend (t-test vs Memento)")
print("="*70)

np.random.seed(42)
n_cells = 180
n_genes = 60

obs = pd.DataFrame({
    "condition": ["Ctrl"] * 90 + ["Disease"] * 90,
    "batch": np.repeat(["b1", "b2"], 90),
})
obs.index = [f"cell_{i}" for i in range(n_cells)]

# Synthetic counts with some active transcription signal in unspliced
base = np.random.poisson(1.5, size=(n_cells, n_genes)).astype(float)
effect_genes = [0, 1, 2]
for g in effect_genes:
    base[90:, g] += np.random.poisson(3.5, size=90)

spliced = np.clip((base * 0.75).astype(int), 0, None)
unspliced = np.clip((base * 0.25 + np.random.poisson(0.8, size=base.shape)).astype(int), 0, None)

adata = ad.AnnData(
    X=spliced.copy().astype(float),
    obs=obs,
    var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)]),
)
adata.layers["spliced"] = spliced
adata.layers["unspliced"] = unspliced
adata.layers["counts"] = (spliced + unspliced).astype(float)

print(f"Created synthetic velocity AnnData: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"Groups: {adata.obs['condition'].value_counts().to_dict()}")

# Store raw counts + original spliced/unspliced layers early
# (before any HVG/normalize/log). This is the recommended early step for scATrans.
# The sidecar snapshot preserves the full-gene counts AND velocity layers so they
# survive later HVG/cell subsetting; recover them with
# scat.restore_raw_counts(adata, full_genes=True).
scat.store_raw_counts(adata, layer="counts", save_raw=False)
print("Called store_raw_counts early (full-gene counts + velocity preserved in snapshot)")

# --- Standard run with traditional DE (t-test) ---
print("\n--- Standard active_score with default / t-test DE backend ---")
adata_t, sig_t, all_t = scat.active_score(
    adata_input=adata.copy(),
    groupby="condition",
    target_group="Disease",
    reference_group="Ctrl",
    de_method="t-test_overestim_var",   # explicit traditional choice
    de_preprocess="auto",
    use_permutation=False,
    show_plot=False,
    min_total_counts=5,
)
print("Traditional DE columns sample:", ["logFC", "p_adj"])
print(all_t.head(3)[["active_score", "logFC", "p_adj", "velocity_residual"]])

# --- Run with Memento as the DE backend inside active_score ---
print("\n--- active_score + Memento as DE backend (use_memento_de=True) ---")
adata_m, sig_m, all_m = scat.active_score(
    adata_input=adata.copy(),
    groupby="condition",
    target_group="Disease",
    reference_group="Ctrl",
    use_memento_de=True,
    memento_capture_rate=0.10,
    memento_num_boot=300,
    de_preprocess="none",   # Memento prefers count scale
    use_permutation=False,
    show_plot=False,
    min_total_counts=5,
)
print("Memento-specific columns present:", [c for c in all_m.columns if c.startswith("memento_")])
print(all_m.head(3)[["active_score", "logFC", "p_adj", "velocity_residual"]])

meta = adata_m.uns.get("scatrans", {})
print("Memento recorded in metadata:", meta.get("use_memento_de"), meta.get("memento_capture_rate"))

# =============================================================================
# PART 2: Additional capability — pure DE when you have NO velocity data
# =============================================================================

print("\n" + "="*70)
print("PART 2: Additional / non-primary use case")
print("Data has NO spliced/unspliced layers at all.")
print("Use differential_expression() + downstream enrichment/plotting.")
print("This is supported but is not the original intended main workflow.")
print("="*70)

# Plain count AnnData (no velocity layers whatsoever)
np.random.seed(123)
n2, g2 = 120, 80
obs2 = pd.DataFrame({"condition": ["Control"] * 60 + ["Treatment"] * 60})
X2 = np.random.negative_binomial(5, 0.3, size=(n2, g2)).astype(float)
X2[60:, 5:10] += np.random.negative_binomial(20, 0.2, size=(60, 5))  # DE genes

plain = ad.AnnData(X=X2, obs=obs2, var=pd.DataFrame(index=[f"GENE_{i}" for i in range(g2)]))
plain.layers["counts"] = X2.copy()

print(f"Plain count-only AnnData (no velocity): {plain.shape}")

# Store raw early (before any HVG/log) -- required for reliable Memento
scat.store_raw_counts(plain, layer="counts", save_raw=False)

# Choose traditional t-test
plain_t, de_t = scat.differential_expression(
    plain.copy(),
    groupby="condition",
    target_group="Treatment",
    reference_group="Control",
    de_method="t-test_overestim_var",
    de_preprocess="none",
    min_total_counts=10,
)
print("\nPure DE with t-test backend. Top by p_adj:")
print(de_t.head(3)[["logFC", "p_adj"]])

# Choose Memento
plain_m, de_m = scat.differential_expression(
    plain.copy(),
    groupby="condition",
    target_group="Treatment",
    reference_group="Control",
    use_memento_de=True,
    memento_capture_rate=0.08,
    memento_num_boot=200,
    de_preprocess="none",
)
print("\nPure DE with Memento backend. Memento columns:", [c for c in de_m.columns if c.startswith("memento_")])
print(de_m.head(3)[["logFC", "p_adj", "memento_dv_coef"]])

# Downstream tools work the same
cands = scat.filter_active_genes(de_m, pval_cutoff=0.05, logfc_cutoff=0.5)
print(f"\nFiltered candidates from Memento DE: {len(cands)}")

if len(cands) > 0:
    try:
        enrich = scat.run_enrichment(
            gene_list=cands.index.tolist(),
            gene_sets="GO_Biological_Process",  # base name — auto uses the correct
            # bundled organism-specific 2026 set (Mm/Hs). Legacy _scATrans names are
            # also auto-mapped for convenience. Pass adata= so universe is auto from store.
            organism="mouse",  # or "human"; only needed for KEGG or to disambiguate
            adata=plain,
            pval_cutoff=0.2,
        )
        print(f"Enrichment ran, got {len(enrich)} terms.")
    except Exception as e:
        print("Enrichment note:", e)

    try:
        scat.pl.volcano_plot(de_m, top_n=6)
        print("volcano_plot() call succeeded on pure DE results.")
    except Exception as e:
        print("Plot call note:", e)

print("\n=== Example completed successfully ===")
print("Remember: The primary, recommended usage remains active_score on velocity data.")
print("differential_expression + Memento is a powerful additional capability for general DE work.")
