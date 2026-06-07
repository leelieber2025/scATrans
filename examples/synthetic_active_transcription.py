#!/usr/bin/env python3
"""
Synthetic example demonstrating the core scATrans workflow.

This script generates random velocity-like data, runs active_score,
and produces publication-style plots (including use of ax= for multi-panel
figures). It uses only synthetic data and requires no external h5ad files.

Run:
    python examples/synthetic_active_transcription.py
"""

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

import scatrans as scat

# Use the improved professional style (inspired by OmicVerse etc.)
scat.pl.set_style()

print("Generating synthetic AnnData with spliced/unspliced layers...")

np.random.seed(42)
n_cells = 120
n_genes = 300

# Simulate count data
X = np.random.negative_binomial(4, 0.45, size=(n_cells, n_genes)).astype(float)

adata = sc.AnnData(X)
adata.obs["condition"] = ["Disease"] * 60 + ["Control"] * 60

# Add velocity layers
adata.layers["spliced"] = X.copy()
adata.layers["unspliced"] = X * 0.55

# Add plausible gene features (so bias correction has something to do)
adata.var["gene_length"] = np.random.randint(700, 5000, n_genes)
adata.var["intron_number"] = np.random.randint(0, 12, n_genes)

print("Running active_score (heuristic + small permutation for demo)...")

adata_res, sig, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",
    target_group="Disease",
    reference_group="Control",
    mode="heuristic",
    use_permutation=True,
    n_perm=20,           # small for speed in example
    show_plot=False,     # we will create custom figure below
    min_total_counts=30,
)

print(f"Found {len(sig)} significant genes (demo data).")
print("Top genes by active score:")
print(sig.head(5))

# ------------------------------------------------------------------
# Custom multi-panel figure demonstrating ax= support
# ------------------------------------------------------------------
print("\nCreating custom multi-panel figure using ax= parameter...")

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

# Left: comet plot embedded in subplot
scat.pl.comet_plot(
    all_results,
    top_n=8,
    title="Comet Plot (ax= demo)",
    ax=axes[0],
)

# Right: simple volcano (or bias diagnostic)
scat.pl.volcano_plot(
    all_results,
    top_n=6,
    title="Volcano Plot (ax= demo)",
    ax=axes[1],
)

plt.suptitle("scATrans synthetic example — custom layout with ax=", fontsize=14, y=1.02)
plt.tight_layout()

# Save the combined figure
out_path = "examples/synthetic_multi_panel.pdf"
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved multi-panel figure to {out_path}")

# Also demonstrate bias diagnostic (standalone)
print("\nGenerating bias diagnostic plot...")
scat.pl.bias_diagnostic_plot(
    all_results,
    title="Bias Correction Diagnostic (synthetic)",
    save_path="examples/synthetic_bias_diagnostic.pdf",
)

print("\nExample complete. Key outputs:")
print("  - examples/synthetic_multi_panel.pdf")
print("  - examples/synthetic_bias_diagnostic.pdf")
print("\nYou can open the PDFs to see the style (vector-friendly, clean, journal-ready).")

if __name__ == "__main__":
    pass  # script is executable directly
