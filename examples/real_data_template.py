#!/usr/bin/env python3
"""
Real-data workflow template for scATrans (recommended practice as of 2026-06).

This is NOT a runnable script out of the box — it shows the full recommended
workflow with heavy comments for users who have their own spliced/unspliced data
(from velocyto, kb_python --lamanno/velocity mode, etc.).

Typical real-data steps covered:
  1. Load AnnData with layers
  2. (Optional but recommended) Standalone QC
  3. Attach gene features for bias correction
  4. Run active_score with diagnostics-friendly settings (permutation on, advanced when appropriate)
  5. Inspect the rich new diagnostics in adata.uns["scatrans"]
  6. Generate publication figures (comet, bias diagnostic) using ax= where helpful
  7. Run enrichment on significant genes

Copy this file, adapt the paths / column names, and run cell-by-cell in a notebook.

See also:
  - README.md "Choosing mode" section and "Quick data quality check"
  - synthetic_active_transcription.py for a fully runnable (but synthetic) demo
"""

import matplotlib.pyplot as plt
import scanpy as sc

import scatrans as scat

# ------------------------------------------------------------------
# 1. Load your data (example: kb_python or velocyto output)
# ------------------------------------------------------------------
# adata = sc.read_h5ad("path/to/your_velocity_object.h5ad")
# Common patterns:
#   - kb_python lamanno/velocity → layers 'mature' + 'nascent' (auto-remapped)
#   - velocyto / scVelo standard → 'spliced' + 'unspliced'
#
# Make sure you have a categorical obs column for the biological contrast, e.g.
# adata.obs["condition"] = ...   # "Disease", "Control", etc.

# For illustration we assume the object is already loaded and has the required layers + obs.
# adata = ...

# ------------------------------------------------------------------
# 2. Standalone pre-flight QC (HIGHLY RECOMMENDED)
# ------------------------------------------------------------------
print("Running global unspliced fraction check...")
ufrac = scat.qc.unspliced_global(adata, warn_threshold=0.5)
print(f"Global unspliced fraction: {ufrac:.2%}")
# If this is > ~0.5 you should investigate library prep / alignment before trusting active transcription signals.

# ------------------------------------------------------------------
# 3. Attach gene features (length + intron count) for bias correction
# ------------------------------------------------------------------
# The package ships mouse tables. For human or custom annotations use:
#   scat.add_gene_features(adata, gene_features_path="my_features.parquet")
# or the CLI: generate-gene-features --gtf genes.gtf --output features.parquet
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# ------------------------------------------------------------------
# 4. Run the analysis — recommended settings for real data
# ------------------------------------------------------------------
print("\nRunning active_score with diagnostics enabled (permutation + heuristic by default)...")

# Heuristic is robust for most real case/control or stimulus experiments.
# Switch to mode="advanced" only when you have sufficient cells and want moments smoothing.
# See the "Choosing mode" section in the README for the full decision guide.

adata_res, significant, all_results = scat.active_score(
    adata_input=adata,
    groupby="condition",          # <-- change to your column
    target_group="Disease",       # <-- your group of interest
    reference_group="Control",    # <-- your reference
    mode="heuristic",             # or "advanced" when appropriate
    use_permutation=True,
    n_perm=200,                   # 100-500 typical; auto-reduced for tiny pseudobulk designs
    show_plot=False,              # we will make nicer figures below with ax=
    min_total_counts=50,
    # advanced_... parameters only matter if you choose mode="advanced"
)

print(f"\nSignificant active genes: {len(significant)}")
print("Top 8:")
print(significant.head(8))

# ------------------------------------------------------------------
# 5. Inspect the rich diagnostics (new in this improvement cycle)
# ------------------------------------------------------------------
meta = adata_res.uns["scatrans"]
print("\n=== Key diagnostics (always inspect these!) ===")
print("Unspiced global fraction :", meta.get("unspliced_global_fraction"))
print("Bias correction info     :", meta["diagnostics"]["bias_correction"])
print("Permutation approximation:", meta.get("permutation_approximation_note"))
print("Full diagnostics dict is at adata_res.uns['scatrans']['diagnostics']")
print("Note: by default effective_gamma and delta_variance are not added to .var (opt-in via flags).")

# Effective gamma per gene (transparency) — only present if show_effective_gamma=True was used
if "effective_gamma" in adata_res.var.columns:
    print("effective_gamma (first 5):", adata_res.var["effective_gamma"].head().tolist())
else:
    print("effective_gamma not exposed (default). Pass show_effective_gamma=True to active_score to include it.")

# ------------------------------------------------------------------
# 6. Publication-quality figures (ax= support for multi-panel)
# ------------------------------------------------------------------
scat.pl.set_style()   # call once for clean vector output (Type 42 fonts etc.)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=150)

# Main result view
scat.pl.comet_plot(
    all_results,
    top_n=10,
    title="Active Transcription Drivers (real data example)",
    ax=axes[0],
)

# Bias correction diagnostic (before/after) — very important to show reviewers
scat.pl.bias_diagnostic_plot(
    all_results,
    title="Bias correction (length + intron number)",
    axes=(axes[1], None),   # or pass two axes
    show_regression=True,
)

plt.tight_layout()
fig.savefig("real_data_comet_and_bias.pdf", dpi=300, bbox_inches="tight")
print("\nSaved real_data_comet_and_bias.pdf")

# Optional: rank plot or volcano
# scat.pl.active_score_rankplot(all_results, top_n=15, save_path="ranks.pdf")
# scat.pl.volcano_plot(all_results, top_n=8, label_genes=some_gene_list, save_path="volcano.pdf")  # ggVolcano-like control

# ------------------------------------------------------------------
# 7. Functional enrichment on the significant genes
# ------------------------------------------------------------------
if len(significant) > 0:
    enrich = scat.run_enrichment(
        gene_list=significant.index.tolist(),
        gene_sets="GO_Biological_Process_2023",  # defaults to bundled scATrans version
        organism="mouse",   # or "human"
        universe=adata.var_names.tolist(),  # or background= (compat); default = conservative intersect like clusterProfiler
        pval_cutoff=0.05,
    )
    # To use a specific Enrichr historical version, just write the full name:
    # enrich = scat.run_enrichment(..., gene_sets="GO_Biological_Process_2021")
    # For KEGG: scat.run_kegg(..., kegg_library="KEGG_2021")
    print("\nTop enrichment terms:")
    print(enrich.head(6))
    print("universe_info:", enrich.attrs.get("universe_info"))
    # scat.pl.enrich_dotplot(enrich, show_terms=10, save_path="enrich.pdf")
    # or scat.pl.enrich_dotplot(enrich, show_terms=["specific term desc", "another GO term"])

print("\n=== Recommended next steps ===")
print("- Look at adata_res.var for 'active_score', 'velocity_residual', 'effective_gamma', logFC, p_adj, ...")
print("- Check the full diagnostics dict for any red flags (small bias fit n, very high unspliced frac, etc.).")
print("- For top genes, you may want to manually inspect their phase portraits (U vs S colored by group).")
print("- Consider running the same contrast with mode='advanced' (if cell number permits) and comparing the resulting significant lists.")
print("- Always report the exact mode, n_perm, prior_weight, and whether permutation was used in your methods.")

if __name__ == "__main__":
    print("\nThis template is meant to be read and adapted, not executed directly.")