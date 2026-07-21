#!/usr/bin/env python3
"""
Real-data template for the lower-level scATrans scoring path (active_score).

NOTE: the recommended primary workflow is `scat.partition_de_by_mechanism(...)`
(DE selects the changed genes → scATrans partitions them into transcription- vs
stabilization-driven; see examples/partition_de_by_mechanism_example.py and the
Quickstart). This template instead demonstrates the **lower-level** `active_score`
scorer and its diagnostics in detail, for users who want the residual / gamma /
bias internals. It is NOT runnable out of the box — adapt paths / column names.

Steps covered:
  1. Load AnnData with spliced/unspliced layers (velocyto, kb_python velocity, …)
  2. (Optional but recommended) Standalone QC
  3. Attach gene features for bias correction
  4. Run active_score with diagnostics-friendly settings
  5. Inspect the diagnostics in adata.uns["scatrans"]
  6. Generate figures (comet, bias diagnostic) using ax= where helpful
  7. Filter (prefer filter_active_genes(select_by="de")) → enrichment

Copy this file, adapt the paths / column names, and run cell-by-cell in a notebook.

See also:
  - examples/partition_de_by_mechanism_example.py — the recommended workflow
  - the Quickstart and User Guide for the DE→mechanism path and DE backends
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
# 2.5 Preserve raw counts + original spliced/unspliced (CRITICAL for Memento DE and correct enrichment)
# ------------------------------------------------------------------
# Call this early, before any HVG, normalize, or log1p. It saves the full post-QC gene set
# raw counts (and raw velocity layers) so that later Memento, PyDESeq2, and enrichment
# automatically use the true measured background instead of whatever is left in .X after HVG.
scat.store_raw_counts(adata, layer="counts", save_raw=False)
# After this, you can safely do standard Scanpy preprocessing for visualization.
# For DE/enrichment on as many genes as possible, use the adata (or a non-HVG-subset copy) when calling those functions.

# ------------------------------------------------------------------
# 3. Attach gene features (length + intron count) for bias correction
# ------------------------------------------------------------------
# The package ships mouse tables by default.
# For human or custom annotations:
#   1. generate-gene-features --gtf genes.gtf --output my_features.parquet --organism human
#   2. adata = scat.add_gene_features(adata, gene_features_path="my_features.parquet")
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# ------------------------------------------------------------------
# 4. Run the analysis — recommended settings for real data
# ------------------------------------------------------------------
print("\nRunning active_score with diagnostics enabled (permutation + heuristic by default)...")

# Heuristic is robust for most real case/control or stimulus experiments.
# Switch to mode="advanced" only when you have sufficient cells and want moments smoothing.
# See the User Guide (advanced options) and recommend_workflow() for a decision guide.

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

fig = plt.figure(figsize=(16, 5), dpi=150)
grid = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.35)
ax_comet = fig.add_subplot(grid[0, 0])
ax_bias_before = fig.add_subplot(grid[0, 1])
ax_bias_after = fig.add_subplot(grid[0, 2])

# Main result view
scat.pl.comet_plot(
    all_results,
    top_n=10,
    title="Active Transcription Drivers (real data example)",
    ax=ax_comet,
)

# Bias correction diagnostic (before/after) — very important to show reviewers
scat.pl.bias_diagnostic_plot(
    all_results,
    title="Bias correction (length + intron number)",
    axes=(ax_bias_before, ax_bias_after),
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
        gene_sets="GO_Biological_Process",  # or "GO_BP" — auto-resolves to the correct
        # organism-specific bundled built-in (Mm/Hs_GO_Biological_Process_2026.txt).
        # No need to specify year or _scATrans suffix.
        organism="mouse",   # or "human"
        # CRITICAL for correct background: pass adata= (the one on which you called
        # store_raw_counts early) so it auto-uses the preserved full measured gene list
        # instead of whatever is left after HVG. This is the most convenient & correct default.
        adata=adata,
        padj_cutoff=0.05,
    )
    # For KEGG the simplest is:
    # kegg = scat.run_kegg(significant.index.tolist(), organism="mouse", adata=adata)
    # To force a historical Enrichr version instead of the built-in: gene_sets="GO_Biological_Process_2021" or kegg_library="KEGG_2021"
    print("\nTop enrichment terms:")
    print(enrich.head(6))
    print("universe_info:", enrich.attrs.get("universe_info"))
    # scat.pl.enrich_dotplot(enrich, show_terms=10, save_path="enrich.pdf")
    # or scat.pl.enrich_dotplot(enrich, show_terms=["specific term desc", "another GO term"])

print("\n=== Recommended next steps ===")
print("- Look at adata_res.var for 'active_score', 'unspliced_excess_residual', 'effective_gamma', logFC, p_adj, ...")
print("- Check the full diagnostics dict for any red flags (small bias fit n, very high unspliced frac, etc.).")
print("- For top genes, you may want to manually inspect their phase portraits (U vs S colored by group).")
print("- Consider running the same contrast with mode='advanced' (if cell number permits) and comparing the resulting significant lists.")
print("- Always report the exact mode, n_perm, prior_weight, and whether permutation was used in your methods.")

if __name__ == "__main__":
    print("\nThis template is meant to be read and adapted, not executed directly.")