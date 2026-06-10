"""
scATrans plotting module.

Utilities for generating clear, publication-suitable figures for active
transcription analysis results. The default style aims for clean vector
output (PDF/SVG), minimal non-data ink, and readable labels suitable for
scientific journals.

The style configuration is inspired by high-quality single-cell visualization
practices (e.g., libraries such as OmicVerse and Scanpy extensions).
"""

import logging
from contextlib import contextmanager
from typing import Any, Iterable, List, Optional, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def set_style(
    fontfamily="sans-serif",
    fonts=["Arial", "Helvetica", "DejaVu Sans"],
    linewidth=1.0,
    labelsize=11,
    titlesize=12,
    ticksize=9,
    legendsize=9,
    dpi_preview=150,
    dpi_save=300,
    **kwargs,
):
    """
    Apply a clean, minimal style suitable for scientific publication figures.

    Key characteristics (inspired by professional single-cell visualization
    libraries such as OmicVerse):
    - Vector-friendly output (Type 42 fonts for easy editing of PDF/SVG in
      Illustrator, Affinity, etc.)
    - Minimal non-data ink (no top/right spines)
    - Consistent, readable sizes for journal figures (typically 11 pt labels)
    - White background, high contrast, no unnecessary grids

    It is recommended to call this once near the beginning of an analysis
    script or notebook. All `scat.pl.*` plotting functions respect these
    settings.

    Parameters
    ----------
    fontfamily, fonts
        Control the font stack. Arial/Helvetica are preferred for journals.
    linewidth
        Base width for axes, ticks and spines.
    labelsize, titlesize, ticksize, legendsize
        Font sizes in points.
    dpi_preview, dpi_save
        DPI used for on-screen display vs. saved files.
    **kwargs
        Any additional matplotlib rcParams (will override the defaults).
    """
    rc_updates = {
        "font.family": fontfamily,
        "font.sans-serif": fonts,
        # Critical for publication: editable vector text
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        # Clean, high-contrast, minimal style
        "axes.linewidth": linewidth,
        "axes.edgecolor": "#1f1f1f",
        "axes.labelcolor": "#1f1f1f",
        "xtick.color": "#1f1f1f",
        "ytick.color": "#1f1f1f",
        "xtick.major.width": linewidth * 0.8,
        "ytick.major.width": linewidth * 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "figure.dpi": dpi_preview,
        "savefig.dpi": dpi_save,
        "savefig.bbox": "tight",
        "savefig.transparent": False,
        "axes.titlesize": titlesize,
        "axes.labelsize": labelsize,
        "xtick.labelsize": ticksize,
        "ytick.labelsize": ticksize,
        "legend.fontsize": legendsize,
        "figure.titlesize": titlesize + 1,
        # Reduce visual clutter (consistent with high-end scRNA-seq figures)
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
    rc_updates.update(kwargs)
    mpl.rcParams.update(rc_updates)

    sns.set_style(
        "white",
        {
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": linewidth,
        },
    )


@contextmanager
def style_context(**kwargs):
    """
    Context manager for temporary style application.
    """
    original_rc = mpl.rcParams.copy()
    original_sns = sns.axes_style()
    try:
        set_style(**kwargs)
        yield
    finally:
        mpl.rcParams.update(original_rc)
        sns.set_style(original_sns)


def comet_plot(
    df,
    top_n=12,
    save_path=None,
    title="Active Transcription Drivers",
    point_scale=1.0,
    figsize=(8, 6),
    dpi=300,
    fontsize=12,
    cmap="coolwarm",
    ax=None,
):
    """
    Comet plot of log fold change vs. bias-corrected unspliced residual.

    Point size and color are mapped to the active score. Designed to produce
    clear figures suitable for scientific publications with minimal further
    editing.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        If provided, plot into this axes instead of creating a new figure.
        Useful for embedding in multi-panel publication figures.
    """
    try:
        from adjustText import adjust_text
    except ImportError:
        raise ImportError(
            "The 'adjusttext' package is required for comet_plot. "
            "Please install it with: pip install adjusttext"
        )

    logger.info("Generating comet plot...")
    set_style()

    plot_df = df.dropna(subset=["logFC", "velocity_residual", "active_score"]).copy()
    plot_df = plot_df[plot_df["logFC"] > 0]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    sizes = np.clip(plot_df["active_score"] ** 1.6 * 35 * point_scale + 12, 8, 220)

    scatter = ax.scatter(
        x=plot_df["logFC"],
        y=plot_df["velocity_residual"],
        c=plot_df["active_score"],
        s=sizes,
        cmap=cmap,
        alpha=0.85,
        edgecolors="#444444",
        linewidth=0.5,
        zorder=3,
    )

    ax.axhline(0, color="#999999", linestyle="--", linewidth=1, alpha=0.5, zorder=1)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=1, alpha=0.5, zorder=1)

    top_genes = plot_df.nlargest(top_n, "active_score")
    texts = []
    for idx, row in top_genes.iterrows():
        txt = ax.text(
            row["logFC"],
            row["velocity_residual"],
            f"{idx}",
            fontsize=max(8, fontsize - 2),
            fontweight="bold",
            color="#111111",
            bbox=dict(boxstyle="square,pad=0.1", fc="none", ec="none"),
        )
        texts.append(txt)

    if texts:
        adjust_text(
            texts,
            x=plot_df["logFC"].values,
            y=plot_df["velocity_residual"].values,
            arrowprops=dict(arrowstyle="-", color="#666666", lw=0.8, alpha=0.8),
            ax=ax,
        )

    ax.set_xlabel("Log2 Fold Change", fontsize=fontsize, fontweight="bold")
    ax.set_ylabel("Bias-corrected Unspliced Residual", fontsize=fontsize, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=15)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.03, aspect=20)
    cbar.set_label(
        "Active Score", fontsize=max(9, fontsize - 1), fontweight="bold", rotation=270, labelpad=15
    )
    cbar.outline.set_visible(False)

    sns.despine(ax=ax, top=True, right=True)

    if _created_fig:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=True)
        logger.info("Comet plot saved to %s", save_path)

    if _created_fig:
        plt.show()
    return fig, ax


def volcano_3d(
    df,
    top_n=8,
    save_path=None,
    point_scale=1.0,
    title="3D Active Volcano Plot",
    figsize=(10, 8),
    dpi=300,
    fontsize=11,
    cmap="coolwarm",
    ax=None,
):
    """
    3D volcano-style view (logFC, -log10(p_adj), velocity residual).

    If `ax` (a 3D axes) is provided, plot into it.
    """
    logger.info("Generating 3D volcano plot...")
    set_style()

    plot_df = df.copy().dropna(subset=["logFC", "p_adj", "velocity_residual", "active_score"])
    plot_df["neg_log_pval"] = -np.log10(plot_df["p_adj"].astype(float) + 1e-300)

    if ax is None:
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection="3d")
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    sizes = np.clip(plot_df["active_score"] ** 1.4 * 18 * point_scale + 8, 10, 180)

    scatter = ax.scatter(
        plot_df["logFC"],
        plot_df["neg_log_pval"],
        plot_df["velocity_residual"],
        c=plot_df["active_score"],
        s=sizes,
        cmap=cmap,
        alpha=0.8,
        edgecolors="#444444",
        linewidth=0.4,
        zorder=3,
    )

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis.line.set_color((1.0, 1.0, 1.0, 0.0))

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"].update({"color": "#E5E5E5", "linestyle": "-"})

    top_genes = plot_df.nlargest(top_n, "active_score")
    x_offset = plot_df["logFC"].max() * 0.1
    z_offset = plot_df["velocity_residual"].max() * 0.15

    for idx, row in top_genes.iterrows():
        px, py, pz = row["logFC"], row["neg_log_pval"], row["velocity_residual"]
        tx, ty, tz = px + x_offset, py, pz + z_offset
        ax.plot([px, tx], [py, ty], [pz, tz], color="#888888", ls=":", lw=1.2, alpha=0.8)
        ax.text(
            tx, ty, tz, f"{idx}", fontsize=max(8, fontsize - 1), fontweight="bold", color="#111111"
        )

    ax.set_xlabel("Log2 Fold Change", fontsize=fontsize, fontweight="bold", labelpad=10)
    ax.set_ylabel("-Log10(adj. P-value)", fontsize=fontsize, fontweight="bold", labelpad=10)
    ax.set_zlabel("Unspliced Residual", fontsize=fontsize, fontweight="bold", labelpad=10)

    if title:
        ax.set_title(title, fontsize=fontsize + 3, fontweight="bold", pad=15)

    ax.view_init(elev=20, azim=-55)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, pad=0.1, aspect=15)
    cbar.set_label(
        "Active Score", fontsize=max(9, fontsize - 2), fontweight="bold", rotation=270, labelpad=15
    )
    cbar.outline.set_visible(False)

    if _created_fig:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=True)
        logger.info("3D Volcano plot saved to %s", save_path)

    if _created_fig:
        plt.show()
    return fig, ax


def enrich_dotplot(
    enrich_df,
    top_n=15,
    show_terms: Optional[Union[int, List[str], Tuple[str, ...]]] = None,
    title="Enrichment Dotplot",
    save_path=None,
    figsize=(7, 8),
    dpi=300,
    fontsize=12,
    x="GeneRatio",
    color_by="Adjusted P-value",
    size_by="Count",
    cmap="viridis_r",
    ax=None,
):
    """
    Dotplot for enrichment results (clusterProfiler style).

    `show_terms` gives clusterProfiler-like flexibility:
      - int: show top N terms (overrides top_n)
      - list/tuple of str: show exactly the matching terms (match on Term or Description;
        order of the list is respected when possible). This is analogous to
        `dotplot(..., showCategory = c("term1", "term2"))`.

    `top_n` is still supported for the common "top N" case (when show_terms is None).
    Supports `ax` for embedding in publication multi-panel figures.
    """
    if enrich_df.empty:
        logger.warning("Enrichment dataframe is empty. Nothing to plot.")
        return None, None

    logger.info("Generating enrichment dotplot...")
    set_style()

    # clusterProfiler-style term selection (show_terms takes precedence)
    if show_terms is not None:
        if isinstance(show_terms, int):
            plot_df = enrich_df.head(show_terms).copy()
        else:
            wanted = {str(x).strip().lower() for x in show_terms}
            def _matches(row):
                t = str(row.get("Term", "")).strip().lower()
                d = str(row.get("Description", "")).strip().lower()
                for w in wanted:
                    if w in t or w in d:
                        return True
                return False
            mask = enrich_df.apply(_matches, axis=1)
            plot_df = enrich_df[mask].copy()
            # Try to preserve caller-specified order
            if not plot_df.empty and len(show_terms) > 0:
                order_map = {str(x).strip().lower(): i for i, x in enumerate(show_terms)}
                def _order_key(row):
                    t = str(row.get("Term", "")).strip().lower()
                    d = str(row.get("Description", "")).strip().lower()
                    return min(order_map.get(t, 10**9), order_map.get(d, 10**9))
                plot_df = plot_df.copy()
                plot_df["_sel_order"] = plot_df.apply(_order_key, axis=1)
                plot_df = plot_df.sort_values("_sel_order").drop(columns=["_sel_order"], errors="ignore")
    else:
        plot_df = enrich_df.head(top_n).copy()

    plot_df = plot_df.iloc[::-1]  # visual: top at top of y axis for horizontal dotplot

    def clean_term(text):
        text = str(text).split(" (GO:")[0].split(" (KEGG")[0]
        return text[:50] + "..." if len(text) > 50 else text

    plot_df["Term_Clean"] = plot_df["Term"].apply(clean_term)

    pval_candidates = ["p.adjust", "Adjusted P-value", "p_adj", "padj", "FDR_qval", "pvalue"]
    pval_col = next((c for c in pval_candidates if c in plot_df.columns), None)
    if pval_col is None:
        pval_col = plot_df.columns[0]

    count_candidates = ["Count", "Size", "leadingEdge_count"]
    size_col = next((c for c in count_candidates if c in plot_df.columns), "Count")

    if x == "-log10(p.adj)" or x == "-log10(p.adjust)":
        if pval_col:
            plot_df["neg_log_padj"] = -np.log10(plot_df[pval_col].astype(float).clip(lower=1e-300))
            x_col = "neg_log_padj"
            x_label = f"-log10({pval_col})"
        else:
            x_col = x if x in plot_df.columns else "GeneRatio"
            x_label = x_col
    else:
        x_col = x if x in plot_df.columns else "GeneRatio"
        x_label = x_col

    if x_col == "GeneRatio" and "FoldEnrichment" in plot_df.columns:
        gene_ratio_range = plot_df["GeneRatio"].max() - plot_df["GeneRatio"].min()
        if gene_ratio_range < 0.08:
            logger.warning(
                "⚠️ GeneRatio values have very low variation. Switching to 'FoldEnrichment'."
            )
            x_col = "FoldEnrichment"
            x_label = "Fold Enrichment"

    if size_by in plot_df.columns:
        size_col = size_by
    if color_by in plot_df.columns:
        color_col = color_by
    else:
        color_col = pval_col if pval_col else plot_df.columns[0]

    sizes = np.clip(plot_df[size_col] * 18 + 30, 20, 400)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    scatter = ax.scatter(
        x=plot_df[x_col],
        y=plot_df["Term_Clean"],
        s=sizes,
        c=plot_df[color_col],
        cmap=cmap,
        edgecolors="#333333",
        linewidth=0.5,
        alpha=0.9,
    )

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight="bold", labelpad=10)
    ax.set_ylabel("", fontsize=fontsize)

    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=20)

    ax.xaxis.grid(True, linestyle="--", color="#DDDDDD", alpha=0.8, zorder=0)
    ax.yaxis.grid(True, linestyle=":", color="#EEEEEE", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.4, pad=0.03, aspect=15)
    cbar_label = color_col
    if color_col == "Adjusted P-value":
        cbar_label = "Adjusted P-value (smaller = more sig.)"
    cbar.set_label(cbar_label, fontsize=fontsize - 1, fontweight="bold", rotation=270, labelpad=20)
    cbar.outline.set_visible(False)

    try:
        handles, labels = scatter.legend_elements(
            prop="sizes", alpha=0.6, num=4, func=lambda s: f"{max(1, int((s - 30) / 18))}"
        )
        ax.legend(
            handles,
            labels,
            title=size_col,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            title_fontsize=fontsize - 1,
        )
    except Exception:
        pass

    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)

    if _created_fig:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=True)
        logger.info("Enrichment dotplot saved to %s", save_path)

    if _created_fig:
        plt.show()
    return fig, ax


def volcano_plot(
    df,
    top_n=10,
    label_genes: Optional[Iterable[str]] = None,
    save_path=None,
    title="Volcano Plot of Active Transcription",
    point_scale=1.0,
    figsize=(8, 6),
    dpi=300,
    fontsize=12,
    cmap="coolwarm",
    logfc_cutoff=0.5,
    pval_cutoff=0.05,
    color_by="active_score",
    ax=None,
):
    """
    2D volcano plot with ggVolcano-inspired flexibility and style options.

    - `top_n`: number of top genes (by active_score or p_adj) to label.
    - `label_genes`: iterable of gene names to *force* label (manual specification,
      even if not in top_n). This + top_n gives the common ggVolcano usage
      pattern (label_number + explicit genes). Duplicates are handled automatically.
    - When color_by != active_score (or not present), falls back to classic
      up / down / ns coloring (red / blue / gray) based on cutoffs — matching
      the popular ggVolcano "beautiful volcano" look.
    - Cutoff lines are drawn with labels.
    - Supports `ax` for embedding.

    Style reference: https://github.com/BioSenior/ggVolcano (label control,
    clean up/down distinction, readable labels with repel).
    """
    try:
        from adjustText import adjust_text
    except ImportError:
        raise ImportError(
            "The 'adjusttext' package is required for volcano_plot. "
            "Please install it with: pip install adjusttext"
        )

    logger.info("Generating 2D volcano plot...")
    set_style()

    plot_df = df.copy().dropna(subset=["logFC", "p_adj"])
    plot_df["neg_log_pval"] = -np.log10(plot_df["p_adj"].astype(float) + 1e-300)

    # ggVolcano-style classic coloring (up/down/ns) when not using active_score
    use_classic = (color_by != "active_score") or ("active_score" not in plot_df.columns)
    if color_by == "active_score" and "active_score" in plot_df.columns:
        color_values = plot_df["active_score"]
        cbar_label = "Active Score"
        colors_for_scatter = None
    else:
        up_mask = (plot_df["logFC"] > logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
        down_mask = (plot_df["logFC"] < -logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
        color_values = np.where(up_mask, 2, np.where(down_mask, 1, 0))  # 2=up, 1=down, 0=ns
        cbar_label = None
        # Use explicit nice colors similar to common ggVolcano / EnhancedVolcano
        colors_for_scatter = ["#808080", "#1f77b4", "#d62728"]  # ns, down, up

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    sizes = np.clip(plot_df.get("active_score", 50) ** 1.3 * 8 * point_scale + 15, 12, 180)

    scatter_kwargs = dict(
        x=plot_df["logFC"],
        y=plot_df["neg_log_pval"],
        s=sizes,
        alpha=0.75,
        edgecolors="#444444",
        linewidth=0.4,
        zorder=3,
    )
    if colors_for_scatter is not None:
        scatter = ax.scatter(c=color_values, cmap=None, color=[colors_for_scatter[int(c)] for c in color_values], **scatter_kwargs)
    else:
        scatter = ax.scatter(c=color_values, cmap=cmap, **scatter_kwargs)

    ax.axhline(
        -np.log10(pval_cutoff),
        color="#d62728",
        linestyle="--",
        linewidth=1.2,
        alpha=0.85,
        label=f"p_adj = {pval_cutoff}",
    )
    ax.axvline(
        logfc_cutoff,
        color="#d62728",
        linestyle="--",
        linewidth=1.2,
        alpha=0.85,
        label=f"logFC = {logfc_cutoff}",
    )
    ax.axvline(-logfc_cutoff, color="#1f77b4", linestyle="--", linewidth=1.0, alpha=0.7)

    # --- ggVolcano-like gene labeling: top_n + manually specified genes ---
    genes_to_label = set()
    if label_genes is not None:
        genes_to_label.update(str(g).strip() for g in label_genes if str(g).strip())

    # Always include the top_n (by active_score when available)
    if "active_score" in plot_df.columns:
        top_df = plot_df.nlargest(top_n, "active_score")
    else:
        top_df = plot_df.nsmallest(top_n, "p_adj")
    for g in top_df.index:
        genes_to_label.add(str(g))

    label_df = plot_df.loc[plot_df.index.astype(str).isin(genes_to_label)].copy() if genes_to_label else pd.DataFrame()

    texts = []
    for idx, row in label_df.iterrows():
        txt = ax.text(
            row["logFC"],
            row["neg_log_pval"],
            str(idx),
            fontsize=max(8, fontsize - 2),
            fontweight="bold",
            color="#111111",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.75),
        )
        texts.append(txt)

    if texts:
        adjust_text(
            texts,
            x=plot_df["logFC"].values,
            y=plot_df["neg_log_pval"].values,
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.7, alpha=0.7),
            ax=ax,
        )

    ax.set_xlabel("Log2 Fold Change", fontsize=fontsize, fontweight="bold")
    ax.set_ylabel("-Log10(adj. P-value)", fontsize=fontsize, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=15)

    if color_by == "active_score" and "active_score" in plot_df.columns:
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02, aspect=20)
        cbar.set_label(
            cbar_label, fontsize=max(9, fontsize - 1), fontweight="bold", rotation=270, labelpad=15
        )
        cbar.outline.set_visible(False)
    else:
        ax.legend(loc="upper left", frameon=False, fontsize=fontsize - 1)

    sns.despine(ax=ax, top=True, right=True)

    if _created_fig:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=True)
        logger.info("2D Volcano plot saved to %s", save_path)

    if _created_fig:
        plt.show()
    return fig, ax


def bias_diagnostic_plot(
    results_df,
    save_path=None,
    title="Bias Correction Diagnostic",
    figsize=(12, 5),
    dpi=300,
    fontsize=11,
    show_regression=True,
    axes=None,
):
    """
    Diagnostic plot showing the effect of gene length / intron number bias
    correction on velocity delta (before vs after).

    Supports external `axes` (tuple of two Axes) for embedding in custom figures.
    """
    logger.info("Generating bias correction diagnostic plot...")
    set_style()

    required = ["velocity_delta_raw", "velocity_residual", "gene_length", "intron_number"]
    if not all(col in results_df.columns for col in required):
        raise ValueError(f"results_df must contain columns: {required}")

    plot_df = results_df.dropna(subset=required).copy()
    if len(plot_df) < 10:
        logger.warning("Too few genes with valid features for diagnostic plot.")
        return None, None

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        _created_fig = True
        ax1 = axes[0]
        ax2 = axes[1]
    else:
        if len(axes) != 2:
            raise ValueError("axes must be a sequence of exactly two matplotlib Axes")
        fig = axes[0].figure
        _created_fig = False
        ax1 = axes[0]
        ax2 = axes[1]

    # Left: Before correction
    x = np.log1p(plot_df["gene_length"])
    y_raw = plot_df["velocity_delta_raw"]
    ax1.scatter(x, y_raw, s=15, alpha=0.5, c="#1f77b4", edgecolors="none")
    if show_regression:
        from scipy.stats import linregress

        try:
            slope, intercept, _, _, _ = linregress(x, y_raw)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax1.plot(
                x_line,
                slope * x_line + intercept,
                "--",
                color="#d62728",
                lw=1.5,
                label="Trend (raw)",
            )
        except Exception:
            pass
    ax1.set_xlabel("log1p(Gene Length)", fontsize=fontsize, fontweight="bold")
    ax1.set_ylabel("Velocity Delta (raw)", fontsize=fontsize, fontweight="bold")
    ax1.set_title("Before Bias Correction", fontsize=fontsize + 1, fontweight="bold")
    ax1.legend(frameon=False)
    sns.despine(ax=ax1)

    # Right: After correction
    y_res = plot_df["velocity_residual"]
    ax2.scatter(x, y_res, s=15, alpha=0.5, c="#2ca02c", edgecolors="none")
    ax2.axhline(0, color="#d62728", linestyle="--", lw=1.2, alpha=0.8)
    ax2.set_xlabel("log1p(Gene Length)", fontsize=fontsize, fontweight="bold")
    ax2.set_ylabel("Velocity Residual (bias-corrected)", fontsize=fontsize, fontweight="bold")
    ax2.set_title("After Bias Correction", fontsize=fontsize + 1, fontweight="bold")
    sns.despine(ax=ax2)

    if title:
        fig.suptitle(title, fontsize=fontsize + 2, fontweight="bold", y=1.02)

    if _created_fig:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=dpi, bbox_inches="tight", transparent=True)
        logger.info("Bias diagnostic plot saved to %s", save_path)

    if _created_fig:
        plt.show()
    return fig, axes


# =============================================================================
# Additional / legacy plotting helpers
# =============================================================================


def enrich_barplot(enrich_df, top_n=15, title="Enrichment Barplot", save_path=None, **kwargs):
    """Barplot wrapper around the dotplot implementation (for API compatibility)."""
    return enrich_dotplot(enrich_df, top_n=top_n, title=title, save_path=save_path, **kwargs)


def active_score_rankplot(results_df, top_n=20, save_path=None, ax=None, **kwargs):
    """
    Simple horizontal rank barplot of top active scores.

    Supports `ax` for embedding. For publication figures prefer
    `pl.comet_plot` or `pl.volcano_plot`.
    """
    logger.info("Generating active score rank plot...")
    set_style()

    if results_df is None or results_df.empty:
        logger.warning("No results to plot.")
        return None, None

    plot_df = results_df.head(top_n).copy()
    plot_df = plot_df.iloc[::-1]  # top at top

    import seaborn as sns

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, max(4, 0.35 * top_n)), dpi=300)
        _created = True
    else:
        fig = ax.figure
        _created = False

    sns.barplot(
        data=plot_df,
        y=plot_df.index,
        x="active_score",
        ax=ax,
        color="#2ca02c",
        edgecolor="#333333",
        linewidth=0.5,
    )
    ax.set_xlabel("Active Score", fontweight="bold")
    ax.set_ylabel("")
    ax.set_title("Top Active Drivers (rank)", fontweight="bold", pad=10)
    sns.despine(ax=ax)

    if _created:
        plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
        logger.info("Rank plot saved to %s", save_path)

    if _created:
        plt.show()
    return fig, ax


def active_genes_heatmap(adata, genes=None, groupby=None, save_path=None, **kwargs):
    """
    Convenience wrapper around scanpy heatmap for the active driver genes.

    Users are encouraged to call scanpy.pl.heatmap directly with the genes
    returned by active_score for full control.
    """
    if genes is None:
        # try to guess from var
        if "active_score" in adata.var.columns:
            genes = adata.var.nlargest(20, "active_score").index.tolist()
        else:
            logger.warning("No genes provided and no active_score column found.")
            return None, None

    logger.info(
        "active_genes_heatmap: delegating to scanpy.pl.heatmap (recommended for full control)"
    )
    try:
        sc = __import__("scanpy", fromlist=["pl"])
        fig = sc.pl.heatmap(
            adata,
            var_names=genes,
            groupby=groupby,
            show=False,
            save=None,
            **kwargs,
        )
        if save_path:
            import matplotlib.pyplot as plt

            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            logger.info("Heatmap saved → %s", save_path)
        # scanpy usually shows; we return what we can
        return fig, None
    except Exception as e:
        logger.warning("active_genes_heatmap could not render via scanpy: %s", e)
        return None, None


def set_nature_style():
    """Legacy alias for set_style() kept for backward compatibility."""
    set_style()


def velocity_phase_portraits(
    adata,
    genes,
    groupby=None,
    spliced_layer="spliced",
    unspliced_layer="unspliced",
    max_genes=6,
    figsize_per_gene=(2.8, 2.4),
    save_path=None,
    **kwargs,
):
    """
    Quick diagnostic grid of unspliced vs spliced (phase-portrait style) for selected genes.

    Useful for visually inspecting whether top active genes show the expected excess
    nascent RNA in the target group. Points are colored by the groupby column when provided.

    This is intentionally lightweight — for full control users are encouraged to write
    their own small U/S scatter functions.

    Parameters
    ----------
    genes : list-like
        Gene names (index of adata.var) to plot.
    groupby : str, optional
        obs column used for coloring (e.g. the same contrast column used in active_score).
    max_genes : int
        Maximum number of genes to plot (grid will be truncated).
    """
    import math

    genes = list(genes)[:max_genes]
    if not genes:
        logger.warning("No genes provided for phase portraits.")
        return None, None

    set_style()

    n = len(genes)
    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(ncols * figsize_per_gene[0], nrows * figsize_per_gene[1]),
        dpi=150,
        squeeze=False,
    )
    axes = axes.flatten()

    for i, g in enumerate(genes):
        ax = axes[i]
        if g not in adata.var_names:
            ax.text(0.5, 0.5, f"{g}\n(not found)", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue

        gidx = adata.var_names.get_loc(g)
        u = adata.layers[unspliced_layer][:, gidx]
        s = adata.layers[spliced_layer][:, gidx]
        if sparse.issparse(u):
            u = u.toarray().ravel()
            s = s.toarray().ravel()

        color = None
        if groupby and groupby in adata.obs:
            # Convert labels to numeric codes for scatter c= (avoids matplotlib error on string arrays)
            try:
                color = pd.Categorical(adata.obs[groupby]).codes
            except Exception:
                color = None

        cval = color if (color is not None and np.asarray(color).size == len(s)) else "#2ca02c"
        ax.scatter(s, u, c=cval, s=8, alpha=0.6, edgecolors="none")
        ax.set_xlabel("Spliced", fontsize=9)
        ax.set_ylabel("Unspliced", fontsize=9)
        ax.set_title(str(g), fontsize=10, fontweight="bold")
        sns.despine(ax=ax)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info("Phase portraits saved to %s", save_path)

    plt.show()
    return fig, axes
