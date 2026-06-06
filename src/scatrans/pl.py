"""
scATrans — Publication-ready plotting module

High-quality, editable vector graphics for active transcription analysis.
Style inspired by omicverse and Nature/Cell journals.
Optimized for Cell / Nature level publications.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pandas as pd
from contextlib import contextmanager
import logging

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

def set_style(fontfamily='sans-serif', fonts=['Arial', 'Helvetica', 'DejaVu Sans'],
              linewidth=1.2, labelsize=12, titlesize=14, ticksize=10, legendsize=10,
              dpi_preview=150, dpi_save=300, **kwargs):
    """
    Apply a clean, publication-ready matplotlib/seaborn style inspired by omicverse and Nature/Cell journals.
    
    Call this at the beginning of your analysis script for consistent, high-quality figures:
        import scatrans as scat
        scat.pl.set_style(fontsize=12, titlesize=14)  # customize as needed
    
    Key features:
    - Editable text in PDF/SVG (fonttype=42) -- professional standard for journals
    - Clean Arial/Helvetica fonts (vector editable in Illustrator/Inkscape)
    - Minimal spines, professional tick directions (outward)
    - White background, no legend frames
    - High-res output (300 DPI default) ready for *Cell*/*Nature*
    - Fully customizable via parameters or per-plot overrides (figsize, fontsize, dpi, cmap etc.)
    
    Parameters
    ----------
    fontfamily, fonts, linewidth, labelsize, titlesize, ticksize, legendsize : style controls
    dpi_preview, dpi_save : resolution settings
    **kwargs : additional rcParams updates (e.g. axes.grid=False)
    """
    rc_updates = {
        'font.family': fontfamily,
        'font.sans-serif': fonts,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'axes.linewidth': linewidth,
        'axes.edgecolor': '#222222',
        'axes.labelcolor': '#222222',
        'xtick.color': '#222222',
        'ytick.color': '#222222',
        'xtick.major.width': linewidth,
        'ytick.major.width': linewidth,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'legend.frameon': False,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'figure.dpi': dpi_preview,
        'savefig.dpi': dpi_save,
        'savefig.bbox': 'tight',
        'savefig.transparent': False,
        'axes.titlesize': titlesize,
        'axes.labelsize': labelsize,
        'xtick.labelsize': ticksize,
        'ytick.labelsize': ticksize,
        'legend.fontsize': legendsize,
        'figure.titlesize': titlesize + 2,
    }
    rc_updates.update(kwargs)
    mpl.rcParams.update(rc_updates)
    
    sns.set_style("white", {
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": linewidth,
    })
    # Also set whitegrid variant without grid for flexibility
    sns.set_style("whitegrid", {
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": linewidth,
    })


@contextmanager
def style_context(**kwargs):
    """
    Context manager for temporary style application.
    
    Example:
        with scat.pl.style_context(fontsize=12):
            # plots here will use the style
            scat.pl.comet_plot(...)
    """
    original_rc = mpl.rcParams.copy()
    original_sns = sns.axes_style()
    try:
        set_style(**kwargs)
        yield
    finally:
        # Restore original
        mpl.rcParams.update(original_rc)
        sns.set_style(original_sns)


def comet_plot(df, top_n=12, save_path=None,
               title="Active Transcription Drivers",
               point_scale=1.0,
               figsize=(8, 6),
               dpi=300,
               fontsize=12,
               cmap='coolwarm'):
    """
    Generate a highly customizable, publication-quality Comet Plot.

    Parameters
    ----------
    df : pandas.DataFrame
        Results from scat.active_score() containing 'logFC', 'velocity_residual', 'active_score'.
    top_n : int, default 12
        Number of top genes to label with text.
    save_path : str, optional
        Path to save the figure (e.g. "Comet_Plot.pdf").
    title : str, optional
        Plot title.
    point_scale : float, default 1.0
        Global scaling factor for point sizes.
    figsize, dpi, fontsize, cmap : plot aesthetics.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes objects.
    """
    try:
        from adjustText import adjust_text
    except ImportError:
        raise ImportError(
            "The 'adjusttext' package is required for comet_plot. "
            "Please install it with: pip install adjusttext"
        )

    logger.info("🎨 Generating publication-quality Comet Plot...")
    set_style()

    plot_df = df.dropna(subset=['logFC', 'velocity_residual', 'active_score']).copy()
    plot_df = plot_df[plot_df['logFC'] > 0]

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Point sizes scaled by active_score
    sizes = np.clip(plot_df['active_score'] ** 1.6 * 35 * point_scale + 12, 8, 220)

    scatter = ax.scatter(
        x=plot_df['logFC'],
        y=plot_df['velocity_residual'],
        c=plot_df['active_score'],
        s=sizes,
        cmap=cmap,
        alpha=0.85,
        edgecolors='#444444',
        linewidth=0.5,
        zorder=3
    )

    ax.axhline(0, color='#999999', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
    ax.axvline(0, color='#999999', linestyle='--', linewidth=1, alpha=0.5, zorder=1)

    # Label top genes
    top_genes = plot_df.nlargest(top_n, 'active_score')
    texts = []
    for idx, row in top_genes.iterrows():
        txt = ax.text(
            row['logFC'], row['velocity_residual'], f"{idx}",
            fontsize=max(8, fontsize - 2),
            fontweight='bold',
            color='#111111',
            bbox=dict(boxstyle="square,pad=0.1", fc="none", ec="none")
        )
        texts.append(txt)

    if texts:
        adjust_text(
            texts,
            x=plot_df['logFC'].values,
            y=plot_df['velocity_residual'].values,
            arrowprops=dict(arrowstyle='-', color='#666666', lw=0.8, alpha=0.8),
            ax=ax
        )

    ax.set_xlabel('Log2 Fold Change', fontsize=fontsize, fontweight='bold')
    ax.set_ylabel('Bias-corrected Unspliced Residual', fontsize=fontsize, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight='bold', pad=15)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.03, aspect=20)
    cbar.set_label('Active Score', fontsize=max(9, fontsize - 1), fontweight='bold', rotation=270, labelpad=15)
    cbar.outline.set_visible(False)

    sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Comet plot saved to → {save_path}")

    plt.show()
    return fig, ax


def volcano_3d(df, top_n=8, save_path=None, point_scale=1.0,
               title="3D Active Volcano Plot",
               figsize=(10, 8),
               dpi=300,
               fontsize=11,
               cmap='coolwarm'):
    """
    Generate a publication-quality 3D Volcano Plot.

    Parameters
    ----------
    df : pandas.DataFrame
        Results from scat.active_score().
    top_n : int, default 8
        Number of top genes to label.
    save_path : str, optional
        Path to save the figure.

    Returns
    -------
    fig, ax
        Matplotlib figure and 3D axes.
    """
    logger.info("🎨 Generating publication-quality 3D Volcano Plot...")
    set_style()

    plot_df = df.copy().dropna(subset=['logFC', 'p_adj', 'velocity_residual', 'active_score'])
    plot_df['neg_log_pval'] = -np.log10(plot_df['p_adj'].astype(float) + 1e-300)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')

    sizes = np.clip(plot_df['active_score'] ** 1.4 * 18 * point_scale + 8, 10, 180)

    scatter = ax.scatter(
        plot_df['logFC'],
        plot_df['neg_log_pval'],
        plot_df['velocity_residual'],
        c=plot_df['active_score'],
        s=sizes,
        cmap=cmap,
        alpha=0.8,
        edgecolors='#444444',
        linewidth=0.4,
        zorder=3
    )

    # Transparent panes for clean 3D look
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Subtle grid
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"].update({"color": "#E5E5E5", "linestyle": "-"})

    # Label top genes with connecting lines
    top_genes = plot_df.nlargest(top_n, 'active_score')
    x_offset = plot_df['logFC'].max() * 0.1
    z_offset = plot_df['velocity_residual'].max() * 0.15

    for idx, row in top_genes.iterrows():
        px, py, pz = row['logFC'], row['neg_log_pval'], row['velocity_residual']
        tx, ty, tz = px + x_offset, py, pz + z_offset
        ax.plot([px, tx], [py, ty], [pz, tz], color='#888888', ls=':', lw=1.2, alpha=0.8)
        ax.text(tx, ty, tz, f"{idx}", fontsize=max(8, fontsize - 1), fontweight='bold', color='#111111')

    ax.set_xlabel('Log2 Fold Change', fontsize=fontsize, fontweight='bold', labelpad=10)
    ax.set_ylabel('-Log10(adj. P-value)', fontsize=fontsize, fontweight='bold', labelpad=10)
    ax.set_zlabel('Unspliced Residual', fontsize=fontsize, fontweight='bold', labelpad=10)

    if title:
        ax.set_title(title, fontsize=fontsize + 3, fontweight='bold', pad=15)

    ax.view_init(elev=20, azim=-55)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.5, pad=0.1, aspect=15)
    cbar.set_label('Active Score', fontsize=max(9, fontsize - 2), fontweight='bold', rotation=270, labelpad=15)
    cbar.outline.set_visible(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ 3D Volcano plot saved to → {save_path}")

    plt.show()
    return fig, ax


def enrich_dotplot(enrich_df, top_n=15, title="Enrichment Dotplot",
                   save_path=None, figsize=(7, 8), dpi=300, fontsize=12,
                   x="GeneRatio", color_by="Adjusted P-value", size_by="Count",
                   cmap="viridis_r"):
    """
    Generate a highly customizable, clusterProfiler-style dotplot for GO/KEGG/GSEA results.

    Parameters
    ----------
    enrich_df : pandas.DataFrame
        Result from scat.run_enrichment() (or gseapy prerank results).
    top_n : int, default 15
        Number of top enriched terms to show (sorted by significance or x).
    title, save_path, figsize, dpi, fontsize : plot aesthetics.
    x : str, default "GeneRatio"
        Column to use for x-axis. Common choices: "GeneRatio", "FoldEnrichment", "-log10(p.adj)".
    color_by : str, default "Adjusted P-value"
        Column for point color (smaller = more significant usually).
    size_by : str, default "Count"
        Column for point size.
    cmap : str
        Colormap for color_by (use _r versions for p-value so dark = significant).

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """
    if enrich_df.empty:
        logger.warning("⚠️ Enrichment dataframe is empty. Nothing to plot.")
        return None, None

    logger.info("🎨 Generating clusterProfiler-style Enrichment Dotplot...")
    set_style()

    plot_df = enrich_df.head(top_n).copy()
    plot_df = plot_df.iloc[::-1]   # Most significant at top (or reverse if needed)

    # Clean term names
    def clean_term(text):
        text = str(text).split(' (GO:')[0].split(' (KEGG')[0]
        return text[:50] + '...' if len(text) > 50 else text

    plot_df['Term_Clean'] = plot_df['Term'].apply(clean_term)

    # Robust column detection for compatibility with run_enrichment (p.adjust / FDR q-val)
    pval_candidates = ["p.adjust", "Adjusted P-value", "p_adj", "padj", "FDR_qval", "pvalue"]
    pval_col = next((c for c in pval_candidates if c in plot_df.columns), None)
    if pval_col is None:
        pval_col = plot_df.columns[0]  # fallback

    count_candidates = ["Count", "Size", "leadingEdge_count"]
    size_col = next((c for c in count_candidates if c in plot_df.columns), "Count")
    if size_col not in plot_df.columns:
        size_col = plot_df.columns[0]

    # Dynamic x label and data
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

    # Smart default: if GeneRatio has very low variation, switch to FoldEnrichment for better visualization
    if x_col == "GeneRatio" and "FoldEnrichment" in plot_df.columns:
        gene_ratio_range = plot_df["GeneRatio"].max() - plot_df["GeneRatio"].min()
        if gene_ratio_range < 0.08:  # very narrow range
            logger.warning(
                "⚠️ GeneRatio values have very low variation (range < 0.08). "
                "Automatically switching x-axis to 'FoldEnrichment' for better visualization."
            )
            x_col = "FoldEnrichment"
            x_label = "Fold Enrichment"

    # Use robust detection, override with user-specified if present
    if size_by in plot_df.columns:
        size_col = size_by
    if color_by in plot_df.columns:
        color_col = color_by
    else:
        color_col = pval_col if pval_col else plot_df.columns[0]

    # Scale sizes nicely
    sizes = np.clip(plot_df[size_col] * 18 + 30, 20, 400)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    scatter = ax.scatter(
        x=plot_df[x_col],
        y=plot_df['Term_Clean'],
        s=sizes,
        c=plot_df[color_col],
        cmap=cmap,
        edgecolors='#333333',
        linewidth=0.5,
        alpha=0.9
    )

    ax.set_xlabel(x_label, fontsize=fontsize, fontweight='bold', labelpad=10)
    ax.set_ylabel('', fontsize=fontsize)

    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight='bold', pad=20)

    ax.xaxis.grid(True, linestyle='--', color='#DDDDDD', alpha=0.8, zorder=0)
    ax.yaxis.grid(True, linestyle=':', color='#EEEEEE', alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.4, pad=0.03, aspect=15)
    cbar_label = color_col
    if color_col == "Adjusted P-value":
        cbar_label = "Adjusted P-value (smaller = more sig.)"
    cbar.set_label(cbar_label, fontsize=fontsize-1, fontweight='bold', rotation=270, labelpad=20)
    cbar.outline.set_visible(False)

    # Legend for dot size
    try:
        handles, labels = scatter.legend_elements(
            prop="sizes", 
            alpha=0.6, 
            num=4,
            func=lambda s: f"{max(1, int((s - 30) / 18))}"
        )
        ax.legend(handles, labels, title=size_col, loc="center left",
                  bbox_to_anchor=(1.02, 0.5), frameon=False, title_fontsize=fontsize-1)
    except Exception:
        # Fallback if legend_elements fails (e.g. too few unique sizes)
        pass

    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Enrichment Dotplot saved to → {save_path}")

    plt.show()
    return fig, ax


def volcano_plot(df, top_n=10, save_path=None,
                 title="Volcano Plot of Active Transcription",
                 point_scale=1.0,
                 figsize=(8, 6),
                 dpi=300,
                 fontsize=12,
                 cmap='coolwarm',
                 logfc_cutoff=0.5,
                 pval_cutoff=0.05,
                 color_by='active_score'):
    """
    Generate a publication-quality 2D Volcano Plot (logFC vs -log10(adj. P-value)),
    colored by Active Score or significance. Inspired by omicverse and standard
    single-cell DE visualization practices.

    This is the most commonly used volcano style in publications.

    Parameters
    ----------
    df : pandas.DataFrame
        Results from scat.active_score() containing at least 'logFC', 'p_adj',
        and preferably 'active_score' and 'velocity_residual'.
    top_n : int, default 10
        Number of top genes (by active_score) to label.
    save_path : str, optional
        Path to save the figure (e.g. "Volcano_Plot.pdf").
    title : str, optional
        Plot title.
    point_scale : float, default 1.0
        Scaling factor for point sizes.
    figsize, dpi, fontsize, cmap : standard aesthetics.
    logfc_cutoff, pval_cutoff : threshold lines (visual only).
    color_by : str, default 'active_score'
        Column to color points by. If 'active_score' not present, falls back to
        significance categories.

    Returns
    -------
    fig, ax
        Matplotlib figure and axes.
    """
    try:
        from adjustText import adjust_text
    except ImportError:
        raise ImportError(
            "The 'adjusttext' package is required for volcano_plot. "
            "Please install it with: pip install adjusttext"
        )

    logger.info("🎨 Generating publication-quality 2D Volcano Plot...")
    set_style()

    plot_df = df.copy().dropna(subset=['logFC', 'p_adj'])
    plot_df['neg_log_pval'] = -np.log10(plot_df['p_adj'].astype(float) + 1e-300)

    # Determine coloring
    if color_by == 'active_score' and 'active_score' in plot_df.columns:
        color_values = plot_df['active_score']
        cbar_label = 'Active Score'
    else:
        # Fallback: color by significance
        sig_up = (plot_df['logFC'] > logfc_cutoff) & (plot_df['p_adj'] < pval_cutoff)
        color_values = sig_up.astype(int)
        cbar_label = 'Significant (up)'

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    sizes = np.clip(plot_df.get('active_score', 50) ** 1.3 * 8 * point_scale + 15, 12, 180)

    scatter = ax.scatter(
        x=plot_df['logFC'],
        y=plot_df['neg_log_pval'],
        c=color_values,
        s=sizes,
        cmap=cmap,
        alpha=0.75,
        edgecolors='#444444',
        linewidth=0.4,
        zorder=3
    )

    # Threshold lines
    ax.axhline(-np.log10(pval_cutoff), color='#d62728', linestyle='--', linewidth=1.2, alpha=0.8, label=f'p_adj = {pval_cutoff}')
    ax.axvline(logfc_cutoff, color='#d62728', linestyle='--', linewidth=1.2, alpha=0.8, label=f'logFC = {logfc_cutoff}')
    ax.axvline(-logfc_cutoff, color='#1f77b4', linestyle='--', linewidth=1.0, alpha=0.6)

    # Label top genes
    if 'active_score' in plot_df.columns:
        top_genes = plot_df.nlargest(top_n, 'active_score')
    else:
        top_genes = plot_df.nsmallest(top_n, 'p_adj')  # fallback

    texts = []
    for idx, row in top_genes.iterrows():
        txt = ax.text(
            row['logFC'], row['neg_log_pval'], str(idx),
            fontsize=max(8, fontsize - 2),
            fontweight='bold',
            color='#111111',
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7)
        )
        texts.append(txt)

    if texts:
        adjust_text(
            texts,
            x=plot_df['logFC'].values,
            y=plot_df['neg_log_pval'].values,
            arrowprops=dict(arrowstyle='-', color='#888888', lw=0.7, alpha=0.7),
            ax=ax
        )

    ax.set_xlabel('Log2 Fold Change', fontsize=fontsize, fontweight='bold')
    ax.set_ylabel('-Log10(adj. P-value)', fontsize=fontsize, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight='bold', pad=15)

    # Colorbar or legend
    if color_by == 'active_score' and 'active_score' in plot_df.columns:
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02, aspect=20)
        cbar.set_label(cbar_label, fontsize=max(9, fontsize - 1), fontweight='bold', rotation=270, labelpad=15)
        cbar.outline.set_visible(False)
    else:
        ax.legend(loc='upper left', frameon=False, fontsize=fontsize-1)

    sns.despine(ax=ax, top=True, right=True)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ 2D Volcano plot saved to → {save_path}")

    plt.show()
    return fig, ax


# =============================================================================
# Additional Publication-Quality Plotting Functions
# =============================================================================

def bias_diagnostic_plot(results_df, save_path=None,
                         title="Bias Correction Diagnostic",
                         figsize=(12, 5),
                         dpi=300,
                         fontsize=11,
                         show_regression=True):
    """
    Publication-quality diagnostic plot showing the effect of gene length /
    intron number bias correction on velocity delta.

    Highly recommended to demonstrate the unique bias-correction feature of scATrans.

    Parameters
    ----------
    results_df : pd.DataFrame
        Output from scat.active_score() (all_results), must contain
        'velocity_delta_raw', 'velocity_residual', 'gene_length', 'intron_number'.
    save_path, title, figsize, dpi, fontsize : standard
    show_regression : bool
        Whether to overlay the fitted Huber regression line (if available in data).

    Returns
    -------
    fig, axes
    """
    logger.info("🎨 Generating publication-quality Bias Correction Diagnostic Plot...")
    set_style()

    required = ['velocity_delta_raw', 'velocity_residual', 'gene_length', 'intron_number']
    if not all(col in results_df.columns for col in required):
        raise ValueError(f"results_df must contain columns: {required}")

    plot_df = results_df.dropna(subset=required).copy()
    if len(plot_df) < 10:
        logger.warning("⚠️ Too few genes with valid features for diagnostic plot.")
        return None, None

    fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi)

    # Left: Before correction (raw delta vs gene_length)
    ax1 = axes[0]
    x = np.log1p(plot_df['gene_length'])
    y_raw = plot_df['velocity_delta_raw']
    ax1.scatter(x, y_raw, s=15, alpha=0.5, c='#1f77b4', edgecolors='none')
    if show_regression:
        # Simple visual guide line (median trend)
        from scipy.stats import linregress
        try:
            slope, intercept, _, _, _ = linregress(x, y_raw)
            x_line = np.linspace(x.min(), x.max(), 100)
            ax1.plot(x_line, slope * x_line + intercept, '--', color='#d62728', lw=1.5, label='Trend (raw)')
        except Exception:
            pass
    ax1.set_xlabel('log1p(Gene Length)', fontsize=fontsize, fontweight='bold')
    ax1.set_ylabel('Velocity Delta (raw)', fontsize=fontsize, fontweight='bold')
    ax1.set_title('Before Bias Correction', fontsize=fontsize+1, fontweight='bold')
    ax1.legend(frameon=False)
    sns.despine(ax=ax1)

    # Right: After correction (residual vs gene_length)
    ax2 = axes[1]
    y_res = plot_df['velocity_residual']
    ax2.scatter(x, y_res, s=15, alpha=0.5, c='#2ca02c', edgecolors='none')
    ax2.axhline(0, color='#d62728', linestyle='--', lw=1.2, alpha=0.8)
    ax2.set_xlabel('log1p(Gene Length)', fontsize=fontsize, fontweight='bold')
    ax2.set_ylabel('Velocity Residual (bias-corrected)', fontsize=fontsize, fontweight='bold')
    ax2.set_title('After Bias Correction', fontsize=fontsize+1, fontweight='bold')
    sns.despine(ax=ax2)

    fig.suptitle(title, fontsize=fontsize+2, fontweight='bold', y=1.02)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Bias diagnostic plot saved to → {save_path}")

    plt.show()
    return fig, axes


def active_score_rankplot(results_df, top_n=15, save_path=None,
                          title="Top Active Transcription Drivers",
                          figsize=(8, 6),
                          dpi=300,
                          fontsize=11,
                          color_by='logFC',
                          cmap='coolwarm'):
    """
    Clean, publication-quality horizontal lollipop / rank plot of top active genes
    ranked by Active Score. Very effective for highlighting key genes.

    Parameters
    ----------
    results_df : pd.DataFrame
        From scat.active_score().
    top_n : int
        Number of top genes to show.
    save_path, title, etc. : standard
    color_by : str
        Column used for coloring the lollipops (e.g. 'logFC', 'velocity_residual').
    cmap : str
        Colormap for coloring.

    Returns
    -------
    fig, ax
    """
    logger.info("🎨 Generating publication-quality Active Score Rank Plot...")
    set_style()

    plot_df = results_df.nlargest(top_n, 'active_score').copy()
    plot_df = plot_df.iloc[::-1]  # highest at top

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    y_pos = np.arange(len(plot_df))
    colors = plot_df[color_by] if color_by in plot_df.columns else plot_df['active_score']

    # Lollipop style
    ax.hlines(y=y_pos, xmin=0, xmax=plot_df['active_score'], color='#555555', linewidth=1.5, alpha=0.7)
    scatter = ax.scatter(plot_df['active_score'], y_pos,
                         c=colors, cmap=cmap, s=120, edgecolors='#333333', linewidth=0.6, zorder=5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df.index, fontsize=fontsize-1)
    ax.set_xlabel('Active Score', fontsize=fontsize, fontweight='bold')
    ax.set_title(title, fontsize=fontsize+2, fontweight='bold', pad=15)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label(color_by if color_by in plot_df.columns else 'Active Score',
                   fontsize=fontsize-1, fontweight='bold', rotation=270, labelpad=12)
    cbar.outline.set_visible(False)

    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)
    ax.invert_yaxis()  # highest on top
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Active Score Rank Plot saved to → {save_path}")

    plt.show()
    return fig, ax


def active_genes_heatmap(results_df, adata=None, top_n=20, save_path=None,
                         title="Top Active Genes Heatmap",
                         figsize=(10, 8),
                         dpi=300,
                         fontsize=10,
                         groupby=None):
    """
    Publication-quality heatmap of top active genes.

    - If adata + groupby is provided: shows mean expression (or velocity) of top genes
      across groups (recommended for biological insight).
    - Otherwise: shows a clean heatmap of key metrics (active_score, logFC, velocity_residual)
      for the top genes.

    Parameters
    ----------
    results_df : pd.DataFrame
        From scat.active_score().
    adata : AnnData, optional
        If provided with groupby, generates a grouped expression heatmap.
    top_n : int
    save_path, title, figsize, dpi, fontsize : standard
    groupby : str, optional
        obs column for grouping (e.g. "condition").

    Returns
    -------
    fig, ax or ClusterGrid
    """
    logger.info("🎨 Generating publication-quality Active Genes Heatmap...")
    set_style()

    top_genes = results_df.nlargest(top_n, 'active_score').index.tolist()

    if adata is not None and groupby is not None and groupby in adata.obs.columns:
        # Use scanpy-style grouped heatmap (most informative)
        try:
            import scanpy as sc
            adata_sub = adata[:, top_genes].copy()
            sc.tl.dendrogram(adata_sub, groupby=groupby)
            fig = sc.pl.heatmap(adata_sub, var_names=top_genes, groupby=groupby,
                                cmap='viridis', dendrogram=True, show=False,
                                figsize=figsize, title=title)
            if save_path:
                plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
                logger.info(f"✅ Active Genes Heatmap saved to → {save_path}")
            plt.show()
            return fig, None
        except Exception as e:
            logger.warning(f"⚠️ scanpy heatmap failed ({e}). Falling back to metric heatmap.")

    # Fallback: metric heatmap (always works)
    metric_cols = ['active_score', 'logFC', 'velocity_residual']
    available = [c for c in metric_cols if c in results_df.columns]
    if not available:
        logger.warning("⚠️ No suitable columns for heatmap.")
        return None, None

    plot_df = results_df.loc[top_genes, available].copy()
    plot_df = (plot_df - plot_df.min()) / (plot_df.max() - plot_df.min() + 1e-8)  # min-max normalize for visualization

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    sns.heatmap(plot_df, cmap='RdYlBu_r', annot=True, fmt='.2f',
                linewidths=0.5, ax=ax, cbar_kws={'label': 'Normalized Value'})

    ax.set_title(title, fontsize=fontsize+2, fontweight='bold', pad=15)
    ax.set_xlabel('')
    ax.set_ylabel('Top Active Genes', fontsize=fontsize, fontweight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=fontsize-1)
    plt.yticks(fontsize=fontsize-1)
    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Active Genes (metric) Heatmap saved to → {save_path}")

    plt.show()
    return fig, ax


def enrich_barplot(enrich_df, top_n=15, title="Enrichment Barplot",
                   save_path=None, figsize=(7, 8), dpi=300, fontsize=12,
                   x="GeneRatio", color_by="Adjusted P-value", cmap="viridis_r"):
    """
    Generate a clusterProfiler-style horizontal barplot for enrichment results.

    Very useful for quickly showing the most enriched terms with effect size or significance.

    Parameters
    ----------
    enrich_df : pandas.DataFrame
        Result from scat.run_enrichment() (or gseapy prerank results).
    top_n : int, default 15
        Number of top terms to display.
    x : str, default "GeneRatio"
        What to plot on x-axis: "GeneRatio", "FoldEnrichment", "Count", or "-log10(p.adj)".
    color_by : str, default "Adjusted P-value"
        Column used for bar color.
    title, save_path, etc. : standard plot params.

    Returns
    -------
    fig, ax
    """
    if enrich_df.empty:
        logger.warning("⚠️ Enrichment dataframe is empty. Nothing to plot.")
        return None, None

    logger.info("🎨 Generating clusterProfiler-style Enrichment Barplot...")
    set_style()

    plot_df = enrich_df.head(top_n).copy()
    plot_df = plot_df.iloc[::-1]  # top at top after reverse? wait, for barh better sort ascending for bottom-to-top

    # Clean terms
    def clean_term(text):
        text = str(text).split(' (GO:')[0].split(' (KEGG')[0]
        return text[:55] + '...' if len(text) > 55 else text

    plot_df['Term_Clean'] = plot_df['Term'].apply(clean_term)

    # Prepare x data
    if x == "-log10(p.adj)" and "Adjusted P-value" in plot_df.columns:
        plot_df["neg_log"] = -np.log10(plot_df["Adjusted P-value"].clip(lower=1e-300))
        x_col = "neg_log"
        xlabel = "-log10(Adjusted P-value)"
    else:
        x_col = x if x in plot_df.columns else "GeneRatio"
        xlabel = x_col.replace("_", " ").title()

    color_col = color_by if color_by in plot_df.columns else list(plot_df.columns)[0]

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Horizontal bars
    bars = ax.barh(
        y=plot_df['Term_Clean'],
        width=plot_df[x_col],
        color=plot_df[color_col] if color_col in plot_df.columns else '#2E86AB',
        edgecolor='#333333',
        linewidth=0.6,
        alpha=0.85
    )

    # Color the bars by color_col if numeric
    if pd.api.types.is_numeric_dtype(plot_df[color_col]):
        # Use a colormap
        norm = plt.Normalize(plot_df[color_col].min(), plot_df[color_col].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        for bar, val in zip(bars, plot_df[color_col]):
            bar.set_color(sm.to_rgba(val))
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label(color_col, fontsize=fontsize-1, fontweight='bold', rotation=270, labelpad=15)
        cbar.outline.set_visible(False)

    ax.set_xlabel(xlabel, fontsize=fontsize, fontweight='bold', labelpad=8)
    ax.set_ylabel('')
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight='bold', pad=15)

    ax.xaxis.grid(True, linestyle='--', color='#CCCCCC', alpha=0.7)
    ax.set_axisbelow(True)

    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Enrichment Barplot saved to → {save_path}")

    plt.show()
    return fig, ax
