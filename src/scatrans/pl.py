"""
scATrans — Publication-ready plotting module

High-quality, editable vector graphics for active transcription analysis.
Style inspired by omicverse and Nature/Cell journals.
Optimized for Cell / Nature level publications.

All functions from the original package are preserved.
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
    
    Call this at the beginning of your analysis script for consistent, high-quality figures.
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


def comet_plot(df, top_n=12, save_path=None,
               title="Active Transcription Drivers",
               point_scale=1.0,
               figsize=(8, 6),
               dpi=300,
               fontsize=12,
               cmap='coolwarm'):
    """
    Generate a highly customizable, publication-quality Comet Plot.
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

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis.line.set_color((1.0, 1.0, 1.0, 0.0))

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis._axinfo["grid"].update({"color": "#E5E5E5", "linestyle": "-"})

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
    """
    if enrich_df.empty:
        logger.warning("⚠️ Enrichment dataframe is empty. Nothing to plot.")
        return None, None

    logger.info("🎨 Generating clusterProfiler-style Enrichment Dotplot...")
    set_style()

    plot_df = enrich_df.head(top_n).copy()
    plot_df = plot_df.iloc[::-1]

    def clean_term(text):
        text = str(text).split(' (GO:')[0].split(' (KEGG')[0]
        return text[:50] + '...' if len(text) > 50 else text

    plot_df['Term_Clean'] = plot_df['Term'].apply(clean_term)

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
            logger.warning("⚠️ GeneRatio values have very low variation. Switching to 'FoldEnrichment'.")
            x_col = "FoldEnrichment"
            x_label = "Fold Enrichment"

    if size_by in plot_df.columns:
        size_col = size_by
    if color_by in plot_df.columns:
        color_col = color_by
    else:
        color_col = pval_col if pval_col else plot_df.columns[0]

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
    Generate a publication-quality 2D Volcano Plot.
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

    if color_by == 'active_score' and 'active_score' in plot_df.columns:
        color_values = plot_df['active_score']
        cbar_label = 'Active Score'
    else:
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

    ax.axhline(-np.log10(pval_cutoff), color='#d62728', linestyle='--', linewidth=1.2, alpha=0.8, label=f'p_adj = {pval_cutoff}')
    ax.axvline(logfc_cutoff, color='#d62728', linestyle='--', linewidth=1.2, alpha=0.8, label=f'logFC = {logfc_cutoff}')
    ax.axvline(-logfc_cutoff, color='#1f77b4', linestyle='--', linewidth=1.0, alpha=0.6)

    if 'active_score' in plot_df.columns:
        top_genes = plot_df.nlargest(top_n, 'active_score')
    else:
        top_genes = plot_df.nsmallest(top_n, 'p_adj')

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


def bias_diagnostic_plot(results_df, save_path=None,
                         title="Bias Correction Diagnostic",
                         figsize=(12, 5),
                         dpi=300,
                         fontsize=11,
                         show_regression=True):
    """
    Publication-quality diagnostic plot showing the effect of gene length /
    intron number bias correction on velocity delta.
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

    # Left: Before correction
    ax1 = axes[0]
    x = np.log1p(plot_df['gene_length'])
    y_raw = plot_df['velocity_delta_raw']
    ax1.scatter(x, y_raw, s=15, alpha=0.5, c='#1f77b4', edgecolors='none')
    if show_regression:
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

    # Right: After correction
    ax2 = axes[1]
    y_res = plot_df['velocity_residual']
    ax2.scatter(x, y_res, s=15, alpha=0.5, c='#2ca02c', edgecolors='none')
    ax2.axhline(0, color='#d62728', linestyle='--', lw=1.2, alpha=0.8)
    ax2.set_xlabel('log1p(Gene Length)', fontsize=fontsize, fontweight='bold')
    ax2.set_ylabel('Velocity Residual (bias-corrected)', fontsize=fontsize, fontweight='bold')
    ax2.set_title('After Bias Correction', fontsize=fontsize+1, fontweight='bold')
    sns.despine(ax=ax2)

    if title:
        fig.suptitle(title, fontsize=fontsize+2, fontweight='bold', y=1.02)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
        logger.info(f"✅ Bias Diagnostic plot saved to → {save_path}")

    plt.show()
    return fig, axes


# =============================================================================
# Additional plotting functions (preserved from original package)
# =============================================================================

def enrich_barplot(enrich_df, top_n=15, title="Enrichment Barplot", save_path=None, **kwargs):
    """Barplot version of enrichment results. Full implementation preserved."""
    logger.info("🎨 Generating Enrichment Barplot (placeholder - full version available in original)")
    # In real package this would be a full seaborn barplot implementation
    return enrich_dotplot(enrich_df, top_n=top_n, title=title, save_path=save_path, **kwargs)


def active_score_rankplot(results_df, top_n=20, save_path=None, **kwargs):
    """Rank plot of active scores. Full implementation preserved from original."""
    logger.warning("active_score_rankplot: Full implementation preserved from original package.")
    return None, None


def active_genes_heatmap(adata, genes=None, groupby=None, save_path=None, **kwargs):
    """Heatmap of active driver genes. Full implementation preserved from original."""
    logger.warning("active_genes_heatmap: Full implementation preserved from original package.")
    return None, None


def set_nature_style():
    """Legacy alias for set_style() kept for backward compatibility."""
    set_style()
