"""
scATrans plotting module.

Utilities for generating clear, publication-suitable figures for active
transcription analysis results (volcano/comet, enrichment dotplots, rank plots,
bias diagnostics, phase portraits, heatmaps, etc.). The default style aims for
clean vector output (PDF/SVG), minimal non-data ink, and readable labels
suitable for scientific journals.

Design and implementation details draw heavily from high-quality patterns in
OmicVerse (https://github.com/omicverse/omicverse) and gseapy, including:
- constrained_layout + careful bbox handling for colorbars + legends (no more
  overlapping "two legend parts")
- direct `s=` fixed-size controls + min/max + diagnostics for dense plots
- gradient coloring and clean ranked barplots
- outward spine offsets, patch legends, balanced gene labeling
- consistent ax= embedding, return (fig, ax), and save_path behavior
- modern adjustText usage and sensible defaults for journal figures

**Figure lifecycle (important for batch use):**
Plot functions that create figures return them. By default we do **not** close them
(`close=False`) so notebooks and interactive use keep the figure. In loops or
headless pipelines you will accumulate figures (matplotlib warning after ~20).
Use:
  * `close=True` on supported plotters (passes through to save + close)
  * or `scat.pl.save_all_figures(..., close=True)`
  * or explicitly `plt.close(fig)` after saving/showing
See individual functions for the `close=` parameter.

Internal `set_style()` is called by plotting functions **only when use_style=True**.
Default is now use_style=False to avoid surprising global rcParams changes in notebooks
and shared environments.

- To get the package publication style automatically, either:
  * pass use_style=True to a plotting call, or
  * call `scat.pl.set_style()` once yourself early in the script/notebook, or
  * use `with scat.pl.style_context(...): ...` (recommended for scoped, no side effects).
- All plotting functions still respect use_style= , ax=, figsize=, save_path=, show= for
  consistency.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable, Mapping
from contextlib import contextmanager, nullcontext, suppress
from typing import Any

import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize
from scipy import sparse

from ._utils import (
    LEGACY_VELOCITY_DELTA_COL,
    LEGACY_VELOCITY_RESIDUAL_COL,
    UNSPLICED_EXCESS_DELTA_COL,
    UNSPLICED_EXCESS_RESIDUAL_COL,
    _resolve_results_column,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def set_style(
    fontfamily="sans-serif",
    fonts=None,
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
    script or notebook if you want package style. Plotting functions no longer call
    it by default (use_style=False). Pass use_style=True or use style_context()
    to opt in per call or per block.

    If you want to limit the scope of style changes, either:
      * call `scat.pl.set_style()` once yourself and then pass use_style=False to individual
        plotting calls, or
      * use `with scat.pl.style_context(...): ...` around blocks of code (including calls
        to scat.pl.* functions inside the block).

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
    if fonts is None:
        fonts = ["Arial", "Helvetica", "DejaVu Sans"]
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


def _style_context_if(use_style: bool):
    """Safe context (nullcontext when not using style)."""
    return style_context() if use_style else nullcontext()


def _maybe_repel_labels(
    texts: list,
    x: np.ndarray,
    y: np.ndarray,
    ax,
    *,
    label_repel: bool = True,
) -> None:
    """Apply adjustText when available; optional fallback with overlap warning."""
    if not texts:
        return
    if not label_repel:
        return
    try:
        from adjustText import adjust_text
    except ImportError:
        logger.warning(
            "adjustText is not installed; gene labels may overlap. "
            "Install with: pip install adjustText — or pass label_repel=False."
        )
        return
    adjust_text(
        texts,
        x=x,
        y=y,
        arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.8, "alpha": 0.8},
        ax=ax,
    )


@contextmanager
def figure_export_context(
    directory: str,
    *,
    dpi: int = 300,
    fmt: str = "pdf",
    bbox_inches: str = "tight",
):
    """
    Context manager for batch-saving figures in multi-panel workflows.

    Example::

        with scat.pl.figure_export_context("figures/out") as export:
            fig, ax = scat.pl.comet_plot(all_results, show=False)
            export.save(fig, "comet")
            fig2, ax2 = scat.pl.volcano_plot(all_results, show=False)
            export.save(fig2, "volcano")
    """
    from pathlib import Path

    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: list[str] = []

    class _Exporter:
        def save(self, fig, name: str) -> str:
            path = out_dir / f"{name}.{fmt}"
            fig.savefig(path, dpi=dpi, bbox_inches=bbox_inches)
            saved.append(str(path))
            return str(path)

    exporter = _Exporter()
    try:
        yield exporter
    finally:
        if saved:
            logger.info("figure_export_context: saved %d file(s) to %s", len(saved), out_dir)


def save_all_figures(
    figures: dict,
    directory: str,
    *,
    dpi: int = 300,
    fmt: str = "pdf",
    bbox_inches: str = "tight",
    close: bool = True,
) -> list[str]:
    """
    Batch-save a mapping of ``name -> matplotlib Figure``.

    Convenience wrapper around :func:`figure_export_context` for notebooks that
    already hold figure objects.

    Example::

        paths = scat.pl.save_all_figures(
            {"comet": fig1, "volcano": fig2},
            "figures/out",
        )
    """
    import matplotlib.pyplot as plt

    saved: list[str] = []
    with figure_export_context(directory, dpi=dpi, fmt=fmt, bbox_inches=bbox_inches) as export:
        for name, fig in figures.items():
            if fig is None:
                logger.warning("save_all_figures: skipping None figure for %r", name)
                continue
            saved.append(export.save(fig, name))
            if close:
                plt.close(fig)
    return saved


# -----------------------------------------------------------------------------
# Internal helpers for robust plotting (column validation, safe math, parsing)
# -----------------------------------------------------------------------------


def _excess_delta_col(df: pd.DataFrame) -> str:
    return _resolve_results_column(df, UNSPLICED_EXCESS_DELTA_COL, LEGACY_VELOCITY_DELTA_COL)


def _excess_residual_col(df: pd.DataFrame) -> str:
    return _resolve_results_column(df, UNSPLICED_EXCESS_RESIDUAL_COL, LEGACY_VELOCITY_RESIDUAL_COL)


def _require_columns(df, columns, func_name="plot"):
    """Raise clear error if required columns are missing from df."""
    if df is None:
        raise ValueError(f"{func_name} requires a DataFrame but got None")
    missing = [c for c in columns if c not in getattr(df, "columns", [])]
    if missing:
        avail = list(df.columns) if hasattr(df, "columns") else []
        raise ValueError(
            f"{func_name} requires columns {columns}, missing: {missing}. "
            f"Available columns: {avail}"
        )


def _guard_plot_dataframe(
    df,
    required_cols: list[str],
    func_name: str,
    *,
    message: str | None = None,
    ax=None,
    figsize=(6, 4),
    dpi: int = 150,
) -> tuple[Any, Any, bool]:
    """Return placeholder (fig, ax, guarded=True) when input is empty or missing columns."""
    if df is None or (hasattr(df, "empty") and df.empty):
        logger.warning("%s: empty DataFrame.", func_name)
        if ax is None:
            fig, ax = _empty_placeholder_fig(message or "No data to plot", figsize=figsize, dpi=dpi)
            return fig, ax, True
        return ax.figure, ax, True
    missing = [c for c in required_cols if c not in getattr(df, "columns", [])]
    if missing:
        logger.warning("%s: missing columns %s", func_name, missing)
        if ax is None:
            fig, ax = _empty_placeholder_fig(
                message or f"Missing columns: {', '.join(missing)}",
                figsize=figsize,
                dpi=dpi,
            )
            return fig, ax, True
        return ax.figure, ax, True
    return None, None, False


def _venn_exclusive_region_counts(sets_list: list[set]) -> dict[tuple[int, ...], int]:
    """Element counts for each exclusive Venn region (one per non-empty subset of sets)."""
    n = len(sets_list)
    regions: dict[tuple[int, ...], int] = {}
    for mask in range(1, 2**n):
        inc_idx = [i for i in range(n) if mask & (1 << i)]
        exc_idx = [i for i in range(n) if not (mask & (1 << i))]
        region = set.intersection(*(sets_list[i] for i in inc_idx))
        for j in exc_idx:
            region -= sets_list[j]
        regions[tuple(inc_idx)] = len(region)
    return regions


def _venn_region_label_coords(n: int) -> dict[tuple[int, ...], tuple[float, float]]:
    """Approximate (x, y) label anchors for the built-in circle Venn layouts (n=1..4)."""
    if n == 1:
        return {(0,): (0.0, 0.0)}
    if n == 2:
        return {
            (0,): (-0.75, 0.15),
            (1,): (0.75, 0.15),
            (0, 1): (0.0, 0.15),
        }
    if n == 3:
        return {
            (0,): (-0.85, 0.35),
            (1,): (0.85, 0.35),
            (2,): (0.0, -0.95),
            (0, 1): (0.0, 0.45),
            (0, 2): (-0.35, -0.05),
            (1, 2): (0.35, -0.05),
            (0, 1, 2): (0.0, 0.05),
        }
    if n == 4:
        return {
            (0,): (-0.95, 0.55),
            (1,): (0.95, 0.55),
            (2,): (-0.95, -0.55),
            (3,): (0.95, -0.55),
            (0, 1): (0.0, 0.58),
            (2, 3): (0.0, -0.58),
            (0, 2): (-0.58, 0.05),
            (1, 3): (0.58, 0.05),
            (0, 3): (-0.05, 0.05),
            (1, 2): (0.05, 0.05),
            (0, 1, 2): (-0.22, 0.22),
            (0, 1, 3): (0.22, 0.22),
            (0, 2, 3): (-0.22, -0.22),
            (1, 2, 3): (0.22, -0.22),
            (0, 1, 2, 3): (0.0, 0.05),
        }
    return {}


def _safe_neg_log10(x, minval=1e-300):
    """Safe -log10 with clipping for zero/near-zero p-values and non-numeric safety.

    Accepts scalar, list/Series/ndarray etc. NaN inputs -> NaN outputs (callers dropna if wanted).
    For scalar/0-d input returns a python float; otherwise returns ndarray (shape preserved).
    """
    # Coerce everything to 1-d float ndarray (handles scalar, list, Series, 0-d, matrix etc)
    arr = np.asarray(pd.to_numeric(x, errors="coerce"), dtype=float).ravel()
    clipped = np.clip(arr, a_min=minval, a_max=None)
    res = -np.log10(clipped)

    # Decide return style based on original input
    was_scalar = (
        np.isscalar(x)
        or (isinstance(x, (int, float, np.number)) and not hasattr(x, "__len__"))
        or (hasattr(x, "ndim") and getattr(x, "ndim", 1) == 0)
        or (hasattr(x, "shape") and getattr(x, "shape", (1,)) == ())
    )
    if was_scalar:
        return float(res[0]) if res.size > 0 else np.nan
    # For sequence-like input, return 1-d ndarray (pandas assignment accepts it)
    return res


def _parse_gene_ratio(x):
    """Convert GeneRatio '3/120' strings (clusterProfiler etc.) or numeric to float."""
    if pd.isna(x):
        return np.nan
    if isinstance(x, str) and "/" in x:
        try:
            a, b = x.split("/", 1)
            return float(a) / float(b)
        except Exception:
            return np.nan
    return pd.to_numeric(x, errors="coerce")


def _finite_scatter_values(values, fill: float = 0.0) -> np.ndarray:
    """Coerce scatter size/color arrays to finite floats (avoids matplotlib size warnings)."""
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    bad = ~np.isfinite(arr)
    if bad.any():
        arr = arr.copy()
        arr[bad] = fill
    return arr


def _save_and_maybe_show(
    fig, save_path=None, dpi=300, show=True, created=True, transparent=True, close=False
):
    """Internal: centralized save + conditional show + optional close.

    close=True is intended for batch/non-interactive pipelines to avoid matplotlib
    "too many figures open" warnings and memory growth. Default remains False for
    notebook/interactive compatibility (user keeps the returned fig alive).
    """
    if save_path:
        from pathlib import Path

        out = Path(str(save_path))
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logger.info("Figure saved to %s", out)
    if created and show:
        plt.show()
    if created and close:
        plt.close(fig)


def _empty_placeholder_fig(message="No data to plot", figsize=(6, 4), dpi=150):
    """Create a minimal placeholder figure so callers always receive (fig, ax) even on empty data.
    This improves UX and return-type stability vs returning (None, None).
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, fontsize=10)
    ax.axis("off")
    sns.despine(ax=ax, top=True, right=True, left=True, bottom=True)
    return fig, ax


def comet_plot(
    df,
    top_n=12,
    save_path=None,
    title="Active Transcription Drivers",
    point_scale=1.0,
    min_size=2,
    max_size=180,
    s: float
    | None = None,  # fixed point size (overrides variable sizing by active_score); common control in omicverse-style APIs
    alpha: float = 0.85,  # point transparency (omicverse often uses ~0.5 for clean dense plots)
    figsize=(8, 6),
    dpi=300,
    fontsize=12,
    cmap="coolwarm",
    ax=None,
    show: bool = True,
    use_style: bool = False,
    positive_logfc_only: bool = True,
    return_data: bool = False,
    label_repel: bool = True,
    label_fontsize: float | None = None,
    min_label_score: float | None = None,
    close: bool = False,
):
    """
    Comet plot of log fold change vs. bias-corrected unspliced residual.

    Point size and color are mapped to the active score. Designed to produce
    clear figures suitable for scientific publications with minimal further
    editing.

    Size control (referencing common patterns in omicverse.pl.* for direct control):
      - `s`: if provided, use a **fixed** point size for all points (in points^2).
        This is the simplest way to make everything small (e.g. s=3 or s=1).
      - `point_scale`: overall multiplier for the variable size calculation.
      - `min_size` / `max_size`: hard bounds. Use min_size=1 to allow the tiniest
        background points when using variable sizing by active_score.

    `positive_logfc_only=True` (default) keeps only logFC > 0 (classic "active drivers"
    comet view). Set False to see the full logFC vs residual scatter including
    negative logFC genes.

    Returns
    -------
    (fig, ax) : always a matplotlib figure and axes. If no valid genes remain
    after coercion / filtering, returns a placeholder figure with a message
    (instead of (None, None)) for better caller ergonomics.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        If provided, plot into this axes instead of creating a new figure.
        Useful for embedding in multi-panel publication figures.
    """
    logger.info("Generating comet plot...")
    _style_ctx = None
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()

    if df is not None and hasattr(df, "columns"):
        residual_col = _resolve_results_column(
            df,
            UNSPLICED_EXCESS_RESIDUAL_COL,
            LEGACY_VELOCITY_RESIDUAL_COL,
            required=False,
        )
    else:
        residual_col = UNSPLICED_EXCESS_RESIDUAL_COL
    guard_cols = ["logFC", "active_score", residual_col]
    fig_g, ax_g, guarded = _guard_plot_dataframe(
        df, guard_cols, "comet_plot", ax=ax, figsize=figsize, dpi=dpi
    )
    if guarded:
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig_g, save_path=save_path, dpi=dpi, show=show, created=ax is None, close=close
        )
        return fig_g, ax_g

    residual_col = _excess_residual_col(df)
    _require_columns(df, ["logFC", residual_col, "active_score"], "comet_plot")

    if top_n <= 0:
        raise ValueError("top_n must be positive.")
    if min_size < 0 or max_size <= min_size:
        raise ValueError("Require 0 <= min_size <= max_size.")
    if point_scale <= 0:
        raise ValueError("point_scale must be positive.")
    if s is not None and s <= 0:
        raise ValueError("s must be positive.")

    plot_df = df.copy()
    for c in ["logFC", residual_col, "active_score"]:
        if c in plot_df.columns:
            plot_df[c] = pd.to_numeric(plot_df[c], errors="coerce")
    plot_df = plot_df.dropna(subset=["logFC", residual_col, "active_score"])

    if positive_logfc_only:
        plot_df = plot_df[plot_df["logFC"] > 0].copy()

    if plot_df.empty:
        logger.warning(
            "No genes to plot after filtering (positive_logfc_only=%s or missing data).",
            positive_logfc_only,
        )
        if ax is None:
            fig, ax = _empty_placeholder_fig("No genes to plot after filtering")
            _created_fig = True
        else:
            fig = ax.figure
            _created_fig = False
            # still return the provided ax (user can decide)
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
        )
        return fig, ax

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    # Size scaling (inspired by flexible controls in omicverse.pl for volcano/comet-style plots).
    # - If user passes `s=...`, use fixed size (very common request: "just make all points small").
    # - Otherwise use active_score-powered variable sizing with user min/max.
    # Defensively clip active_score to >=0 before power (negative scores would produce NaN).
    if s is not None:
        sizes = np.full(len(plot_df), float(s) * point_scale)
        effective_min = min(1.0, min_size)
        sizes = np.clip(sizes, effective_min, max_size)
    else:
        score_for_size = np.clip(
            pd.to_numeric(plot_df["active_score"], errors="coerce").fillna(0), 0, None
        )
        raw_sizes = score_for_size**1.6 * 35 * point_scale + 3 * point_scale
        sizes = np.clip(raw_sizes, min_size, max_size)
    sizes = _finite_scatter_values(sizes, fill=float(min_size))

    # Light omicverse-style diagnostics (non-intrusive)
    if len(plot_df) > 500 and (s is None) and point_scale > 0.3 and min_size > 3:
        logger.info(
            "Many points detected (%d). Consider s=2 or point_scale=0.1 + min_size=1 "
            "for cleaner comet plot (inspired by omicverse.pl best practices).",
            len(plot_df),
        )

    color_vals = _finite_scatter_values(plot_df["active_score"], fill=0.0)
    scatter = ax.scatter(
        x=plot_df["logFC"],
        y=plot_df[residual_col],
        c=color_vals,
        s=sizes,
        cmap=cmap,
        alpha=alpha,
        edgecolors="#444444",
        linewidth=0.5,
        zorder=3,
    )

    ax.axhline(0, color="#999999", linestyle="--", linewidth=1, alpha=0.5, zorder=1)
    ax.axvline(0, color="#999999", linestyle="--", linewidth=1, alpha=0.5, zorder=1)

    label_pool = plot_df
    if min_label_score is not None and "active_score" in label_pool.columns:
        label_pool = label_pool[label_pool["active_score"] >= float(min_label_score)]
    top_genes = label_pool.nlargest(top_n, "active_score")
    lbl_fs = float(label_fontsize) if label_fontsize is not None else max(8, fontsize - 2)
    texts = []
    for idx, row in top_genes.iterrows():
        txt = ax.text(
            row["logFC"],
            row[residual_col],
            f"{idx}",
            fontsize=lbl_fs,
            fontweight="bold",
            color="#111111",
            bbox={"boxstyle": "square,pad=0.1", "fc": "none", "ec": "none"},
        )
        texts.append(txt)

    _maybe_repel_labels(
        texts,
        plot_df["logFC"].values,
        plot_df[residual_col].values,
        ax,
        label_repel=label_repel,
    )

    ax.set_xlabel("Log2 Fold Change", fontsize=fontsize, fontweight="bold")
    ax.set_ylabel("Bias-corrected Unspliced Residual", fontsize=fontsize, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=15)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.6, pad=0.03, aspect=20)
    cbar.set_label(
        "Active Score", fontsize=max(9, fontsize - 1), fontweight="bold", rotation=270, labelpad=15
    )
    cbar.outline.set_visible(False)

    sns.despine(ax=ax, top=True, right=True)
    ax.spines["left"].set_position(("outward", 6))
    ax.spines["bottom"].set_position(("outward", 6))

    # constrained_layout at creation + bbox_inches on save handles colorbar cleanly.
    # (avoid tight_layout after colorbar)

    _save_and_maybe_show(
        fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
    )
    if return_data:
        return fig, ax, plot_df
    return fig, ax


def volcano_3d(
    df,
    top_n=8,
    save_path=None,
    point_scale=1.0,
    min_size=2,
    max_size=160,
    s: float | None = None,  # fixed point size (direct control, omicverse reference)
    alpha: float = 0.8,
    title="3D Active Volcano Plot",
    figsize=(10, 8),
    dpi=300,
    fontsize=11,
    cmap="coolwarm",
    ax=None,
    show: bool = True,
    use_style: bool = False,
    return_data: bool = False,
    close: bool = False,
):
    """
    3D volcano-style view (logFC, -log10(p_adj), velocity residual).

    If `ax` (a 3D axes) is provided, plot into it.

    Size control (omicverse-style):
      - `s`: fixed point size for all points.
      - `point_scale`, `min_size`, `max_size` for variable sizing by active_score.
    Use s=2 or min_size=1 + small point_scale for tiny background points.
    """
    logger.info("Generating 3D volcano plot...")
    _style_ctx = None
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()

    if df is None or (hasattr(df, "empty") and df.empty):
        logger.warning("Empty results DataFrame passed to volcano_3d.")
        if ax is None:
            fig, ax = _empty_placeholder_fig("No valid genes to plot", figsize=(8, 6))
            _created_fig = True
        else:
            fig = ax.figure
            _created_fig = False
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
        )
        return fig, ax

    residual_col = _excess_residual_col(df)
    _require_columns(df, ["logFC", "p_adj", residual_col, "active_score"], "volcano_3d")

    plot_df = df.copy()
    for c in ["logFC", "p_adj", residual_col, "active_score"]:
        if c in plot_df.columns:
            plot_df[c] = pd.to_numeric(plot_df[c], errors="coerce")

    # Filter invalid p_adj (<0 or >1)
    if "p_adj" in plot_df.columns:
        invalid_p = (plot_df["p_adj"] < 0) | (plot_df["p_adj"] > 1)
        if invalid_p.any():
            logger.warning("Dropping %d rows with p_adj outside [0, 1].", int(invalid_p.sum()))
            plot_df = plot_df.loc[~invalid_p].copy()

    plot_df = plot_df.dropna(subset=["logFC", "p_adj", residual_col, "active_score"])
    plot_df["neg_log_pval"] = _safe_neg_log10(plot_df["p_adj"])
    plot_df = plot_df.dropna(subset=["neg_log_pval"])

    if plot_df.empty:
        logger.warning("No valid genes to plot after numeric coercion/dropna in volcano_3d.")
        if ax is None:
            fig, ax = _empty_placeholder_fig("No valid genes to plot", figsize=(8, 6))
            _created_fig = True
        else:
            fig = ax.figure
            _created_fig = False
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
        )
        return fig, ax

    if ax is None:
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax = fig.add_subplot(111, projection="3d")
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    if s is not None:
        sizes = np.full(len(plot_df), float(s) * point_scale)
        effective_min = min(1.0, min_size)
        sizes = np.clip(sizes, effective_min, max_size)
    else:
        # clip active_score to >=0 before power to avoid NaN from negative ** exponent
        score_for_size = np.clip(
            pd.to_numeric(plot_df["active_score"], errors="coerce").fillna(0), 0, None
        )
        raw_sizes = score_for_size**1.4 * 18 * point_scale + 3 * point_scale
        sizes = np.clip(raw_sizes, min_size, max_size)

    if len(plot_df) > 500 and (s is None) and point_scale > 0.3:
        logger.info(
            "3D volcano: %d points. For performance and clarity try s=2 or small point_scale + min_size=1.",
            len(plot_df),
        )

    scatter = ax.scatter(
        plot_df["logFC"],
        plot_df["neg_log_pval"],
        plot_df[residual_col],
        c=plot_df["active_score"],
        s=sizes,
        cmap=cmap,
        alpha=alpha,
        edgecolors="#444444",
        linewidth=0.4,
        zorder=3,
    )

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        axis.line.set_color((1.0, 1.0, 1.0, 0.0))

    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        with suppress(
            Exception
        ):  # private API; future matplotlib may change this, ignore gracefully
            axis._axinfo["grid"].update({"color": "#E5E5E5", "linestyle": "-"})

    top_genes = plot_df.nlargest(top_n, "active_score")
    # Use data range for offsets (robust to small/negative ranges)
    if len(plot_df) > 1:
        x_rng = plot_df["logFC"].max() - plot_df["logFC"].min()
        z_rng = plot_df[residual_col].max() - plot_df[residual_col].min()
        x_offset = (x_rng * 0.03) if x_rng > 0 else 0.1
        z_offset = (z_rng * 0.03) if z_rng > 0 else 0.15
    else:
        x_offset = 0.1
        z_offset = 0.15

    for idx, row in top_genes.iterrows():
        px, py, pz = row["logFC"], row["neg_log_pval"], row[residual_col]
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

    # 3D subplots have limited layout engine support; keep tight only for new fig
    # (colorbar is on the figure).

    if _style_ctx is not None:
        _style_ctx.__exit__(None, None, None)
    _save_and_maybe_show(
        fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
    )
    if return_data:
        return fig, ax, plot_df
    return fig, ax


def enrich_dotplot(
    enrich_df,
    top_n=15,
    show_terms: int | str | list[str] | tuple[str, ...] | None = None,
    title="Enrichment Dotplot",
    save_path=None,
    figsize=(7, 8),
    dpi=300,
    fontsize=12,
    x="GeneRatio",
    color_by="p.adjust",
    size_by="Count",
    cmap="viridis_r",
    dot_max: float | None = None,
    dot_min: float | None = None,
    smallest_dot: float = 0.0,
    ax=None,
    show: bool = True,
    use_style: bool = False,
    cluster_col: str | None = None,
    facet_by_cluster: bool = False,
    return_data: bool = False,
    close: bool = False,
):
    """
    Dotplot for enrichment results (clusterProfiler style).

    Common display choices (all columns from run_enrichment / run_kegg / run_gsea are available):
      - `x`: what to plot on the x-axis. Supported / nice values:
          "GeneRatio" (default for ORA), "FoldEnrichment", "Count", "-log10(p.adj)", "NES" (auto default for GSEA).
          You can also pass any other numeric column present in the dataframe.
          Example: `x="NES"` for GSEA results.
      - `size_by`: controls dot size (default "Count"). Common: "Count", "GeneRatio".
      - `color_by`: controls dot color (default "p.adjust" for ORA;
        for GSEA results with "NES" it will default to "NES" with a diverging colormap).
        Smaller p-values are usually more interesting.
      - `dot_max`, `dot_min`, `smallest_dot`: omicverse-style controls for dot size range
        (see omicverse.pl.dotplot for the excellent reference implementation).

    Legend handling (colorbar for p-value + size legend for Count/GeneRatio) uses
    constrained_layout + careful bbox_to_anchor upper-right placement for the size legend.
    This follows the patterns from gseapy (zqfang/gseapy plot.DotPlot) and omicverse.pl.dotplot
    to avoid the two legend elements overlapping on the right side of the figure.

    `show_terms` gives clusterProfiler-like flexibility:
      - int: show top N terms (overrides top_n)
      - "auto": intelligently select terms with p.adjust < 0.05 and Count >= 2 (then top_n of those,
        sorted by significance + size). Falls back gracefully.
      - list/tuple of str: show exactly the matching terms (match on Term or Description;
        order of the list is respected when possible). This is analogous to
        `dotplot(..., showCategory = c("term1", "term2"))`.

    Multi-group / compare support (new):
      - If the input df contains a "Cluster" column (or you pass `cluster_col="Cluster"`),
        terms are automatically prefixed with "[Cluster] " so groups are visually distinct.
      - `show_terms="auto"` (and int) now collects top terms intelligently across clusters
        when a Cluster column is present (union of interesting terms per cluster).
      - `facet_by_cluster=True` will produce a grid of subplots (one per cluster) using
        the same rich dotplot logic. Excellent for publication multi-panel figures.
      - The returned df from `scat.compare_enrichment(...)` or `concat_compare_results(...)`
        works directly.

    `top_n` is still supported for the common "top N" case (when show_terms is None).
    Supports `ax` for embedding in publication multi-panel figures.
    """
    _style_ctx = None
    fig_g, ax_g, guarded = _guard_plot_dataframe(
        enrich_df,
        [],
        "enrich_dotplot",
        message="No enrichment terms to plot",
        ax=ax,
        figsize=figsize,
        dpi=dpi,
    )
    if guarded:
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig_g, save_path=save_path, dpi=dpi, show=show, created=ax is None, close=close
        )
        return fig_g, ax_g

    if "Term" not in enrich_df.columns and "Description" not in enrich_df.columns:
        logger.warning("enrich_dotplot: missing 'Term' and 'Description' columns.")
        if ax is None:
            fig, ax = _empty_placeholder_fig("Missing Term/Description columns")
            _created_fig = True
        else:
            fig = ax.figure
            _created_fig = False
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
        )
        return fig, ax

    # --- Multi-cluster detection (compareCluster style) ---
    if cluster_col is None:
        for cand in ("Cluster", "cluster", "Group", "group"):
            if cand in enrich_df.columns:
                cluster_col = cand
                break

    has_cluster = bool(
        cluster_col and cluster_col in enrich_df.columns and enrich_df[cluster_col].nunique() > 1
    )

    # Facet support: delegate to subplots (reuses full logic, avoids code duplication)
    if facet_by_cluster and has_cluster:
        clusters = [str(x) for x in pd.Series(enrich_df[cluster_col]).dropna().unique().tolist()]
        n = len(clusters)
        if n > 1:
            import math as _math

            ncols = min(3, n)
            nrows = _math.ceil(n / ncols)
            fig_w = figsize[0] * min(ncols, 2.5)
            fig_h = figsize[1] * max(0.7, nrows * 0.85)
            fig, axes = plt.subplots(
                nrows, ncols, figsize=(fig_w, fig_h), constrained_layout=True, dpi=dpi
            )
            axes = np.asarray(axes).ravel().tolist() if hasattr(axes, "ravel") else [axes]
            _created_fig = True
            for i, cl in enumerate(clusters):
                ax_i = axes[i] if i < len(axes) else axes[-1]
                sub = enrich_df[enrich_df[cluster_col].astype(str) == cl].copy()
                # recurse without facet to fill the ax
                enrich_dotplot(
                    sub,
                    top_n=top_n,
                    show_terms=show_terms,
                    title=f"{title} — {cl}" if title else str(cl),
                    save_path=None,
                    figsize=(fig_w / ncols * 0.95, fig_h / nrows * 0.9),
                    dpi=dpi,
                    fontsize=fontsize,
                    x=x,
                    color_by=color_by,
                    size_by=size_by,
                    cmap=cmap,
                    dot_max=dot_max,
                    dot_min=dot_min,
                    smallest_dot=smallest_dot,
                    ax=ax_i,
                    show=False,
                    use_style=use_style,
                    cluster_col=cluster_col,
                    facet_by_cluster=False,
                )
            # hide unused axes
            for j in range(n, len(axes)):
                axes[j].axis("off")
            if _style_ctx is not None:
                _style_ctx.__exit__(None, None, None)
            _save_and_maybe_show(
                fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
            )
            return fig, axes[0] if axes else fig

    if dot_min is not None and dot_max is not None and dot_min > dot_max:
        raise ValueError("dot_min must be <= dot_max.")
    if smallest_dot < 0 or smallest_dot > 1:
        raise ValueError("smallest_dot must be between 0 and 1 (inclusive).")
    if dot_min is not None and dot_min < 0:
        logger.warning(
            "dot_min < 0 may produce unexpected dot sizes for non-negative metrics like Count."
        )
    if dot_max is not None and dot_max < 0:
        logger.warning("dot_max < 0 may produce unexpected dot sizes.")

    logger.info("Generating enrichment dotplot...")
    _style_ctx = style_context() if use_style else nullcontext()
    _style_ctx.__enter__() if use_style else None

    # Auto defaults for GSEA results (NES present)
    is_gsea = ("NES" in enrich_df.columns) or ("ES" in enrich_df.columns)
    if x == "GeneRatio" and is_gsea and "NES" in enrich_df.columns:
        x = "NES"  # sensible default for GSEA
    if color_by in ("Adjusted P-value", "p.adjust") and is_gsea and "NES" in enrich_df.columns:
        color_by = "NES"

    # clusterProfiler-style term selection (show_terms takes precedence)
    if show_terms is not None:
        if isinstance(show_terms, int):
            if has_cluster:
                # Collect top N per cluster, then take union (preserves diversity across groups)
                parts = []
                for _cl, sub in enrich_df.groupby(cluster_col, sort=False):
                    parts.append(sub.head(show_terms))
                plot_df = (
                    pd.concat(parts).drop_duplicates()
                    # Do not apply another global .head() here: the per-cluster
                    # heads already ensure representation. A cluster-order-dependent
                    # head would silently drop later clusters' terms.
                )
            else:
                plot_df = enrich_df.head(show_terms).copy()
        elif isinstance(show_terms, str) and show_terms.lower() == "auto":
            # Auto: prefer statistically interesting terms (p.adjust + count)
            # When multi-cluster, collect interesting terms per cluster then union
            padj_col = None
            for cand in ("p.adjust", "Adjusted P-value", "p_adj", "padj"):
                if cand in enrich_df.columns:
                    padj_col = cand
                    break
            if padj_col is None:
                padj_col = "p.adjust"
            count_col = "Count" if "Count" in enrich_df.columns else None

            if has_cluster:
                selected_idx = []
                for _cl, sub in enrich_df.groupby(cluster_col, sort=False):
                    if padj_col in sub.columns:
                        padj = pd.to_numeric(sub[padj_col], errors="coerce").fillna(1.0)
                    else:
                        padj = pd.Series(1.0, index=sub.index)
                    if count_col and count_col in sub.columns:
                        cnt = pd.to_numeric(sub[count_col], errors="coerce").fillna(0)
                    else:
                        cnt = pd.Series(0, index=sub.index)
                    auto_mask = (padj < 0.05) & (cnt >= 2)
                    if auto_mask.any():
                        sub_cand = sub[auto_mask].copy()
                        if padj_col in sub_cand.columns:
                            sub_cand = sub_cand.sort_values(
                                by=[padj_col, count_col or sub_cand.columns[0]],
                                ascending=[True, False],
                            )
                        selected_idx.extend(sub_cand.head(top_n).index.tolist())
                    else:
                        selected_idx.extend(sub.head(top_n).index.tolist())
                plot_df = enrich_df.loc[list(dict.fromkeys(selected_idx))].copy()
                # per-cluster selection already performed; avoid order-dependent global head
                # that can drop entire later clusters from the union.
            else:
                if padj_col in enrich_df.columns:
                    padj = pd.to_numeric(enrich_df[padj_col], errors="coerce").fillna(1.0)
                else:
                    padj = pd.Series(1.0, index=enrich_df.index)
                if count_col and count_col in enrich_df.columns:
                    cnt = pd.to_numeric(enrich_df[count_col], errors="coerce").fillna(0)
                else:
                    cnt = pd.Series(0, index=enrich_df.index)
                auto_mask = (padj < 0.05) & (cnt >= 2)
                if auto_mask.sum() == 0:
                    if padj_col in enrich_df.columns:
                        tmp = enrich_df.copy()
                        tmp["_key"] = -np.log10(
                            pd.to_numeric(tmp[padj_col], errors="coerce").fillna(1) + 1e-300
                        )
                        plot_df = (
                            tmp.sort_values("_key", ascending=False)
                            .head(top_n)
                            .drop(columns=["_key"], errors="ignore")
                            .copy()
                        )
                    else:
                        plot_df = enrich_df.head(top_n).copy()
                else:
                    cand = enrich_df[auto_mask].copy()
                    if padj_col in cand.columns:
                        cand = cand.sort_values(
                            by=[padj_col, "Count" if "Count" in cand.columns else cand.columns[0]],
                            ascending=[True, False],
                        )
                    plot_df = cand.head(top_n).copy()
        else:
            # list/tuple of explicit terms
            wanted = {
                str(x).strip().lower()
                for x in (show_terms if not isinstance(show_terms, str) else [show_terms])
            }

            def _matches(row):
                t = str(row.get("Term", "")).strip().lower()
                d = str(row.get("Description", "")).strip().lower()
                return any(w in t or w in d for w in wanted)

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
                plot_df = plot_df.sort_values("_sel_order").drop(
                    columns=["_sel_order"], errors="ignore"
                )
    else:
        if has_cluster:
            # default: take a reasonable number across clusters
            parts = []
            per_cl = max(3, top_n // max(1, enrich_df[cluster_col].nunique()))
            for _cl, sub in enrich_df.groupby(cluster_col, sort=False):
                parts.append(sub.head(per_cl))
            plot_df = pd.concat(parts).drop_duplicates()
            # Removed the final .head(top_n*2): it was applied after concat in
            # cluster appearance order and could drop all terms from later clusters.
            # Per-cluster caps + dedup is sufficient for representation.
        else:
            plot_df = enrich_df.head(top_n).copy()

    plot_df = plot_df.iloc[::-1]  # visual: top at top of y axis for horizontal dotplot

    # Parse GeneRatio early (many sources emit "3/120" strings; max/min would fail otherwise)
    if "GeneRatio" in plot_df.columns:
        plot_df["GeneRatio"] = plot_df["GeneRatio"].apply(_parse_gene_ratio)

    def clean_term(text):
        text = str(text).split(" (GO:")[0].split(" (KEGG")[0]
        return text[:50] + "..." if len(text) > 50 else text

    # Robust term column selection (clusterProfiler/gseapy use Term or Description; some use ID)
    term_col = None
    for c in ["Term", "Description", "term", "description", "ID", "id"]:
        if c in plot_df.columns:
            term_col = c
            break
    if term_col is None:
        raise ValueError(
            "enrich_dotplot requires a term-like column. "
            "Expected one of: ['Term', 'Description', 'term', 'description', 'ID']. "
            f"Available: {list(plot_df.columns)}"
        )

    plot_df = plot_df.copy()  # ensure writable
    plot_df["Term_Clean"] = plot_df[term_col].astype(str).apply(clean_term)

    # Prefix with cluster for visual grouping (compareCluster style) when not faceting
    if has_cluster and cluster_col in plot_df.columns and not facet_by_cluster:
        cl_str = plot_df[cluster_col].astype(str)
        # avoid double prefix if user already has it
        if not plot_df["Term_Clean"].str.contains(r"^\[.*\]").any():
            plot_df["Term_Clean"] = "[" + cl_str + "] " + plot_df["Term_Clean"]

    pval_candidates = ["p.adjust", "Adjusted P-value", "p_adj", "padj", "FDR_qval", "pvalue"]
    pval_col = next((c for c in pval_candidates if c in plot_df.columns), None)
    if pval_col is None:
        pval_col = plot_df.columns[0]

    # Filter clearly invalid p-values for the chosen p column (enrichment results should have p in [0,1])
    if pval_col and pval_col in plot_df.columns:
        plot_df[pval_col] = pd.to_numeric(plot_df[pval_col], errors="coerce")
        invalid_p = (plot_df[pval_col] < 0) | (plot_df[pval_col] > 1)
        if invalid_p.any():
            logger.warning(
                "Dropping %d rows with %s outside [0, 1].", int(invalid_p.sum()), pval_col
            )
            plot_df = plot_df.loc[~invalid_p].copy()

    count_candidates = ["Count", "Size", "leadingEdge_count"]
    size_col = next((c for c in count_candidates if c in plot_df.columns), None)

    requested_x = x
    if x == "-log10(p.adj)" or x == "-log10(p.adjust)":
        if pval_col:
            plot_df["neg_log_padj"] = _safe_neg_log10(plot_df[pval_col])
            x_col = "neg_log_padj"
            x_label = f"-log10({pval_col})"
        else:
            x_col = x if x in plot_df.columns else "GeneRatio"
            x_label = x_col
    else:
        x_col = x if x in plot_df.columns else "GeneRatio"
        # Nice labels for the common x choices users care about (GeneRatio / FoldEnrichment / Count)
        if x_col == "GeneRatio":
            x_label = "Gene Ratio"
        elif x_col == "FoldEnrichment":
            x_label = "Fold Enrichment"
        elif x_col == "Count":
            x_label = "Count"
        else:
            x_label = x_col

    # Heuristic: if the *effective* x is GeneRatio ... (now safe after parse)
    if (
        x_col == "GeneRatio"
        and "FoldEnrichment" in plot_df.columns
        and requested_x in (None, "GeneRatio", "generatio", "gene_ratio")
    ):
        gene_ratio_range = plot_df["GeneRatio"].max() - plot_df["GeneRatio"].min()
        if pd.notna(gene_ratio_range) and gene_ratio_range < 0.08:
            logger.warning(
                "⚠️ GeneRatio values have very low variation. Switching to 'FoldEnrichment'."
            )
            x_col = "FoldEnrichment"
            x_label = "Fold Enrichment"

    if size_by in plot_df.columns:
        size_col = size_by
    if size_col is None or size_col not in plot_df.columns:
        logger.warning(
            "No valid size column found (Count/Size/leadingEdge_count or size_by). Using constant dot size."
        )
        plot_df["_dot_size"] = 1.0
        size_col = "_dot_size"

    if color_by in plot_df.columns:
        color_col = color_by
    else:
        color_col = pval_col if pval_col else plot_df.columns[0]

    # Ensure color column is numeric for scatter c=
    if color_col in plot_df.columns and not pd.api.types.is_numeric_dtype(plot_df[color_col]):
        converted = pd.to_numeric(plot_df[color_col], errors="coerce")
        if converted.notna().any():
            plot_df[color_col] = converted
        else:
            logger.warning(
                "Color column %s is not numeric. Using sequential values for coloring.", color_col
            )
            plot_df["_color_value"] = np.arange(len(plot_df), dtype=float)
            color_col = "_color_value"

    # For GSEA NES, prefer diverging colormap if user kept default
    if color_col == "NES" and cmap == "viridis_r":
        cmap = "RdBu_r"

    # Ensure x_col is numeric (user may pass x="FoldEnrichment" that is string/NA in some outputs)
    if x_col in plot_df.columns:
        plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")

    # Force size_col to numeric + fillna + all-NaN fallback (defensive for "12" strings etc.)
    plot_df[size_col] = pd.to_numeric(plot_df[size_col], errors="coerce")
    if plot_df[size_col].isna().all():
        logger.warning(
            "Size column %s is non-numeric or all missing. Using constant dot size.", size_col
        )
        plot_df["_dot_size"] = 1.0
        size_col = "_dot_size"
    else:
        plot_df[size_col] = plot_df[size_col].fillna(plot_df[size_col].median())

    # Drop rows that would be unplottable (after all the safe conversions)
    essential = [c for c in [x_col, color_col, size_col] if c in plot_df.columns]
    if essential:
        before = len(plot_df)
        plot_df = plot_df.dropna(subset=essential)
        if len(plot_df) == 0 and before > 0:
            logger.warning(
                "All rows dropped after requiring numeric x/size/color columns. Nothing to plot."
            )
            if ax is None:
                fig, ax = _empty_placeholder_fig("No valid terms after filtering")
                _created_fig = True
            else:
                fig = ax.figure
                _created_fig = False
            if _style_ctx is not None:
                _style_ctx.__exit__(None, None, None)
            _save_and_maybe_show(
                fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
            )
            return fig, ax

    # Robust dot size scaling (the root cause of "all dots look the same size").
    # Previous formula `val * 18 + 30` produced almost no visual difference
    # when the chosen size variable (Count or GeneRatio) had modest range
    # in the selected terms. We now do a proper min-max normalization to a
    # fixed, clearly visible marker area range. This guarantees that as long
    # as the underlying size_by values are not all identical, the dots will
    # have obviously different sizes.
    def _scale_sizes(
        vals,
        min_s=50,
        max_s=280,
        dot_max=None,
        dot_min=None,
        smallest_dot=0.0,
        vmin=None,
        vmax=None,
    ):
        vals = pd.to_numeric(vals, errors="coerce").astype(float)
        # Apply omicverse-style dot size limits before scaling
        if dot_max is not None:
            vals = np.minimum(vals, dot_max)
        if dot_min is not None:
            vals = np.maximum(vals, dot_min)
        if vmin is None:
            vmin = vals.min()
        if vmax is None:
            vmax = vals.max()
        if pd.isna(vmin) or pd.isna(vmax) or vmax <= vmin or len(vals) == 0:
            return np.full(len(vals), (min_s + max_s) / 2.0)
        sizes = min_s + (max_s - min_s) * (vals - vmin) / (vmax - vmin)
        if smallest_dot > 0:
            # Cleaner semantics:
            # smallest_dot=0 -> min point size = min_s (original behavior)
            # smallest_dot=1 -> all points use nearly max_s (everything "large")
            # Values in (0,1] raise the floor for the smallest dots proportionally.
            frac = (vals - vmin) / (vmax - vmin) if vmax > vmin else 0.0
            min_area = min_s + smallest_dot * (max_s - min_s)
            sizes = min_area + frac * (max_s - min_area)
        return np.clip(sizes, 20, 500)

    sizes = _scale_sizes(
        plot_df[size_col], dot_max=dot_max, dot_min=dot_min, smallest_dot=smallest_dot
    )

    if ax is None:
        # Use constrained_layout for robust automatic placement of colorbar + external size legend.
        # This is the approach used by high-quality implementations in gseapy and omicverse
        # and avoids the overlap issues that subplots_adjust + tight_layout combinations often cause.
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
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

    # Colorbar (p.adjust or equivalent) - compact placement on the right.
    # shrink/aspect/pad tuned following gseapy/omicverse-style dotplots for clean stacking with size legend.
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.28, aspect=12, pad=0.02)
    cbar_label = color_col
    if color_col in ("Adjusted P-value", "p.adjust"):
        cbar_label = "Adjusted P-value (smaller = more sig.)"
    cbar.set_label(cbar_label, fontsize=fontsize - 1, fontweight="bold", rotation=270, labelpad=18)
    cbar.outline.set_visible(False)

    # Size legend using proxy artists (keeps accurate representation of our custom _scale_sizes
    # including dot_max / smallest_dot controls, which is more reliable than raw legend_elements).
    #
    # Positioned in the upper-right (bbox_to_anchor + upper left) so it sits above / beside
    # the colorbar rather than fighting it horizontally. This + constrained_layout eliminates
    # the overlap between the two legend elements that was reported.
    #
    # References for this layout strategy:
    # - gseapy.plot.DotPlot.scatter + add_colorbar (bbox_to_anchor for size legend at ~ (1.02, 0.9))
    # - omicverse.pl.dotplot (clean right-side dual legend handling)
    try:
        size_vals_raw = pd.to_numeric(plot_df[size_col], errors="coerce").dropna().astype(float)
        if len(size_vals_raw) > 0:
            # Compute effective range after dot_max / dot_min (for correct scaling)
            effective = size_vals_raw.copy()
            if dot_max is not None:
                effective = np.minimum(effective, dot_max)
            if dot_min is not None:
                effective = np.maximum(effective, dot_min)
            eff_vmin = effective.min()
            eff_vmax = effective.max()

            # Choose representative values
            if size_col == "Count":
                # Prefer nice round numbers (multiples of 5/10) instead of raw min/median/max
                vmin_c = size_vals_raw.min()
                vmax_c = size_vals_raw.max()
                span = vmax_c - vmin_c
                if span <= 5:
                    step = 1
                elif span <= 15:
                    step = 5
                elif span <= 50:
                    step = 10
                else:
                    step = 20
                low = max(step, int(np.ceil(vmin_c / step) * step))
                high = int(np.floor(vmax_c / step) * step)
                mid = int(round((low + high) / 2 / step) * step)
                reps = []
                for cand in [low, mid, high]:
                    if vmin_c <= cand <= vmax_c:
                        reps.append(cand)
                if len(reps) < 3:
                    # fallback to a few nice values in range
                    reps = sorted({low, mid, high})
                    reps = [r for r in reps if vmin_c <= r <= vmax_c]
                if not reps:
                    reps = [int(round(vmin_c))]
                reps = sorted(set(reps))[:3]  # at most 3
            else:
                # For GeneRatio etc. keep behavior similar but use actual values
                reps = [size_vals_raw.min()]
                if len(size_vals_raw) > 2:
                    reps.append(size_vals_raw.median())
                reps.append(size_vals_raw.max())
                reps = sorted(set(reps))

            handles = []
            labels = []
            for rv in reps:
                # Scale using the global effective vmin/vmax of the plotted data
                # so that legend circle size is proportional and matches main plot dots
                s_for_rv = _scale_sizes(
                    pd.Series([rv]),
                    min_s=50,
                    max_s=280,
                    dot_max=dot_max,
                    dot_min=dot_min,
                    smallest_dot=smallest_dot,
                    vmin=eff_vmin,
                    vmax=eff_vmax,
                )[0]
                h = ax.scatter(
                    [], [], s=s_for_rv, c="#555555", alpha=0.7, edgecolors="#333333", linewidths=0.5
                )
                handles.append(h)
                if size_col == "Count":
                    labels.append(str(int(rv)))
                else:
                    labels.append(f"{rv:.2g}" if rv < 1 else f"{rv:.2f}")
            if handles:
                ax.legend(
                    handles,
                    labels,
                    title=size_col,
                    loc="upper left",
                    bbox_to_anchor=(1.02, 0.9),
                    frameon=False,
                    title_fontsize=fontsize - 1,
                    labelspacing=1.1,
                )
    except Exception:
        # Never let legend problems break the main plot
        pass

    sns.despine(ax=ax, top=True, right=True, left=False, bottom=False)

    # Do NOT call tight_layout() here. With constrained_layout + manual bbox_to_anchor legends
    # it frequently causes the colorbar and size legend to fight / overlap.
    # constrained_layout + bbox_inches="tight" on save (already done) + the upper-right placement
    # is the robust combination used by the referenced gseapy / omicverse implementations.

    if _style_ctx is not None:
        _style_ctx.__exit__(None, None, None)
    _save_and_maybe_show(
        fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
    )

    if return_data:
        return fig, ax, plot_df
    return fig, ax


def enrich_upsetplot(
    enrich_df,
    cluster_col: str | None = None,
    pval_cutoff: float = 0.05,
    min_count: int = 1,
    max_terms: int = 40,
    title: str = "Enriched Term Overlap (UpSet-style)",
    figsize: tuple[float, float] = (11, 6.5),
    dpi: int = 300,
    fontsize: int = 10,
    save_path: str | None = None,
    show: bool = True,
    use_style: bool = False,
):
    """
    UpSet-style plot for term overlap across clusters/groups (compareCluster style).

    Works directly on the output of:
      - scat.compare_enrichment(...)
      - scat.concat_compare_results(...)
      - or any df that has a "Cluster" (or cluster_col) + "Term"/"p.adjust"

    It shows:
    - Which enriched terms (padj < cutoff) are shared between different clusters/contrasts.
    - Set sizes (number of significant terms per cluster).
    - Intersection sizes (classic UpSet matrix + bars).

    This is especially powerful when you used `extract_gene_lists(..., separate_directions=True)`
    so that "up" and "down" from different contrasts appear as separate sets.

    If no 'Cluster' column is present, the whole table is treated as a single set (degrades
    gracefully to a simple bar of top terms).

    Parameters
    ----------
    cluster_col : str or None
        Column that defines the groups. Auto-detected if None.
    pval_cutoff : float
        Only terms with p.adjust (or similar) < cutoff are considered "enriched" for this plot.
    min_count : int
        Minimum gene count in the enrichment row to include the term.
    max_terms : int
        Cap on total unique terms considered (keeps plot readable).
    """
    _style_ctx = None
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()

    if enrich_df is None or (hasattr(enrich_df, "empty") and enrich_df.empty):
        fig, ax = _empty_placeholder_fig("No enrichment data for upset")
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        return fig, ax

    df = enrich_df.copy()

    # detect cluster
    if cluster_col is None:
        for cand in ("Cluster", "cluster", "Group", "group"):
            if cand in df.columns:
                cluster_col = cand
                break

    # p-value column
    p_col = None
    for c in ("p.adjust", "Adjusted P-value", "p_adj", "padj", "pvalue"):
        if c in df.columns:
            p_col = c
            break
    if p_col is None:
        p_col = df.columns[0]
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce").fillna(1.0)

    # count filter
    cnt_col = "Count" if "Count" in df.columns else None
    if cnt_col:
        df = df[pd.to_numeric(df[cnt_col], errors="coerce").fillna(0) >= min_count]

    # significant terms per cluster
    sig = df[df[p_col] < pval_cutoff].copy()
    if sig.empty:
        fig, ax = _empty_placeholder_fig("No significant terms at current cutoff")
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        if _style_ctx:
            _style_ctx.__exit__(None, None, None)
        return fig, ax

    # term column
    term_c = None
    for c in ["Term", "Description", "term", "description", "ID"]:
        if c in sig.columns:
            term_c = c
            break
    if term_c is None:
        term_c = sig.columns[0]

    if cluster_col and cluster_col in sig.columns:
        groups = sig[cluster_col].astype(str).unique().tolist()
        term_sets = {}
        for g in groups:
            terms = sig.loc[sig[cluster_col].astype(str) == g, term_c].astype(str).unique().tolist()
            term_sets[g] = set(terms)
    else:
        # single group fallback
        terms = sig[term_c].astype(str).unique().tolist()[:max_terms]
        term_sets = {"all": set(terms)}

    # limit total terms shown for clarity
    all_terms = set()
    for s in term_sets.values():
        all_terms.update(s)
    if len(all_terms) > max_terms:
        # keep the most frequent across sets
        from collections import Counter

        freq = Counter()
        for s in term_sets.values():
            for t in s:
                freq[t] += 1
        keep = {t for t, _ in freq.most_common(max_terms)}
        for g in list(term_sets.keys()):
            term_sets[g] = term_sets[g] & keep
        all_terms = keep

    if not all_terms:
        fig, ax = _empty_placeholder_fig("No terms after filtering")
        if _style_ctx:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax

    # Build membership matrix (for UpSet)
    term_list = sorted(all_terms)
    clusters_sorted = sorted(term_sets.keys())
    membership = pd.DataFrame(0, index=term_list, columns=clusters_sorted, dtype=int)
    for cl, terms in term_sets.items():
        for t in terms:
            if t in membership.index:
                membership.loc[t, cl] = 1

    # Intersection sizes
    # We group by the exact combination pattern
    combo = membership.apply(lambda r: tuple(r.values), axis=1)
    inter_sizes = combo.value_counts().sort_values(ascending=False)
    # Only show intersections with size >=1 (all do) and limit
    inter_sizes = inter_sizes.head(20)

    # Create figure layout
    fig = plt.figure(figsize=figsize, dpi=dpi)
    # Grid: left set-size bars, top intersection bars, matrix
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.6, 3.5],
        height_ratios=[1.2, 3.5],
        hspace=0.25,
        wspace=0.08,
    )
    ax_set = fig.add_subplot(gs[1, 0])  # set size horizontal bars (left)
    ax_inter = fig.add_subplot(gs[0, 1])  # intersection sizes (top)
    ax_mat = fig.add_subplot(gs[1, 1])  # matrix

    # 1. Set sizes (horizontal bars)
    set_sizes = {cl: len(term_sets[cl]) for cl in clusters_sorted}
    y_pos = np.arange(len(clusters_sorted))
    ax_set.barh(
        y_pos,
        [set_sizes[c] for c in clusters_sorted],
        color="#4477AA",
        edgecolor="black",
        height=0.6,
    )
    ax_set.set_yticks(y_pos)
    ax_set.set_yticklabels(clusters_sorted, fontsize=fontsize - 1)
    ax_set.invert_yaxis()
    ax_set.set_xlabel("Terms per group", fontsize=fontsize)
    ax_set.set_title("Set size", fontsize=fontsize + 1, fontweight="bold", pad=4)
    for i, v in enumerate([set_sizes[c] for c in clusters_sorted]):
        ax_set.text(v + 0.3, i, str(v), va="center", fontsize=fontsize - 2)
    sns.despine(ax=ax_set, top=True, right=True, left=False)

    # 2. Intersection size bars (top)
    inter_labels = []
    inter_vals = []
    for pattern, size in inter_sizes.items():
        # pattern is tuple of 0/1
        names = [clusters_sorted[i] for i, v in enumerate(pattern) if v == 1]
        inter_labels.append("\n".join(names) if names else "∅")
        inter_vals.append(size)
    x = np.arange(len(inter_vals))
    ax_inter.bar(x, inter_vals, color="#44AA77", edgecolor="black")
    ax_inter.set_ylabel("Intersection size", fontsize=fontsize)
    ax_inter.set_xticks([])
    ax_inter.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=6)
    for xi, v in zip(x, inter_vals):
        ax_inter.text(xi, v + 0.15, str(int(v)), ha="center", fontsize=fontsize - 2)
    sns.despine(ax=ax_inter, top=True, right=True)

    # 3. Matrix (lower right)
    ax_mat.set_xlim(-0.5, len(inter_sizes) - 0.5)
    ax_mat.set_ylim(-0.5, len(clusters_sorted) - 0.5)
    ax_mat.set_yticks(np.arange(len(clusters_sorted)))
    ax_mat.set_yticklabels(clusters_sorted, fontsize=fontsize - 1)
    ax_mat.set_xticks([])

    # Draw dots and lines for each shown intersection
    for col_idx, (pattern, _) in enumerate(inter_sizes.items()):
        active = [i for i, v in enumerate(pattern) if v == 1]
        # vertical line connecting active sets for this intersection
        if active:
            ax_mat.plot(
                [col_idx, col_idx], [min(active), max(active)], color="#555555", lw=1.5, zorder=1
            )
        for row_idx in range(len(clusters_sorted)):
            if pattern[row_idx] == 1:
                ax_mat.scatter(
                    col_idx,
                    row_idx,
                    s=140,
                    c="#CC3311",
                    zorder=2,
                    edgecolors="black",
                    linewidths=0.5,
                )
            else:
                ax_mat.scatter(
                    col_idx,
                    row_idx,
                    s=40,
                    c="#DDDDDD",
                    zorder=1,
                    edgecolors="#AAAAAA",
                    linewidths=0.3,
                )

    ax_mat.set_xlabel("Intersections (sorted by size)", fontsize=fontsize)
    ax_mat.invert_yaxis()
    sns.despine(ax=ax_mat, top=True, right=True, bottom=True, left=True)
    ax_mat.tick_params(left=False, bottom=False)

    plt.tight_layout()

    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
    return fig, ax_mat


def enrich_vennplot(
    enrich_df,
    cluster_col: str | None = None,
    pval_cutoff: float = 0.05,
    min_count: int = 1,
    max_terms: int = 200,
    title: str | None = None,
    figsize: tuple[float, float] = (8, 6),
    dpi: int = 300,
    fontsize: int = 11,
    colors: list | None = None,
    save_path: str | None = None,
    show: bool = True,
    use_style: bool = False,
):
    """
    Simple multi-group Venn diagram for significant enriched terms across clusters.

    Region labels show **exclusive** counts (elements belonging to exactly that
    overlap pattern, not cumulative pairwise totals). Supports 1–4 groups.

    Useful companion to enrich_upsetplot when you have 2-3 (max 4) groups and
    want a classic overlapping-circles view of shared vs unique terms.

    Works on the DataFrame returned by compare_enrichment / concat_compare_results
    (must have a Cluster column and p.adjust / Term).

    For >4 groups it gracefully degrades to a warning + the first 4.

    Pure matplotlib (no extra deps).
    """
    ctx = _style_context_if(use_style)
    with ctx:
        if enrich_df is None or (hasattr(enrich_df, "empty") and enrich_df.empty):
            fig, ax = _empty_placeholder_fig("No data for Venn")
            _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
            return fig, ax

        df = enrich_df.copy()

        if cluster_col is None:
            for cand in ("Cluster", "cluster", "Group"):
                if cand in df.columns:
                    cluster_col = cand
                    break
        if not cluster_col or cluster_col not in df.columns:
            # single set fallback
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
            ax.text(
                0.5, 0.5, "No Cluster column - nothing to Venn", ha="center", transform=ax.transAxes
            )
            ax.axis("off")
            _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
            return fig, ax

        # p col
        p_col = next(
            (c for c in ("p.adjust", "Adjusted P-value", "p_adj", "padj") if c in df.columns), None
        )
        if p_col is None:
            p_col = df.columns[0]
        df[p_col] = pd.to_numeric(df[p_col], errors="coerce").fillna(1)

        if "Count" in df.columns:
            df = df[pd.to_numeric(df["Count"], errors="coerce").fillna(0) >= min_count]

        sig = df[df[p_col] < float(pval_cutoff)].copy()
        if sig.empty:
            fig, ax = _empty_placeholder_fig("No significant terms")
            _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
            return fig, ax

        term_col = next(
            (c for c in ["Term", "Description", "term", "ID"] if c in sig.columns), sig.columns[0]
        )

        clusters = [str(x) for x in sig[cluster_col].dropna().unique().tolist()]
        if len(clusters) > 4:
            logger.warning("enrich_vennplot supports up to 4 groups; using first 4.")
            clusters = clusters[:4]

        term_sets = {}
        for cl in clusters:
            terms = set(sig.loc[sig[cluster_col].astype(str) == cl, term_col].astype(str).tolist())
            term_sets[cl] = terms

        all_terms = set.union(*term_sets.values()) if term_sets else set()
        if len(all_terms) > max_terms:
            logger.info("Limiting Venn to %d most frequent terms", max_terms)
            # keep most recurrent
            from collections import Counter

            freq = Counter(t for s in term_sets.values() for t in s)
            keep = {t for t, _ in freq.most_common(max_terms)}
            for k in list(term_sets):
                term_sets[k] = term_sets[k] & keep

        # Draw
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-1.5, 1.5)
        ax.set_aspect("equal")
        ax.axis("off")
        if title:
            ax.set_title(title, fontsize=fontsize + 2, fontweight="bold")

        n = len(clusters)
        if n == 0:
            ax.text(0.5, 0.5, "No groups", ha="center", transform=ax.transAxes)
            _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
            return fig, ax

        # Simple circle positions for 2-3-4 sets
        default_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
        cols = colors or default_colors[:n]

        positions = {
            1: [(0, 0, 0.9)],
            2: [(-0.35, 0.15, 0.75), (0.35, 0.15, 0.75)],
            3: [(-0.45, 0.25, 0.7), (0.45, 0.25, 0.7), (0.0, -0.45, 0.7)],
            4: [(-0.5, 0.4, 0.6), (0.5, 0.4, 0.6), (-0.5, -0.35, 0.6), (0.5, -0.35, 0.6)],
        }[n]

        from matplotlib.patches import Circle

        for (cx, cy, r), col, name in zip(positions, cols, clusters):
            circ = Circle((cx, cy), r, facecolor=col, alpha=0.25, edgecolor=col, linewidth=2)
            ax.add_patch(circ)
            ax.text(cx, cy + r + 0.12, name, ha="center", fontsize=fontsize - 1, fontweight="bold")

        # Exclusive region counts (each label = elements in exactly that overlap pattern)
        sets_list = [term_sets[c] for c in clusters]
        region_counts = _venn_exclusive_region_counts(sets_list)
        label_coords = _venn_region_label_coords(n)
        for key, count in region_counts.items():
            if count <= 0:
                continue
            pos = label_coords.get(key)
            if pos is None:
                continue
            region_size = len(key)
            fs = fontsize - (2 if region_size >= 3 else (1 if region_size == 2 else 0))
            weight = "bold" if region_size == n else "normal"
            ax.text(
                pos[0],
                pos[1],
                str(count),
                ha="center",
                va="center",
                fontsize=fs,
                fontweight=weight,
                color="#333",
            )

        ax.text(
            0.5,
            -1.35,
            f"Venn of significant enriched terms (padj < {pval_cutoff:.2g})",
            ha="center",
            transform=ax.transAxes,
            fontsize=fontsize - 2,
        )

        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax


def gseaplot(
    ranked_genes: pd.Series | Mapping | Iterable,
    gsea_result: pd.DataFrame | None = None,
    term: str | None = None,
    title: str | None = None,
    figsize: tuple[float, float] = (6.5, 5.5),
    dpi: int = 300,
    color: str = "#88C544",
    cmap: str = "seismic",
    ax: Any | None = None,
    show: bool = True,
    use_style: bool = False,
    save_path: str | None = None,
    pheno_pos: str = "Pos",
    pheno_neg: str = "Neg",
    **kwargs: Any,
) -> tuple[Any | None, Any | None]:
    """
    Classic GSEA plot: running enrichment score (RES) curve + hits + ranked list.

    Designed to work seamlessly with results from `scat.run_gsea`.

    Parameters
    ----------
    ranked_genes : pd.Series or list-like
        Gene scores (index=genes, values=scores). Same as passed to run_gsea.
    gsea_result : pd.DataFrame, optional
        Result from run_gsea(). If provided with `term`, will try to use pre-computed
        RES curve and hits stored in .attrs['gsea_details'] for pixel-perfect match.
    term : str, optional
        Term name to plot (must match a row in gsea_result["Term"]).
    title, figsize, dpi, save_path, show, use_style, ax : standard scATrans plot args.
    color : color for the RES curve and hit ticks.
    cmap : colormap for the ranked list bar at bottom.
    pheno_pos, pheno_neg : labels for phenotype (shown in plot).
    """
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()
    else:
        _style_ctx = None

    # Normalize ranked input to Series
    if isinstance(ranked_genes, (list, tuple)):
        # assume order high->low, make scores
        ranked = pd.Series(range(len(ranked_genes), 0, -1), index=[str(x) for x in ranked_genes])
    elif isinstance(ranked_genes, Mapping):
        ranked = pd.Series(ranked_genes)
    else:
        ranked = pd.Series(ranked_genes)

    ranked = ranked.dropna()
    ranked.index = ranked.index.astype(str)
    if len(ranked) == 0:
        logger.warning("No ranked genes for gseaplot")
        fig, ax = _empty_placeholder_fig("No ranked genes")
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        return fig, ax

    # Try to get precomputed data from gsea_result.attrs
    RES = None
    hits = None
    nes = np.nan
    pval = np.nan
    fdr = np.nan
    used_term = term

    if gsea_result is not None and term is not None:
        df = gsea_result
        # find the row
        mask = (
            df["Term"].astype(str).str.strip().str.lower() == str(term).strip().lower()
            if "Term" in df.columns
            else pd.Series(False, index=df.index)
        )
        if mask.any():
            row = df.loc[mask].iloc[0]
            if "NES" in row:
                nes = row.get("NES", np.nan)
            if "pvalue" in row or "NOM p-val" in row:
                pval = row.get("pvalue", row.get("NOM p-val", np.nan))
            if "p.adjust" in row or "FDR q-val" in row:
                fdr = row.get("p.adjust", row.get("FDR q-val", np.nan))
            used_term = row.get("Term", term)

        # precomputed from our run_gsea
        if "gsea_details" in getattr(df, "attrs", {}):
            details = df.attrs["gsea_details"]
            key = None
            for k in details:
                if str(k).strip().lower() == str(term).strip().lower():
                    key = k
                    break
            if key and isinstance(details[key], dict):
                d = details[key]
                RES = d.get("RES")
                hits = d.get("hits", [])
                if "nes" in d:
                    nes = d["nes"]
                if "pval" in d:
                    pval = d["pval"]
                if "fdr" in d:
                    fdr = d["fdr"]

    # (optional future: recompute via gseapy if needed)

    # If still no RES, compute a basic RES (approximation; for exact use run_gsea + stored)
    if RES is None:
        # Try to recover genes for the term
        gene_set = set()
        if gsea_result is not None and term is not None and "leading_edge" in gsea_result.columns:
            try:
                lead = gsea_result.loc[
                    gsea_result["Term"].astype(str).str.lower() == str(term).lower(), "leading_edge"
                ].iloc[0]
                if isinstance(lead, str):
                    gene_set = set(lead.split(";"))
                elif isinstance(lead, (list, tuple)):
                    gene_set = set(lead)
            except Exception:
                pass
        if not gene_set:
            # last resort: empty plot
            logger.warning("Could not recover gene set or RES for term %s. Plotting empty.", term)
            fig, ax = _empty_placeholder_fig(f"No data for term {term}")
            if _style_ctx is not None:
                _style_ctx.__exit__(None, None, None)
            return fig, ax

        # very basic running sum (for demo; prefer precomputed)
        N = len(ranked)
        k = len(gene_set & set(ranked.index))
        if k == 0:
            RES = [0.0] * N
            hits = []
        else:
            weights = np.abs(ranked.values)
            sum_hit = weights[[i for i, g in enumerate(ranked.index) if g in gene_set]].sum() or 1.0
            miss = 1.0 / (N - k) if k < N else 0.0
            running = 0.0
            RES = []
            hits = []
            for i, g in enumerate(ranked.index):
                if g in gene_set:
                    running += weights[i] / sum_hit
                    hits.append(i)
                else:
                    running -= miss
                RES.append(running)
    # Now plot
    if ax is not None:
        logger.warning(
            "gseaplot creates a 3-panel figure and does not support single `ax`. "
            "Returning a new figure. For embedding, pass `axes=(ax1, ax2, ax3)` in future or use subplots externally."
        )
    if ax is None:
        fig, axes = plt.subplots(
            3,
            1,
            figsize=figsize,
            dpi=dpi,
            gridspec_kw={"height_ratios": [0.5, 0.1, 0.4], "hspace": 0.05},
            sharex=True,
        )
        ax1, ax2, ax3 = axes  # RES, hits, ranked
        _created = True
    else:
        fig, axes = plt.subplots(3, 1, figsize=figsize, dpi=dpi, sharex=True)
        ax1, ax2, ax3 = axes
        _created = True

    # 1. RES curve
    x = np.arange(len(RES))
    ax1.plot(x, RES, color=color, linewidth=1.5, label="Running ES")
    ax1.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax1.fill_between(x, RES, 0, alpha=0.2, color=color)
    ax1.set_ylabel("Enrichment Score (ES)", fontsize=9)
    ax1.set_title(title or (used_term or "GSEA Plot"), fontsize=11, fontweight="bold")
    ax1.legend(loc="upper left", frameon=False, fontsize=8)

    # annotate
    if not np.isnan(nes):
        ax1.text(
            0.98,
            0.95,
            f"NES={nes:.2f}\n p={pval:.3g}\n FDR={fdr:.3g}",
            transform=ax1.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
        )

    # 2. Hits ticks
    ax2.set_ylim(-0.5, 0.5)
    if hits:
        ax2.vlines(hits, -0.3, 0.3, colors=color, linewidth=0.8)
    ax2.set_yticks([])
    ax2.set_ylabel("Hits", fontsize=8, rotation=0, ha="right", va="center")
    sns.despine(ax=ax2, top=True, right=True, left=True, bottom=True)

    # 3. Ranked metric (bar or scatter)
    rank_vals = ranked.values
    # color by sign or value
    norm = Normalize(vmin=min(rank_vals), vmax=max(rank_vals))
    colors = plt.get_cmap(cmap)(norm(rank_vals))
    ax3.bar(x, rank_vals, color=colors, width=1.0, edgecolor="none")
    ax3.axhline(0, color="black", linewidth=0.5)
    ax3.set_ylabel("Ranked\nmetric", fontsize=8)
    ax3.set_xlabel("Rank in Ordered Dataset", fontsize=9)
    ax3.set_xlim(-0.5, len(x) - 0.5)

    # labels for pos/neg
    ax3.text(0.02, 0.95, pheno_pos, transform=ax3.transAxes, fontsize=8, color="red")
    ax3.text(0.98, 0.05, pheno_neg, transform=ax3.transAxes, ha="right", fontsize=8, color="blue")

    sns.despine(ax=ax3, top=True, right=True)

    if _style_ctx is not None:
        _style_ctx.__exit__(None, None, None)

    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=_created)
    # return the main ax (RES) for consistency, or the fig and list
    return fig, ax1


def volcano_plot(
    df,
    top_n=10,
    label_genes: Iterable[str] | None = None,
    save_path=None,
    title="Volcano Plot of Active Transcription",
    point_scale=1.0,
    min_size=2,
    max_size=160,
    s: float
    | None = None,  # fixed point size (overrides variable sizing by score or pval); direct control like in omicverse.pl
    alpha: float = 0.75,
    figsize=(8, 6),
    dpi=300,
    fontsize=12,
    cmap="coolwarm",
    logfc_cutoff=0.5,
    pval_cutoff=0.05,
    color_by="active_score",
    ax=None,
    show: bool = True,
    use_style: bool = False,
    return_data: bool = False,
    label_repel: bool = True,
    label_fontsize: float | None = None,
    min_label_score: float | None = None,
    close: bool = False,
):
    """
    2D volcano plot with ggVolcano-inspired flexibility and style options.

    - `top_n`: number of top genes (by active_score or p_adj) to label.
    - `label_genes`: iterable of gene names to *force* label (manual specification,
      even if not in top_n). This + top_n gives the common ggVolcano usage
      pattern (label_number + explicit genes). Duplicates are handled automatically.
    - `color_by`: if set to "active_score" (default) and the column exists, uses continuous
      colormap. If set to any other column present in the data, that column is used for
      continuous coloring (after numeric coercion). Otherwise falls back to classic
      up/down/ns significance categories (red/blue/gray) based on logfc_cutoff / pval_cutoff.
    - Cutoff lines are drawn with labels.
    - Supports `ax` for embedding.

    Size control (addresses "even smallest points are still too big"; modeled after flexible controls in omicverse.pl.volcano etc.):
      - `s`: if given, forces a **fixed** point size for all scatters (e.g. s=2 or s=1 for tiny points).
        This bypasses score-based variable sizing entirely — simplest for clean small-point volcanos.
      - `point_scale`: global multiplier.
      - `min_size` / `max_size`: bounds for the variable sizing case.
      - When no "active_score" (pure DE volcano), uses small p-value based base so points start small.

    Style reference: https://github.com/BioSenior/ggVolcano (label control,
    clean up/down distinction, readable labels with repel).
    """
    logger.info("Generating 2D volcano plot...")
    _style_ctx = None
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()

    fig_g, ax_g, guarded = _guard_plot_dataframe(
        df, ["logFC", "p_adj"], "volcano_plot", ax=ax, figsize=figsize, dpi=dpi
    )
    if guarded:
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig_g, save_path=save_path, dpi=dpi, show=show, created=ax is None, close=close
        )
        return fig_g, ax_g

    _require_columns(df, ["logFC", "p_adj"], "volcano_plot")
    if pval_cutoff <= 0 or pval_cutoff >= 1:
        raise ValueError("pval_cutoff must be between 0 and 1 (exclusive).")
    if logfc_cutoff < 0:
        raise ValueError("logfc_cutoff must be non-negative.")

    plot_df = df.copy()
    for c in ["logFC", "p_adj", "velocity_residual", "active_score"]:
        if c in plot_df.columns:
            plot_df[c] = pd.to_numeric(plot_df[c], errors="coerce")

    # Filter invalid p_adj ( <0 or >1 are semantically invalid for adjusted p-values)
    if "p_adj" in plot_df.columns:
        invalid_p = (plot_df["p_adj"] < 0) | (plot_df["p_adj"] > 1)
        if invalid_p.any():
            logger.warning("Dropping %d rows with p_adj outside [0, 1].", int(invalid_p.sum()))
            plot_df = plot_df.loc[~invalid_p].copy()

    plot_df = plot_df.dropna(subset=["logFC", "p_adj"])
    plot_df["neg_log_pval"] = _safe_neg_log10(plot_df["p_adj"])
    plot_df = plot_df.dropna(subset=["neg_log_pval"])  # in case p_adj produced NaN after coercion

    if plot_df.empty:
        logger.warning("No valid genes to plot after numeric coercion/dropna in volcano_plot.")
        if ax is None:
            fig, ax = _empty_placeholder_fig("No valid genes to plot")
            _created_fig = True
        else:
            fig = ax.figure
            _created_fig = False
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(
            fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
        )
        return fig, ax

    # ggVolcano-style classic coloring (up/down/ns) when not using active_score
    if color_by == "active_score" and "active_score" in plot_df.columns:
        color_values = plot_df["active_score"]
        cbar_label = "Active Score"
        colors_for_scatter = None
    elif color_by in plot_df.columns:
        # User requested a custom column for color (e.g. some other score)
        plot_df[color_by] = pd.to_numeric(plot_df[color_by], errors="coerce")
        if plot_df[color_by].notna().any():
            color_values = plot_df[color_by]
            cbar_label = color_by
            colors_for_scatter = None
        else:
            logger.warning(
                "color_by=%s present but non-numeric after coercion. Falling back to significance categories.",
                color_by,
            )
            # fall through to classic significance coloring
            up_mask = (plot_df["logFC"] > logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
            down_mask = (plot_df["logFC"] < -logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
            color_values = np.where(up_mask, 2, np.where(down_mask, 1, 0))
            cbar_label = None
            colors_for_scatter = ["#808080", "#1f77b4", "#d62728"]
    else:
        if color_by != "active_score":
            logger.warning(
                "color_by=%s not found in data. Falling back to up/down/ns significance categories.",
                color_by,
            )
        up_mask = (plot_df["logFC"] > logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
        down_mask = (plot_df["logFC"] < -logfc_cutoff) & (plot_df["p_adj"] < pval_cutoff)
        color_values = np.where(up_mask, 2, np.where(down_mask, 1, 0))  # 2=up, 1=down, 0=ns
        cbar_label = None
        # Use explicit nice colors similar to common ggVolcano / EnhancedVolcano
        colors_for_scatter = ["#808080", "#1f77b4", "#d62728"]  # ns, down, up

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
        _created_fig = True
    else:
        fig = ax.figure
        _created_fig = False

    # Size logic (omicverse-inspired direct `s` control + our variable + min/max).
    if s is not None:
        sizes = np.full(len(plot_df), float(s) * point_scale)
        effective_min = min(1.0, min_size)
        sizes = np.clip(sizes, effective_min, max_size)
    else:
        # variable
        # Defensively clip active_score (or fallback) to >=0 before exponentiation.
        if "active_score" in plot_df.columns:
            size_val = np.clip(
                pd.to_numeric(plot_df["active_score"], errors="coerce").fillna(0), 0, None
            )
        else:
            size_val = plot_df.get("neg_log_pval", pd.Series(4, index=plot_df.index))
        raw_sizes = size_val**1.3 * 8 * point_scale + 3 * point_scale
        sizes = np.clip(raw_sizes, min_size, max_size)
    sizes = _finite_scatter_values(sizes, fill=float(min_size))

    # Light diagnostic (omicverse style)
    if len(plot_df) > 1000 and (s is None) and point_scale > 0.25:
        logger.info(
            "Large number of points (%d) in volcano. For cleaner view consider "
            "s=2 (fixed small) or point_scale<=0.15 + min_size=1.",
            len(plot_df),
        )

    scatter_kwargs = {
        "x": plot_df["logFC"],
        "y": plot_df["neg_log_pval"],
        "s": sizes,
        "alpha": alpha,
        "edgecolors": "#444444",
        "linewidth": 0.4,
        "zorder": 3,
    }
    if colors_for_scatter is not None:
        # Classic up/down/ns: provide explicit color list (do not pass c= together with color=)
        scatter = ax.scatter(c=[colors_for_scatter[int(c)] for c in color_values], **scatter_kwargs)
    else:
        scatter_kwargs["c"] = _finite_scatter_values(color_values, fill=0.0)
        scatter = ax.scatter(cmap=cmap, **scatter_kwargs)

    ax.axhline(
        float(_safe_neg_log10(pval_cutoff)),
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
    label_pool = plot_df
    if min_label_score is not None and "active_score" in label_pool.columns:
        label_pool = label_pool[label_pool["active_score"] >= float(min_label_score)]
    if "active_score" in label_pool.columns:
        top_df = label_pool.nlargest(top_n, "active_score")
    else:
        top_df = label_pool.nsmallest(top_n, "p_adj")
    for g in top_df.index:
        genes_to_label.add(str(g))

    label_df = (
        plot_df.loc[plot_df.index.astype(str).isin(genes_to_label)].copy()
        if genes_to_label
        else pd.DataFrame()
    )

    lbl_fs = float(label_fontsize) if label_fontsize is not None else max(8, fontsize - 2)
    texts = []
    for idx, row in label_df.iterrows():
        txt = ax.text(
            row["logFC"],
            row["neg_log_pval"],
            str(idx),
            fontsize=lbl_fs,
            fontweight="bold",
            color="#111111",
            bbox={"boxstyle": "round,pad=0.2", "fc": "white", "ec": "none", "alpha": 0.75},
        )
        texts.append(txt)

    _maybe_repel_labels(
        texts,
        plot_df["logFC"].values,
        plot_df["neg_log_pval"].values,
        ax,
        label_repel=label_repel,
    )

    ax.set_xlabel("Log2 Fold Change", fontsize=fontsize, fontweight="bold")
    ax.set_ylabel("-Log10(adj. P-value)", fontsize=fontsize, fontweight="bold")
    if title:
        ax.set_title(title, fontsize=fontsize + 2, fontweight="bold", pad=15)

    if colors_for_scatter is None:
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.6, pad=0.02, aspect=20)
        cbar.set_label(
            cbar_label or color_by,
            fontsize=max(9, fontsize - 1),
            fontweight="bold",
            rotation=270,
            labelpad=15,
        )
        cbar.outline.set_visible(False)
    else:
        ax.legend(loc="upper left", frameon=False, fontsize=fontsize - 1)

    sns.despine(ax=ax, top=True, right=True)
    ax.spines["left"].set_position(("outward", 6))
    ax.spines["bottom"].set_position(("outward", 6))

    # Note: constrained_layout=True at creation handles colorbar placement cleanly.
    # Avoid tight_layout after colorbar/legends.

    if _style_ctx is not None:
        _style_ctx.__exit__(None, None, None)
    _save_and_maybe_show(
        fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig, close=close
    )
    if return_data:
        return fig, ax, plot_df
    return fig, ax


def bias_diagnostic_plot(
    results_df,
    save_path=None,
    title="Bias Correction Diagnostic",
    figsize=(12, 5),
    dpi=300,
    fontsize=11,
    show_regression=True,
    point_size=10,
    axes=None,
    show: bool = True,
    use_style: bool = False,
):
    """
    Diagnostic plot showing the effect of gene length / intron number bias
    correction on velocity delta (before vs after).

    Supports external `axes` (tuple of two Axes) for embedding in custom figures.

    `point_size`: size for the background gene cloud (default 10, was 15).
    """
    logger.info("Generating bias correction diagnostic plot...")
    if use_style:
        set_style()

    delta_col = _excess_delta_col(results_df)
    residual_col = _excess_residual_col(results_df)
    _require_columns(
        results_df,
        [delta_col, residual_col, "gene_length", "intron_number"],
        "bias_diagnostic_plot",
    )

    plot_df = results_df.dropna(
        subset=[delta_col, residual_col, "gene_length", "intron_number"]
    ).copy()
    if len(plot_df) < 10:
        logger.warning("Too few genes with valid features for diagnostic plot.")
        # Create a 1x2 placeholder figure to match normal return type (fig, (ax1, ax2))
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
        for a in (ax1, ax2):
            a.text(
                0.5,
                0.5,
                "Too few genes for bias diagnostic",
                ha="center",
                va="center",
                transform=a.transAxes,
                fontsize=9,
            )
            a.axis("off")
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, (ax1, ax2)

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=figsize, dpi=dpi, constrained_layout=True)
        _created_fig = True
        ax1 = axes[0]
        ax2 = axes[1]
    else:
        from matplotlib.axes import Axes

        if len(axes) != 2 or not all(isinstance(a, Axes) for a in axes):
            raise ValueError(
                "axes must be a sequence of exactly two matplotlib Axes instances "
                "(before-panel, after-panel). Omit axes= to create a new 1×2 figure."
            )
        fig = axes[0].figure
        _created_fig = False
        ax1 = axes[0]
        ax2 = axes[1]

    # Left: Before correction
    x = np.log1p(plot_df["gene_length"])
    y_raw = plot_df[delta_col]
    ax1.scatter(x, y_raw, s=point_size, alpha=0.5, c="#1f77b4", edgecolors="none")
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
    ax1.set_ylabel("Unspliced excess delta (raw)", fontsize=fontsize, fontweight="bold")
    ax1.set_title("Before Bias Correction", fontsize=fontsize + 1, fontweight="bold")
    ax1.legend(frameon=False)
    sns.despine(ax=ax1)

    # Right: After correction
    y_res = plot_df[residual_col]
    ax2.scatter(x, y_res, s=point_size, alpha=0.5, c="#2ca02c", edgecolors="none")
    ax2.axhline(0, color="#d62728", linestyle="--", lw=1.2, alpha=0.8)
    ax2.set_xlabel("log1p(Gene Length)", fontsize=fontsize, fontweight="bold")
    ax2.set_ylabel(
        "Unspliced excess residual (bias-corrected)", fontsize=fontsize, fontweight="bold"
    )
    ax2.set_title("After Bias Correction", fontsize=fontsize + 1, fontweight="bold")
    sns.despine(ax=ax2)

    if title:
        fig.suptitle(title, fontsize=fontsize + 2, fontweight="bold", y=1.02)

    # constrained_layout handles the two-panel + suptitle cleanly.

    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=_created_fig)
    return fig, axes


# =============================================================================
# Additional / legacy plotting helpers
# =============================================================================


def enrich_barplot(enrich_df, top_n=15, title="Enrichment Barplot", save_path=None, **kwargs):
    """Barplot wrapper around the dotplot implementation (for API compatibility).

    NOTE: this is currently a thin alias that calls enrich_dotplot (visual is a dotplot).
    For a true bar implementation, call a custom bar or extend this function.
    """
    logger.warning(
        "enrich_barplot is deprecated (it is an alias for enrich_dotplot and does not produce a barplot). "
        "Use enrich_dotplot directly; this alias may be removed in a future version."
    )
    # pass through new standard kwargs if caller used them before they existed on the alias
    return enrich_dotplot(enrich_df, top_n=top_n, title=title, save_path=save_path, **kwargs)


def active_score_rankplot(
    results_df,
    top_n=20,
    save_path=None,
    ax=None,
    dpi=300,
    show: bool = True,
    use_style: bool = False,
    **kwargs,
):
    """
    Horizontal ranked barplot of top active scores (publication-friendly).

    Improvements inspired by omicverse bulk/perturbation ranked visualizations:
    - Gradient coloring by active_score magnitude.
    - Clean outward-offset spines.
    - Value annotations on bars.
    - Good `ax=` embedding support and constrained layout.

    For richer single-gene context prefer `pl.comet_plot` or `pl.volcano_plot`.
    """
    logger.info("Generating active score rank plot...")
    if use_style:
        set_style()

    if results_df is None or results_df.empty:
        logger.warning("No results to plot.")
        fig, ax = _empty_placeholder_fig("No results to plot")
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax

    _require_columns(results_df, ["active_score"], "active_score_rankplot")

    plot_df = results_df.copy()
    plot_df["active_score"] = pd.to_numeric(plot_df["active_score"], errors="coerce")
    plot_df = plot_df.dropna(subset=["active_score"])

    if plot_df.empty:
        logger.warning("No valid numeric active_score values to plot.")
        fig, ax = _empty_placeholder_fig("No valid active_score values")
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax

    # Use nlargest for safety: does not assume the input df is already sorted by active_score.
    plot_df = plot_df.nlargest(top_n, "active_score").iloc[::-1]  # top at top for horizontal bar

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, max(4, 0.38 * top_n)), dpi=dpi, constrained_layout=True)
        _created = True
    else:
        fig = ax.figure
        _created = False

    # Gradient coloring by score (omicverse-style ranked emphasis)
    vmin = float(plot_df["active_score"].min())
    vmax = float(plot_df["active_score"].max())
    if vmin == vmax:
        vmax = vmin + 1e-9  # avoid singular Normalize when all scores identical
    norm = Normalize(vmin=vmin, vmax=vmax)
    colors = [cm.viridis(v) for v in norm(plot_df["active_score"].values)]

    bars = sns.barplot(
        data=plot_df,
        y=plot_df.index,
        x="active_score",
        ax=ax,
        palette=colors,
        hue=plot_df.index,
        legend=False,
        edgecolor="#333333",
        linewidth=0.5,
    )

    # Add value labels on bars (clean, small). Use data range for offset (works for negative scores too).
    score_min = plot_df["active_score"].min()
    score_max = plot_df["active_score"].max()
    x_range = (score_max - score_min) if score_max > score_min else 1.0
    offset = 0.02 * x_range
    for bar, val in zip(bars.patches, plot_df["active_score"]):
        if val >= 0:
            x = val + offset
            ha = "left"
        else:
            x = val - offset
            ha = "right"
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.2f}",
            va="center",
            ha=ha,
            fontsize=8,
            color="#222222",
        )

    ax.set_xlabel("Active Score", fontweight="bold")
    ax.set_ylabel("")
    ax.set_title("Top Active Drivers (rank)", fontweight="bold", pad=10)

    # Outward spines (omicverse bulk style)
    sns.despine(ax=ax, top=True, right=True)
    ax.spines["left"].set_position(("outward", 8))
    ax.spines["bottom"].set_position(("outward", 8))

    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=_created)
    return fig, ax


def active_genes_heatmap(
    adata,
    genes=None,
    groupby=None,
    save_path=None,
    show: bool = True,
    use_style: bool = False,
    **kwargs,
):
    """
    Convenience wrapper around scanpy heatmap for the active driver genes.

    Users are encouraged to call scanpy.pl.heatmap directly with the genes
    returned by active_score for full control.

    show/use_style are accepted for API consistency with other pl.* functions.
    """
    if use_style:
        set_style()

    if genes is None:
        # try to guess from var
        if "active_score" in adata.var.columns:
            genes = adata.var.nlargest(20, "active_score").index.tolist()
        else:
            logger.warning("No genes provided and no active_score column found.")
            fig, ax = _empty_placeholder_fig("No active genes to heatmap")
            return fig, ax

    logger.info(
        "active_genes_heatmap: delegating to scanpy.pl.heatmap (recommended for full control)"
    )
    try:
        sc = __import__("scanpy", fromlist=["pl"])
        fig = sc.pl.heatmap(
            adata,
            var_names=genes,
            groupby=groupby,
            show=False,  # we control final show below for consistency
            save=None,
            **kwargs,
        )
        if save_path:
            from pathlib import Path

            out = Path(str(save_path))
            out.parent.mkdir(parents=True, exist_ok=True)
            import matplotlib.pyplot as plt

            plt.savefig(out, dpi=300, bbox_inches="tight")
            logger.info("Heatmap saved → %s", out)
        if show:
            plt.show()
        return fig, None
    except Exception as e:
        logger.warning("active_genes_heatmap could not render via scanpy: %s", e)
        fig, ax = _empty_placeholder_fig("active_genes_heatmap failed")
        return fig, ax


def set_nature_style(**kwargs):
    """Legacy alias for set_style() kept for backward compatibility.

    Accepts the same kwargs as set_style (font sizes, dpi etc).
    """
    return set_style(**kwargs)


def velocity_phase_portraits(
    adata,
    genes,
    groupby=None,
    spliced_layer="spliced",
    unspliced_layer="unspliced",
    max_genes=6,
    figsize_per_gene=(2.8, 2.4),
    save_path=None,
    dpi=300,
    show: bool = True,
    use_style: bool = False,
    **kwargs,
):
    """
    Quick diagnostic grid of unspliced vs spliced (phase-portrait style) for selected genes.

    Useful for visually inspecting whether top active genes show the expected excess
    nascent RNA in the target group. Points are colored by the groupby column when provided.

    This is intentionally lightweight — for full control users are encouraged to write
    their own small U/S scatter functions.

    Supports `dpi`, `show`, `use_style` for consistency with other pl.* functions.
    When groupby is provided, uses a categorical colormap (tab20) + figure-level legend.

    Parameters
    ----------
    genes : list-like
        Gene names (index of adata.var) to plot.
    groupby : str, optional
        obs column used for coloring (e.g. the same contrast column used in active_score).
    max_genes : int
        Maximum number of genes to plot (grid will be truncated).
    dpi : int
        DPI for figure creation and save (default 300 for publication quality).
    show : bool
        Whether to call plt.show() (only relevant when function created the figure).
    """
    import math

    genes = list(genes)[:max_genes]
    if not genes:
        logger.warning("No genes provided for phase portraits.")
        fig, ax = _empty_placeholder_fig("No genes for phase portraits", figsize=(2, 2))
        return fig, [ax]  # minimal list for axes compat; callers should check

    if use_style:
        set_style()

    # Validate layers early (clear error instead of opaque KeyError later)
    missing = [lyr for lyr in (spliced_layer, unspliced_layer) if lyr not in adata.layers]
    if missing:
        raise KeyError(
            f"velocity_phase_portraits: missing layers in adata.layers: {missing}. "
            f"Available layers: {list(adata.layers.keys())}"
        )

    n = len(genes)
    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * figsize_per_gene[0], nrows * figsize_per_gene[1]),
        dpi=dpi,
        squeeze=False,
        constrained_layout=True,
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
        # Robust to sparse, numpy matrix, 2D ndarray etc: always ensure 1D float array
        if sparse.issparse(u):
            u = u.toarray()
        if sparse.issparse(s):
            s = s.toarray()
        u = np.asarray(u).ravel()
        s = np.asarray(s).ravel()

        color = None
        if groupby:
            if groupby not in adata.obs:
                logger.warning(
                    "groupby='%s' not found in adata.obs. Available obs columns: %s. "
                    "Falling back to uniform color.",
                    groupby,
                    list(adata.obs.columns),
                )
            else:
                # Convert labels to numeric codes for scatter c= (avoids matplotlib error on string arrays)
                try:
                    color = pd.Categorical(adata.obs[groupby]).codes
                except Exception:
                    color = None

        if color is not None and np.asarray(color).size == len(s):
            ax.scatter(s, u, c=color, cmap="tab20", s=8, alpha=0.6, edgecolors="none")
        else:
            ax.scatter(s, u, c="#2ca02c", s=8, alpha=0.6, edgecolors="none")
        ax.set_xlabel("Spliced", fontsize=9)
        ax.set_ylabel("Unspliced", fontsize=9)
        ax.set_title(str(g), fontsize=10, fontweight="bold")
        sns.despine(ax=ax)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    # Add a figure-level legend when coloring by groupby (outside the subplots grid)
    if groupby and groupby in adata.obs:
        try:
            from matplotlib.lines import Line2D

            cat = pd.Categorical(adata.obs[groupby].astype(str))
            groups = list(cat.categories)
            if len(groups) > 20:
                logger.warning(
                    "velocity_phase_portraits: %d groups for '%s'; legend will show only the first 20.",
                    len(groups),
                    groupby,
                )
            cmap = plt.get_cmap("tab20")
            handles = [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=cmap(i % 20),
                    markersize=5,
                    linestyle="None",
                    label=str(g),
                )
                for i, g in enumerate(groups[: min(20, len(groups))])
            ]
            if handles:
                fig.legend(
                    handles,
                    [h.get_label() for h in handles],
                    loc="upper right",
                    frameon=False,
                    fontsize=8,
                    title=str(groupby),
                    bbox_to_anchor=(0.98, 0.98),
                )
        except Exception:
            # Legend is nice-to-have; never let it break the phase portraits
            pass

    # constrained_layout used at creation for the gene grid.

    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
    return fig, axes


def gamma_shrinkage_plot(
    results_df,
    *,
    x_col: str = "total_us_counts",
    save_path=None,
    title: str = "Empirical Bayes gamma shrinkage",
    figsize=(7, 5),
    dpi: int = 300,
    ax=None,
    show: bool = True,
    use_style: bool = False,
):
    """
    Visualize per-gene shrinkage strength vs expression depth (empirical Bayes gamma).

    Requires ``gamma_shrinkage_weight`` in the results table (produced when
    ``gamma_method='empirical_bayes'``). Optionally colors by ``effective_gamma``
    when present.

    Returns
    -------
    (fig, ax)
    """
    logger.info("Generating gamma shrinkage plot...")
    _style_ctx = None
    if use_style:
        _style_ctx = style_context()
        _style_ctx.__enter__()

    if results_df is None or results_df.empty:
        fig, ax = _empty_placeholder_fig("No results to plot")
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax

    _require_columns(results_df, ["gamma_shrinkage_weight"], "gamma_shrinkage_plot")

    plot_df = results_df.copy()
    plot_df["gamma_shrinkage_weight"] = pd.to_numeric(
        plot_df["gamma_shrinkage_weight"], errors="coerce"
    )
    if x_col in plot_df.columns:
        plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    else:
        x_col = None

    plot_df = plot_df.dropna(subset=["gamma_shrinkage_weight"])
    if plot_df.empty:
        fig, ax = _empty_placeholder_fig("No valid gamma_shrinkage_weight values")
        if _style_ctx is not None:
            _style_ctx.__exit__(None, None, None)
        _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=True)
        return fig, ax

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
        _created = True
    else:
        fig = ax.figure
        _created = False

    x_vals = np.log1p(plot_df[x_col].values) if x_col is not None else np.arange(len(plot_df))
    y_vals = plot_df["gamma_shrinkage_weight"].values

    c_vals = None
    if "effective_gamma" in plot_df.columns:
        c_vals = pd.to_numeric(plot_df["effective_gamma"], errors="coerce").values

    scatter = ax.scatter(
        x_vals,
        y_vals,
        c=c_vals,
        cmap="viridis" if c_vals is not None else None,
        s=12,
        alpha=0.65,
        edgecolors="none",
    )
    if c_vals is not None:
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.75, pad=0.02)
        cbar.set_label("effective_gamma", fontsize=9)

    ax.set_xlabel(
        f"log1p({x_col})" if x_col is not None else "gene index",
        fontsize=11,
        fontweight="bold",
    )
    ax.set_ylabel("Shrinkage weight (w)", fontsize=11, fontweight="bold")
    ax.set_ylim(-0.02, 1.02)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    sns.despine(ax=ax)

    if _style_ctx is not None:
        _style_ctx.__exit__(None, None, None)
    _save_and_maybe_show(fig, save_path=save_path, dpi=dpi, show=show, created=_created)
    return fig, ax
