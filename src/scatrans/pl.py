"""Plotting functions for scATrans.

NOTE: In the current release, only a basic diagnostic scatter plot is
generated inside `active_score(show_plot=True)`. 

Full publication-ready functions (comet_plot, volcano_3d, bias_diagnostic_plot,
enrich_dotplot, etc.) are planned for a future release.

This module currently provides:
- set_style() : consistent matplotlib style
- Placeholder functions that raise clear NotImplementedError with guidance.
"""

import matplotlib.pyplot as plt
import warnings
from typing import Optional, Any


def set_style():
    """Apply a clean seaborn-like style for scATrans figures."""
    try:
        plt.style.use('seaborn-v0_8-whitegrid')
    except Exception:
        plt.style.use('seaborn-whitegrid')


def _not_implemented_plot(name: str):
    def wrapper(*args, **kwargs):
        raise NotImplementedError(
            f"scat.pl.{name}() is not yet implemented in this version of scATrans.\n\n"
            "Current options:\n"
            "  1. Use the built-in diagnostic plot returned by active_score(show_plot=True)\n"
            "  2. Use the returned `all_results` DataFrame to create custom plots with seaborn/matplotlib/plotly\n"
            "  3. For comet/volcano-style plots, see examples in the scATrans paper or implement using\n"
            "     scatter + annotation of top genes.\n\n"
            "We welcome contributions! See GitHub issues."
        )
    return wrapper


# Placeholder functions matching README expectations
comet_plot = _not_implemented_plot("comet_plot")
volcano_plot = _not_implemented_plot("volcano_plot")
volcano_3d = _not_implemented_plot("volcano_3d")
bias_diagnostic_plot = _not_implemented_plot("bias_diagnostic_plot")
enrich_dotplot = _not_implemented_plot("enrich_dotplot")
enrich_barplot = _not_implemented_plot("enrich_barplot")


# You can uncomment and implement real versions later
# def comet_plot(all_results, top_n=12, save_path=None, **kwargs):
#     ...
