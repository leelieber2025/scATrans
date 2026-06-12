"""
Bias correction (Huber regression on gene length + intron number).

The actual implementation lives in _utils._fit_huber_bias_correction so it can be
used from both the main analysis path and from permutation tasks without duplication.

Enhanced return: (residual, bias_info_dict) with fit diagnostics for transparency.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ._utils import _fit_huber_bias_correction as _raw_fit


def fit_huber_bias_correction(*args, **kwargs) -> tuple[np.ndarray, dict[str, Any]]:
    """Public/internal wrapper that returns (residual, bias_info)."""
    return _raw_fit(*args, **kwargs)


__all__ = ["fit_huber_bias_correction"]
