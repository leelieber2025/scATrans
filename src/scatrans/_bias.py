"""
Bias correction (Huber regression on gene length + intron number).

The actual implementation lives in _utils._fit_huber_bias_correction so it can be
used from both the main analysis path and from permutation tasks without duplication.
"""

from __future__ import annotations

from ._utils import _fit_huber_bias_correction as fit_huber_bias_correction

__all__ = ["fit_huber_bias_correction"]
