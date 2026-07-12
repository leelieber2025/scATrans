"""scatrans.tl.filter — internal package module."""

from __future__ import annotations

import logging
import math
import warnings
from typing import Any

import pandas as pd

from .._utils import (
    LEGACY_VELOCITY_RESIDUAL_COL,
    UNSPLICED_EXCESS_FDR_COL,
    UNSPLICED_EXCESS_RESIDUAL_COL,
)
from ._common import (
    _NOT_PROVIDED,
    HEURISTIC_FILTER_DEFAULTS,
    PSEUDOBULK_FILTER_DEFAULTS,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _read_filter_context(results: pd.DataFrame) -> dict[str, Any]:
    ctx = results.attrs.get("scatrans_filter_context")
    return dict(ctx) if isinstance(ctx, dict) else {}


def _builtin_significant_mask(
    var_df: pd.DataFrame,
    *,
    use_permutation: bool,
    extra_metadata: dict[str, Any],
) -> pd.Series:
    """Replicate the built-in ``significant`` mask from :func:`active_score`.

    Uses :data:`PSEUDOBULK_FILTER_DEFAULTS` when ``extra_metadata['is_pseudobulk']``
    is true (residual / active_score scales differ after sample aggregation);
    otherwise :data:`HEURISTIC_FILTER_DEFAULTS`. Explicit cutoffs in
    ``extra_metadata`` (``pval_cutoff``, ``logfc_cutoff``,
    ``unspliced_excess_fdr_cutoff``) still override the scale defaults.
    """
    index = var_df.index
    if not use_permutation or UNSPLICED_EXCESS_FDR_COL not in var_df.columns:
        return pd.Series(False, index=index)

    residual_col = (
        UNSPLICED_EXCESS_RESIDUAL_COL
        if UNSPLICED_EXCESS_RESIDUAL_COL in var_df.columns
        else LEGACY_VELOCITY_RESIDUAL_COL
    )
    scale_defaults = (
        PSEUDOBULK_FILTER_DEFAULTS
        if extra_metadata.get("is_pseudobulk")
        else HEURISTIC_FILTER_DEFAULTS
    )
    ue_fdr_cutoff = extra_metadata.get(
        "unspliced_excess_fdr_cutoff",
        scale_defaults["unspliced_excess_fdr_cutoff"],
    )
    sig_pval = extra_metadata.get("pval_cutoff", scale_defaults["pval_cutoff"])
    sig_logfc = extra_metadata.get("logfc_cutoff", scale_defaults["logfc_cutoff"])
    sig_resid = scale_defaults["unspliced_excess_residual_cutoff"]
    sig_active = scale_defaults["active_score_cutoff"]
    sig_active_fdr = scale_defaults["active_score_fdr_cutoff"]

    mask = pd.Series(True, index=index)
    if "p_adj" in var_df.columns:
        mask &= var_df["p_adj"] < sig_pval
    if "logFC" in var_df.columns:
        mask &= var_df["logFC"] > sig_logfc
    # MixedLM: p_adj tests mixedlm_coef, not sample-aware logFC. Require positive
    # coefficient so significance and effect direction agree (excludes sign-discordant genes).
    if "mixedlm_coef" in var_df.columns:
        coef = pd.to_numeric(var_df["mixedlm_coef"], errors="coerce")
        mask &= coef.notna() & (coef > 0)
    if residual_col in var_df.columns:
        mask &= var_df[residual_col] > sig_resid
    if "active_score" in var_df.columns:
        mask &= var_df["active_score"] >= sig_active
    if "valid_expr" in var_df.columns:
        mask &= var_df["valid_expr"]

    if extra_metadata.get("use_fdr_for_significance", True):
        mask &= var_df[UNSPLICED_EXCESS_FDR_COL] < ue_fdr_cutoff
        if "active_score_fdr" in var_df.columns and sig_active_fdr is not None:
            mask &= var_df["active_score_fdr"].notna() & (
                var_df["active_score_fdr"] < sig_active_fdr
            )

    if extra_metadata.get("use_delta_variance_pval") and "delta_var_pval" in var_df.columns:
        mask &= var_df["delta_var_pval"] < extra_metadata.get("delta_var_pval_cutoff", 0.05)

    return mask


def filter_active_genes(
    results: pd.DataFrame,
    *,
    preset: str | None = None,
    active_score_cutoff: Any = _NOT_PROVIDED,
    pval_cutoff: Any = _NOT_PROVIDED,
    padj_cutoff: Any = _NOT_PROVIDED,
    unspliced_excess_residual_cutoff: Any = _NOT_PROVIDED,
    logfc_cutoff: Any = _NOT_PROVIDED,
    logfc_direction: str = "up",
    active_score_fdr_cutoff: Any = _NOT_PROVIDED,
    unspliced_excess_fdr_cutoff: Any = _NOT_PROVIDED,
    effective_gamma_min: Any = _NOT_PROVIDED,
    effective_gamma_max: Any = _NOT_PROVIDED,
    delta_variance_min: Any = _NOT_PROVIDED,
    return_mask: bool = False,
    inplace: bool = False,
    **deprecated_kwargs: Any,
) -> pd.DataFrame | pd.Series:
    """Apply custom post-filtering to a results DataFrame (from `active_score` or `differential_expression`).

    This helper works for both:
    - Full `active_score` output (has `active_score` + unspliced_excess_residual).
    - Pure DE results from `differential_expression` (only logFC / p_adj + optional memento columns).

    It standardizes the common workflow:
    1. Run `active_score(...)` or `differential_expression(...)`.
    2. Use this function on the returned table to derive a final gene list.

    The function supports `preset` to automatically select reasonable default thresholds
    for different analysis modes:

    - preset="heuristic": stricter defaults suitable for typical single-cell data with default weights
      (active_score >= 55, unspliced_excess_residual > 1.0, etc.).
    - preset="pseudobulk": more lenient defaults that account for the much smaller
      magnitude of unspliced_excess_residual and active_score after sample-level aggregation
      (active_score >= 5, unspliced_excess_residual > 0.05, logFC > 0.2, etc.).
    - preset="significant" (aliases: "builtin", "active_score_significant"): exactly
      reproduces the built-in ``significant`` list from :func:`active_score` using metadata
      stored in ``all_results.attrs['scatrans_filter_context']``. Requires
      ``use_permutation=True`` on the upstream run.
    - preset=None (or "permissive"/"none"): apply only explicitly provided cutoffs; this is the most
      permissive / backward-compatible mode and returns nearly the full table (subject to any
      user-supplied thresholds).

    When ``all_results`` carries ``scatrans_filter_context`` with
    ``use_fdr_for_significance=False`` (common for pseudobulk with few samples),
    preset-based FDR cutoffs are skipped automatically unless you pass explicit
    ``unspliced_excess_fdr_cutoff`` / ``active_score_fdr_cutoff``.

    Presets are oriented toward target-group "activated" / upregulated signals (positive
    logFC + positive unspliced_excess_residual; direction defaults to "up").
    For downregulated or two-sided selection, pass ``logfc_direction="down"`` or
    ``"both"`` (with your desired logfc_cutoff). Residual magnitude cutoffs then
    follow the same direction (positive / negative / |residual|). Note that
    ``unspliced_excess_fdr`` is one-sided for positive residual and is skipped
    automatically when ``logfc_direction`` is not ``"up"``.

    If you explicitly pass any cutoff parameter, it takes precedence over the preset.

    Calling with no arguments (or only the DataFrame) and no preset returns the full
    `all_results` (fully permissive).

    Only filters corresponding to columns present in the DataFrame are applied.
    This is safe whether or not `use_permutation=True` or `use_mixed_model=True` was used.

    New options for power users:
    - return_mask=True: return the boolean mask (pd.Series) instead of the filtered DataFrame.
      Useful to combine with other logic or apply yourself.
    - inplace=True: mutate the input `results` DataFrame in-place (keeps only passing rows
      and re-sorts). Returns the (mutated) DataFrame. Use with caution; the input reference
      will be modified. Ignored when return_mask=True.

    Parameters
    ----------
    results : pd.DataFrame
        The `all_results` table returned as the third element of `active_score`.
    preset : str or None
        One of "heuristic", "pseudobulk", "significant", "permissive", "none".
        When provided, supplies recommended cutoff values for that analysis style
        for any parameters you did not explicitly pass.
    active_score_cutoff : float
        Minimum composite active transcription score (0-100).
    pval_cutoff : float
        Legacy name for the **adjusted** p-value cutoff (filters ``p_adj`` when present,
        else nominal ``p_val``). Prefer ``padj_cutoff`` for clarity.
    padj_cutoff : float
        Preferred name for the adjusted p-value cutoff. When both ``padj_cutoff`` and
        ``pval_cutoff`` are provided, ``padj_cutoff`` wins.
    unspliced_excess_residual_cutoff : float
        Magnitude threshold for bias-corrected unspliced (nascent) excess residual.
        Interpreted with ``logfc_direction``: ``up`` → residual > cutoff; ``down`` →
        residual < -cutoff; ``both`` → residual sign concordant with logFC
        (positive residual when logFC>0, negative when logFC<0). A cutoff of
        ``-inf`` disables the residual magnitude filter (permissive default).
    logfc_cutoff : float
        Magnitude threshold for logFC (treated as non-negative). See logfc_direction.
    logfc_direction : {"up", "down", "both"}
        Direction filter applied when a "logFC" column is present:
        - "up" (default): keep if logFC > logfc_cutoff (upregulated in target)
        - "down": keep if logFC < -logfc_cutoff (downregulated in target)
        - "both": keep if |logFC| > logfc_cutoff (differentially expressed either way)
        Residual magnitude filters follow the same direction (see above).
        Presets remain "active"/up-biased by default.
    active_score_fdr_cutoff : float or None
        If the column exists, max permutation FDR on the composite active_score (ranking aid).
    unspliced_excess_fdr_cutoff : float or None
        If the column exists (use_permutation=True), max permutation FDR on
        ``unspliced_excess_residual`` (recommended for final gene lists when
        ``logfc_direction="up"``). Skipped for ``"down"`` / ``"both"`` because the
        permutation residual test is one-sided for positive excess.
    effective_gamma_min : float
        Minimum reference-group effective gamma. See README section on effective_gamma.
    effective_gamma_max : float or None
        Optional upper bound on effective_gamma.
    delta_variance_min : float or None
        If the column exists (use_mixed_model=True), minimum variance fraction
        explained by condition.

    Returns
    -------
    pd.DataFrame
        Subset of the input. Sorted by active_score (desc) when present;
        otherwise by p_adj then logFC (direction-aware for pure DE tables).
    """
    if not isinstance(results, pd.DataFrame):
        raise ValueError("results must be the all_results DataFrame returned by active_score")

    velocity_residual_cutoff = _NOT_PROVIDED
    if "velocity_residual_cutoff" in deprecated_kwargs:
        velocity_residual_cutoff = deprecated_kwargs.pop("velocity_residual_cutoff")
        warnings.warn(
            "velocity_residual_cutoff is deprecated; use unspliced_excess_residual_cutoff.",
            DeprecationWarning,
            stacklevel=2,
        )
    if deprecated_kwargs:
        raise TypeError(
            "filter_active_genes() got unexpected keyword argument(s): "
            f"{', '.join(sorted(deprecated_kwargs))}"
        )

    df = results.copy()

    # Resolve values from preset + explicit overrides
    if preset is not None:
        p = preset.lower()
        if p in ("significant", "builtin", "active_score_significant"):
            ctx = _read_filter_context(df)
            if not ctx:
                raise ValueError(
                    "preset='significant' requires all_results from active_score with "
                    "attrs['scatrans_filter_context'] (re-run active_score on a recent version)."
                )
            mask = _builtin_significant_mask(
                df,
                use_permutation=bool(ctx.get("use_permutation")),
                extra_metadata=ctx,
            )
            if return_mask:
                return mask
            filtered = df.loc[mask].copy()
            if "active_score" in filtered.columns:
                filtered = filtered.sort_values("active_score", ascending=False)
            elif "p_adj" in filtered.columns:
                filtered = filtered.sort_values("p_adj", ascending=True)
            if inplace:
                # Drop non-passing rows; re-sort surviving rows (preserves identity + attrs)
                to_drop = results.index.difference(filtered.index)
                if len(to_drop) > 0:
                    results.drop(index=to_drop, inplace=True)
                if "active_score" in results.columns:
                    results.sort_values("active_score", ascending=False, inplace=True)
                elif "p_adj" in results.columns:
                    results.sort_values("p_adj", ascending=True, inplace=True)
                return results
            return filtered
        if p in ("heuristic", "single_cell", "default"):
            preset_vals = {
                **HEURISTIC_FILTER_DEFAULTS,
                "velocity_residual_cutoff": HEURISTIC_FILTER_DEFAULTS[
                    "unspliced_excess_residual_cutoff"
                ],
            }
        elif p in ("pseudobulk", "bulk"):
            preset_vals = {
                **PSEUDOBULK_FILTER_DEFAULTS,
                "velocity_residual_cutoff": PSEUDOBULK_FILTER_DEFAULTS[
                    "unspliced_excess_residual_cutoff"
                ],
            }
        elif p in ("permissive", "none", "all", "no_filter"):
            preset_vals = {
                "active_score_cutoff": 0.0,
                "pval_cutoff": float("inf"),
                "velocity_residual_cutoff": float("-inf"),
                "unspliced_excess_residual_cutoff": float("-inf"),
                "logfc_cutoff": float("inf"),
                "active_score_fdr_cutoff": float("inf"),
                "unspliced_excess_fdr_cutoff": float("inf"),
                "effective_gamma_min": float("-inf"),
                "effective_gamma_max": None,
                "delta_variance_min": None,
            }
        else:
            raise ValueError(
                f"Unknown preset '{preset}'. "
                "Valid presets: 'heuristic', 'pseudobulk', 'significant' "
                "(aliases: 'builtin', 'active_score_significant'), 'permissive'."
            )
    else:
        preset_vals = {}

    user_set_ue_fdr = unspliced_excess_fdr_cutoff is not _NOT_PROVIDED
    user_set_as_fdr = active_score_fdr_cutoff is not _NOT_PROVIDED

    # Apply preset only where user did not explicitly provide a value
    def _resolve(name: str, current: Any, default: Any) -> Any:
        if current is not _NOT_PROVIDED:
            return current
        return preset_vals.get(name, default)

    def _coerce_numeric_cutoff(val: Any, default: float, name: str) -> float:
        """Coerce a required numeric cutoff; reject None / non-numeric with ValueError."""
        if val is _NOT_PROVIDED:
            return float(default)
        if val is None:
            raise ValueError(f"{name} must be numeric, got None")
        try:
            return float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be numeric, got {val!r}") from exc

    def _coerce_optional_numeric_cutoff(val: Any, name: str) -> float | None:
        """Coerce optional cutoff: None stays None (skip filter); else require numeric."""
        if val is None:
            return None
        try:
            return float(val)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} must be numeric or None, got {val!r}") from exc

    active_score_cutoff = _coerce_numeric_cutoff(
        _resolve("active_score_cutoff", active_score_cutoff, 0.0),
        0.0,
        "active_score_cutoff",
    )
    # padj_cutoff is the preferred name (filters p_adj); pval_cutoff is the legacy alias.
    user_set_padj = padj_cutoff is not _NOT_PROVIDED
    user_set_pval = pval_cutoff is not _NOT_PROVIDED
    if user_set_padj and user_set_pval:
        logger.warning(
            "Both padj_cutoff and pval_cutoff were provided to filter_active_genes; "
            "using padj_cutoff (pval_cutoff is a legacy alias for adjusted p)."
        )
        pval_cutoff = padj_cutoff
    elif user_set_padj:
        pval_cutoff = padj_cutoff
    pval_cutoff = _coerce_numeric_cutoff(
        _resolve("pval_cutoff", pval_cutoff, float("inf")), float("inf"), "pval_cutoff"
    )
    if pval_cutoff < 0 or (not math.isfinite(pval_cutoff) and not math.isinf(pval_cutoff)):
        raise ValueError(
            "padj_cutoff / pval_cutoff must be non-negative, finite, or +inf (permissive)."
        )
    if (
        velocity_residual_cutoff is not _NOT_PROVIDED
        and unspliced_excess_residual_cutoff is _NOT_PROVIDED
    ):
        unspliced_excess_residual_cutoff = velocity_residual_cutoff
    velocity_residual_cutoff = _coerce_numeric_cutoff(
        _resolve("velocity_residual_cutoff", velocity_residual_cutoff, float("-inf")),
        float("-inf"),
        "velocity_residual_cutoff",
    )
    unspliced_excess_residual_cutoff = _coerce_numeric_cutoff(
        _resolve(
            "unspliced_excess_residual_cutoff",
            unspliced_excess_residual_cutoff,
            velocity_residual_cutoff,
        ),
        float(velocity_residual_cutoff),
        "unspliced_excess_residual_cutoff",
    )
    logfc_cutoff = _coerce_numeric_cutoff(
        _resolve("logfc_cutoff", logfc_cutoff, float("inf")), float("inf"), "logfc_cutoff"
    )
    # logfc_direction is not preset-driven (presets remain up-biased); normalize here
    dir_raw = str(logfc_direction).lower() if logfc_direction is not None else "up"
    if dir_raw in {"up", "positive", "pos", "u"}:
        direction = "up"
    elif dir_raw in {"down", "negative", "neg", "d"}:
        direction = "down"
    elif dir_raw in {"both", "two_sided", "twosided", "abs", "absolute", "either", "any", "b"}:
        direction = "both"
    else:
        raise ValueError(
            f'logfc_direction={logfc_direction!r} not recognized. Use one of: "up", "down", "both".'
        )
    # FDR cutoffs: None means skip that filter (documented); otherwise numeric / ±inf
    active_score_fdr_cutoff = _coerce_optional_numeric_cutoff(
        _resolve("active_score_fdr_cutoff", active_score_fdr_cutoff, float("inf")),
        "active_score_fdr_cutoff",
    )
    unspliced_excess_fdr_cutoff = _coerce_optional_numeric_cutoff(
        _resolve("unspliced_excess_fdr_cutoff", unspliced_excess_fdr_cutoff, float("inf")),
        "unspliced_excess_fdr_cutoff",
    )
    # Optional gamma bounds: None means "no filter" (presets default to None)
    effective_gamma_min = _coerce_optional_numeric_cutoff(
        _resolve("effective_gamma_min", effective_gamma_min, None),
        "effective_gamma_min",
    )
    effective_gamma_max = _coerce_optional_numeric_cutoff(
        _resolve("effective_gamma_max", effective_gamma_max, None),
        "effective_gamma_max",
    )
    delta_variance_min = _coerce_optional_numeric_cutoff(
        _resolve("delta_variance_min", delta_variance_min, None),
        "delta_variance_min",
    )

    ctx = _read_filter_context(df)
    if ctx and not ctx.get("use_fdr_for_significance", True):
        if not user_set_ue_fdr:
            unspliced_excess_fdr_cutoff = float("inf")
        if not user_set_as_fdr:
            active_score_fdr_cutoff = float("inf")
        logger.info(
            "Permutation FDR was disabled for the built-in significant list (%s); "
            "skipping preset FDR cutoffs in filter_active_genes. "
            "Use preset='significant' to match the built-in list exactly, or pass explicit "
            "FDR cutoffs to override.",
            ctx.get("perm_disabled_reason", "small_permutation_space"),
        )

    mask = pd.Series(True, index=df.index)

    # Core filters
    # Prefer adjusted p-value when present (consistent with active_score internal significant mask
    # and common user expectation). Fall back to nominal p_val only if p_adj is absent.
    if "active_score" in df.columns:
        active_vals = pd.to_numeric(df["active_score"], errors="coerce")
        mask &= active_vals.notna() & (active_vals >= active_score_cutoff)
    if "p_adj" in df.columns:
        padj_vals = pd.to_numeric(df["p_adj"], errors="coerce")
        mask &= padj_vals.notna() & (padj_vals < pval_cutoff)
    elif "p_val" in df.columns:
        pval_vals = pd.to_numeric(df["p_val"], errors="coerce")
        mask &= pval_vals.notna() & (pval_vals < pval_cutoff)
    residual_col = (
        UNSPLICED_EXCESS_RESIDUAL_COL
        if UNSPLICED_EXCESS_RESIDUAL_COL in df.columns
        else LEGACY_VELOCITY_RESIDUAL_COL
    )
    if residual_col in df.columns:
        resid_vals = pd.to_numeric(df[residual_col], errors="coerce")
        rc = float(unspliced_excess_residual_cutoff)
        # -inf disables residual magnitude filtering (permissive default).
        if math.isinf(rc) and rc < 0:
            mask &= resid_vals.notna() & (resid_vals > rc)
        else:
            mag = abs(rc) if math.isfinite(rc) else rc
            if direction == "up":
                mask &= resid_vals.notna() & (resid_vals > mag)
            elif direction == "down":
                mask &= resid_vals.notna() & (resid_vals < -mag)
            else:
                # both: require residual sign concordant with logFC when available
                # (up+positive residual OR down+negative residual); else |residual|.
                if "logFC" in df.columns:
                    logfc_for_res = pd.to_numeric(df["logFC"], errors="coerce")
                    concordant = ((logfc_for_res > 0) & (resid_vals > mag)) | (
                        (logfc_for_res < 0) & (resid_vals < -mag)
                    )
                    mask &= resid_vals.notna() & logfc_for_res.notna() & concordant
                else:
                    mask &= resid_vals.notna() & (resid_vals.abs() > mag)
    if "logFC" in df.columns:
        lc = logfc_cutoff
        if isinstance(lc, (int, float)) and math.isinf(lc):
            # permissive preset: no logFC threshold
            pass
        else:
            if not math.isfinite(lc):
                raise ValueError("logfc_cutoff must be finite or +inf (permissive).")
            if lc < 0:
                lc = -lc  # always use positive magnitude
            logfc_vals = pd.to_numeric(df["logFC"], errors="coerce")
            if direction == "up":
                mask &= logfc_vals.notna() & (logfc_vals > lc)
            elif direction == "down":
                mask &= logfc_vals.notna() & (logfc_vals < -lc)
            else:  # both
                mask &= logfc_vals.notna() & (logfc_vals.abs() > lc)
    # MixedLM: keep effect direction consistent with what p_adj tests (mixedlm_coef).
    if "mixedlm_coef" in df.columns and direction in ("up", "down"):
        coef_vals = pd.to_numeric(df["mixedlm_coef"], errors="coerce")
        if direction == "up":
            mask &= coef_vals.notna() & (coef_vals > 0)
        else:
            mask &= coef_vals.notna() & (coef_vals < 0)

    # Permutation FDR on composite score (optional ranking filter)
    if active_score_fdr_cutoff is not None and "active_score_fdr" in df.columns:
        fdr = pd.to_numeric(df["active_score_fdr"], errors="coerce")
        if math.isinf(active_score_fdr_cutoff):
            # Permissive: NaN FDR = not computed / no permutation — do not drop.
            mask &= fdr.isna() | (fdr < active_score_fdr_cutoff)
        else:
            mask &= fdr.notna() & (fdr < active_score_fdr_cutoff)

    # Permutation FDR on unspliced excess residual (recommended significance filter).
    # One-sided for positive residual only — skip when user selected down/both.
    if unspliced_excess_fdr_cutoff is not None and UNSPLICED_EXCESS_FDR_COL in df.columns:
        apply_ue_fdr = True
        if direction != "up" and not (
            isinstance(unspliced_excess_fdr_cutoff, float)
            and math.isinf(unspliced_excess_fdr_cutoff)
        ):
            logger.warning(
                "logfc_direction=%r: unspliced_excess_fdr is a one-sided test for "
                "positive residual; skipping residual FDR filter. Use residual "
                "magnitude cutoffs (direction-aware) for down/both selection.",
                direction,
            )
            apply_ue_fdr = False
        if apply_ue_fdr:
            fdr = pd.to_numeric(df[UNSPLICED_EXCESS_FDR_COL], errors="coerce")
            if math.isinf(unspliced_excess_fdr_cutoff):
                mask &= fdr.isna() | (fdr < unspliced_excess_fdr_cutoff)
            else:
                mask &= fdr.notna() & (fdr < unspliced_excess_fdr_cutoff)

    # effective_gamma bounds (only when user/preset sets a finite min/max)
    if "effective_gamma" in df.columns and (
        effective_gamma_min is not None or effective_gamma_max is not None
    ):
        gamma = pd.to_numeric(df["effective_gamma"], errors="coerce")
        if effective_gamma_min is not None:
            mask &= gamma.notna() & (gamma > effective_gamma_min)
        if effective_gamma_max is not None:
            mask &= gamma.notna() & (gamma < effective_gamma_max)

    # Delta variance
    if delta_variance_min is not None and "delta_variance" in df.columns:
        mask &= df["delta_variance"] >= delta_variance_min

    filtered = df[mask].copy()

    # Sorting: prefer active_score when present (velocity + DE composite),
    # otherwise fall back to p_adj, then logFC according to direction
    # (strongest up first, strongest down first, or largest |logFC| for both).
    if "active_score" in filtered.columns:
        filtered = filtered.sort_values("active_score", ascending=False)
    elif "p_adj" in filtered.columns:
        sort_cols = ["p_adj"]
        ascending = [True]
        if "logFC" in filtered.columns:
            if direction == "down":
                sort_cols.append("logFC")
                ascending.append(True)  # most negative first
            elif direction == "both":
                filtered = filtered.assign(_logfc_abs=filtered["logFC"].abs())
                sort_cols.append("_logfc_abs")
                ascending.append(False)
            else:
                sort_cols.append("logFC")
                ascending.append(False)
        filtered = filtered.sort_values(sort_cols, ascending=ascending)
        if "_logfc_abs" in filtered.columns:
            filtered = filtered.drop(columns=["_logfc_abs"])
    else:
        filtered = filtered.sort_values(filtered.columns[0])

    if return_mask:
        return mask

    if inplace:
        # Mutate caller's DataFrame: drop non-matching rows, keep sorted order
        to_drop = results.index.difference(filtered.index)
        if len(to_drop) > 0:
            results.drop(index=to_drop, inplace=True)
        # Re-apply the sort in place on the surviving rows
        if "active_score" in results.columns:
            results.sort_values("active_score", ascending=False, inplace=True)
        elif "p_adj" in results.columns:
            sort_cols = ["p_adj"]
            ascending = [True]
            if "logFC" in results.columns:
                if direction == "down":
                    sort_cols.append("logFC")
                    ascending.append(True)
                elif direction == "both":
                    # Use direct assignment so inplace=True actually mutates the caller's
                    # DataFrame (assign() would rebind the local variable and lose the effect).
                    results["_logfc_abs"] = results["logFC"].abs()
                    sort_cols.append("_logfc_abs")
                    ascending.append(False)
                else:
                    sort_cols.append("logFC")
                    ascending.append(False)
            results.sort_values(sort_cols, ascending=ascending, inplace=True)
            if "_logfc_abs" in results.columns:
                results.drop(columns=["_logfc_abs"], inplace=True)
        return results

    return filtered
