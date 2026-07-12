"""scatrans.tl.design — internal package module."""

from __future__ import annotations

import logging
from typing import Any

import pandas as pd

from .. import qc as _qc
from .._utils import (
    _normalize_group_label,
    comb,
)
from ._common import (
    _PERM_FDR_MIN_SUCCESS,
    MIXED_MODEL_MIN_SAMPLES_PER_GROUP,
    _resolve_velocity_layer_keys,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


WORKFLOW_PRESETS: dict[str, dict[str, Any]] = {
    "explore": {
        "label": "Quick exploration (ranking only, no permutation)",
        "active_score_kwargs": {
            "use_permutation": False,
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "heuristic",
    },
    "report": {
        "label": "Manuscript reporting with permutation FDR on unspliced excess",
        "active_score_kwargs": {
            "use_permutation": True,
            "n_perm": 500,
            "perm_de_backend": "same",
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "heuristic",
    },
    "pseudobulk_report": {
        "label": "Multi-replicate pseudobulk DE + permutation FDR",
        "active_score_kwargs": {
            "use_pseudobulk": True,
            "pseudobulk_de_backend": "pydeseq2",
            "use_permutation": True,
            "n_perm": 200,
            "perm_de_backend": "same",
            "mode": "heuristic",
            "ranking_mode": "composite",
        },
        "filter_preset": "pseudobulk",
    },
    "nascent_focus": {
        "label": "Rank by bias-corrected nascent excess residual only",
        "active_score_kwargs": {
            "ranking_mode": "nascent_excess",
            "use_permutation": False,
            "mode": "heuristic",
        },
        "filter_preset": "heuristic",
    },
}


def _permutation_power_guidance(
    n_cells_target: int,
    n_cells_reference: int,
    n_genes: int,
    *,
    n_samples_target: int | None = None,
    n_samples_reference: int | None = None,
    n_perm: int = 100,
    n_jobs: int = -1,
) -> dict[str, Any]:
    """Rough runtime / permutation-space guidance for diagnose_design."""
    import os

    is_pseudobulk = (
        n_samples_target is not None
        and n_samples_reference is not None
        and min(n_samples_target, n_samples_reference) >= 2
    )
    max_exact: int | None = None
    if is_pseudobulk:
        n_t = max(1, int(n_samples_target or 1))
        n_r = max(1, int(n_samples_reference or 1))
        if n_t + n_r <= 30:
            max_exact = max(1, comb(n_t + n_r, n_t) - 1)

    cpu = os.cpu_count() or 4
    effective_cores = cpu if n_jobs == -1 else max(1, min(n_jobs, cpu))
    sec_per_perm = max(0.03, (n_genes / 5000.0) * 0.10)
    est_minutes = n_perm * sec_per_perm / effective_cores / 60.0

    notes: list[str] = [
        "Estimates assume heuristic mode on a typical laptop/workstation; "
        "advanced mode, Memento, or perm_de_backend='same' with PyDESeq2 can be slower.",
    ]
    if is_pseudobulk and max_exact is not None and max_exact < n_perm:
        notes.append(
            f"Pseudobulk design allows at most {max_exact} exact label permutations; "
            f"auto_adjust_n_perm=True will cap n_perm accordingly."
        )
    if min(n_cells_target, n_cells_reference) < 50:
        notes.append(
            "Very small cell counts: permutation FDR on unspliced excess will have low power "
            "even if runtime is acceptable."
        )

    return {
        "n_genes": int(n_genes),
        "default_n_perm": int(n_perm),
        "is_pseudobulk_context": is_pseudobulk,
        "max_exact_permutations_pseudobulk": max_exact,
        "estimated_runtime_minutes_heuristic": round(est_minutes, 1),
        "effective_cores_assumed": effective_cores,
        "notes": notes,
    }


def diagnose_design(
    adata_input: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None = None,
    _min_cells_per_sample: int = 10,  # reserved/internal; not yet used (will affect sample filtering/power in future)
    *,
    copy_input: bool = True,
) -> dict[str, Any]:
    """
    Analyze the experimental design and provide guidance on suitable analysis choices
    and expected power/limitations.

    This is intended as a pre-flight or post-subset diagnostic to help users
    interpret warnings and choose between single-cell, pseudobulk, or mixed-model paths.

    Returns a dictionary with keys:
      - n_cells_target, n_cells_reference
      - n_samples_target, n_samples_reference (if sample_col provided)
      - unspliced_global_fraction
      - recommendations: list of human-readable strings
      - warnings: list of human-readable strings
      - suggested_preset: filter_active_genes preset ("heuristic" or "pseudobulk")
      - power_summary: permutation runtime / power guidance dict
      - workflow_preset: recommended entry from WORKFLOW_PRESETS

    Note: ``_min_cells_per_sample`` is reserved for future use and currently ignored.

    copy_input : bool, default True
        If True (default), a full deep copy of the input AnnData is made before
        reading. Set False for a zero-copy read-only diagnostic when the caller
        guarantees the input will not be mutated (saves large amounts of memory
        and time on big datasets with many layers).
    """
    adata = adata_input.copy() if copy_input else adata_input

    if groupby not in adata.obs.columns:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs")

    norm_groups = adata.obs[groupby].map(_normalize_group_label)
    target_mask = norm_groups == _normalize_group_label(target_group)
    ref_mask = norm_groups == _normalize_group_label(reference_group)

    n_t = int(target_mask.sum())
    n_r = int(ref_mask.sum())

    result: dict[str, Any] = {
        "n_cells_target": n_t,
        "n_cells_reference": n_r,
        "n_samples_target": None,
        "n_samples_reference": None,
        "unspliced_global_fraction": None,
        "recommendations": [],
        "warnings": [],
        "suggested_preset": None,
        "workflow_preset": "explore",
        "power_summary": None,
    }

    # Global unspliced fraction (important technical QC)
    layer_keys = _resolve_velocity_layer_keys(adata)
    if layer_keys is not None:
        spliced_key, unspliced_key = layer_keys
        try:
            ufrac = _qc.unspliced_global(
                adata,
                spliced_key=spliced_key,
                unspliced_key=unspliced_key,
                warn_threshold=0.5,
            )
            result["unspliced_global_fraction"] = float(ufrac)
            if ufrac > 0.5:
                result["warnings"].append(
                    f"Global unspliced fraction is high ({ufrac:.1%}). "
                    "This often indicates nuclear enrichment or gDNA contamination and can "
                    "reduce the reliability of velocity-based signals."
                )
        except Exception as e:
            logger.debug("Could not compute global unspliced fraction: %s", e)

    # Sample structure
    if sample_col:
        if sample_col not in adata.obs.columns:
            result["warnings"].append(
                f"sample_col={sample_col!r} was provided but is not in adata.obs columns. "
                f"Available (first 20): {list(map(str, adata.obs.columns[:20]))}. "
                "Design guidance will treat cells as independent until a valid sample_col is set."
            )
            result["recommendations"].append(
                f"Fix sample_col (got {sample_col!r}) so replicate structure can be detected."
            )
        else:
            n_s_t = adata.obs.loc[target_mask, sample_col].nunique()
            n_s_r = adata.obs.loc[ref_mask, sample_col].nunique()
            result["n_samples_target"] = int(n_s_t)
            result["n_samples_reference"] = int(n_s_r)

            if min(n_s_t, n_s_r) < 3:
                result["warnings"].append(
                    f"Very few biological samples per group (target={n_s_t}, reference={n_s_r}). "
                    "Pseudobulk aggregation will have extremely low power for velocity delta. "
                    "Prefer use_pseudobulk=True with PyDESeq2 for DE ranking; "
                    "keep use_permutation=False and interpret velocity residuals cautiously."
                )
            elif min(n_s_t, n_s_r) < MIXED_MODEL_MIN_SAMPLES_PER_GROUP:
                result["warnings"].append(
                    f"Small number of biological samples per group (target={n_s_t}, reference={n_s_r}). "
                    "Power for detecting differential nascent RNA excess will be limited. "
                    "Permutation-based FDR is unreliable with so few label shuffles — prefer "
                    "use_pseudobulk=True without permutation for ranking, then filter_active_genes "
                    "(preset='pseudobulk') or DE p_adj for significance."
                )
                result["suggested_preset"] = "pseudobulk"

            if min(n_s_t, n_s_r) >= MIXED_MODEL_MIN_SAMPLES_PER_GROUP:
                result["recommendations"].append(
                    "With sufficient biological replicates, pseudobulk (PyDESeq2) and "
                    "use_mixed_model=True are both viable. See the small-sample guidance "
                    "in the documentation."
                )
            else:
                result["recommendations"].append(
                    f"With only {min(n_s_t, n_s_r)} sample(s) per group, use use_pseudobulk=True "
                    f"(pseudobulk_de_backend='pydeseq2') rather than use_mixed_model=True "
                    f"(requires >={MIXED_MODEL_MIN_SAMPLES_PER_GROUP} samples per group)."
                )
    else:
        result["recommendations"].append(
            "No sample_col provided. The analysis will treat cells as independent. "
            "If cells come from multiple biological replicates, consider providing sample_col "
            "and using either use_pseudobulk=True or use_mixed_model=True to avoid "
            "pseudoreplication."
        )

    # Very small total cell numbers
    if min(n_t, n_r) < 50:
        result["warnings"].append(
            f"Very small number of cells in at least one group (target={n_t}, reference={n_r}). "
            "Velocity delta estimation and any downstream permutation testing will have low power."
        )

    if n_r < 80:
        result["recommendations"].append(
            f"Reference group has {n_r} cells — consider gamma_method='empirical_bayes' "
            "for more stable per-gene reference gamma (robust log-ratio shrinkage)."
        )

    # Suggest filter_active_genes preset
    if result["suggested_preset"] is None:
        n_samp_t = result.get("n_samples_target")
        if sample_col and n_samp_t is not None and int(n_samp_t) >= 5:
            result["suggested_preset"] = "pseudobulk"
        else:
            result["suggested_preset"] = "heuristic"

    # Workflow preset (ties to WORKFLOW_PRESETS for active_score kwargs).
    # "pseudobulk_report" enables use_pseudobulk (+ optional perm in the preset
    # kwargs). recommend_workflow() then disables use_permutation when the exact
    # shuffle space is too small (< _PERM_FDR_MIN_SUCCESS), so diagnose warnings
    # about weak permutation and the preset stay consistent at the public API.
    min_samples = min(result.get("n_samples_target") or 0, result.get("n_samples_reference") or 0)
    if sample_col and sample_col in adata.obs.columns and min_samples >= 3:
        result["workflow_preset"] = "pseudobulk_report"
    elif min(n_t, n_r) >= 100 and not result["warnings"]:
        result["workflow_preset"] = "report"
    else:
        result["workflow_preset"] = "explore"

    result["power_summary"] = _permutation_power_guidance(
        n_t,
        n_r,
        adata.n_vars,
        n_samples_target=result.get("n_samples_target"),
        n_samples_reference=result.get("n_samples_reference"),
    )

    # General advice
    result["recommendations"].append(
        "After running active_score, always inspect adata.uns['scatrans']['diagnostics'] "
        "and the distributions in the returned all_results DataFrame before applying cutoffs."
    )

    # Print a concise user-facing summary
    logger.info("Design diagnosis:")
    logger.info("  Cells — target: %d | reference: %d", n_t, n_r)
    if result["n_samples_target"] is not None:
        logger.info(
            "  Samples — target: %d | reference: %d",
            result["n_samples_target"],
            result["n_samples_reference"],
        )
    if result["unspliced_global_fraction"] is not None:
        logger.info(
            "  Global unspliced fraction: %.1f%%", result["unspliced_global_fraction"] * 100
        )

    for w in result["warnings"]:
        logger.warning("  [WARNING] %s", w)
    for r in result["recommendations"]:
        logger.info("  [RECOMMENDATION] %s", r)

    ps = result.get("power_summary") or {}
    if ps.get("estimated_runtime_minutes_heuristic") is not None:
        logger.info(
            "  Permutation runtime (heuristic, n_perm=%d): ~%.1f min on %d cores",
            ps.get("default_n_perm", 100),
            ps["estimated_runtime_minutes_heuristic"],
            ps.get("effective_cores_assumed", 1),
        )

    return result


def recommend_workflow(
    adata_input: Any,
    groupby: str,
    target_group: str,
    reference_group: str,
    sample_col: str | None = None,
) -> dict[str, Any]:
    """
    High-level recommendation for analysis path based on experimental design.

    This is a thin, user-friendly wrapper around diagnose_design that returns
    actionable preset + backend suggestions.

    Returns keys:
      - workflow_preset: key into ``WORKFLOW_PRESETS`` (e.g. ``"pseudobulk_report"``)
      - preset_config: full preset dict (label, active_score_kwargs, filter_preset)
      - recommended_preset: ``filter_active_genes`` preset name
      - de_backend: ``"scanpy"`` | ``"pydeseq2"`` | ...
      - suggested_kwargs: merged kwargs for :func:`active_score`
      - warnings, recommendations, power_summary
      - full_diagnosis: raw dict from :func:`diagnose_design`

    Example:
        rec = scat.recommend_workflow(adata, "condition", "GA", "Ctrl", sample_col="sample")
        adata, sig, res = scat.active_score(
            adata, groupby="condition", target_group="GA", reference_group="Ctrl",
            **rec["suggested_kwargs"],
        )
        candidates = scat.filter_active_genes(res, preset=rec["filter_preset"])
    """
    diag = diagnose_design(
        adata_input,
        groupby=groupby,
        target_group=target_group,
        reference_group=reference_group,
        sample_col=sample_col,
        copy_input=True,  # public entry point: be safe by default
    )

    workflow_key = diag.get("workflow_preset") or "explore"
    preset_config = WORKFLOW_PRESETS.get(workflow_key, WORKFLOW_PRESETS["explore"])
    filter_preset = preset_config.get("filter_preset", diag.get("suggested_preset") or "heuristic")

    suggested_kwargs: dict[str, Any] = dict(preset_config.get("active_score_kwargs", {}))
    if sample_col and sample_col in getattr(adata_input, "obs", pd.DataFrame()).columns:
        suggested_kwargs.setdefault("sample_col", sample_col)
    if suggested_kwargs.get("use_pseudobulk"):
        de_backend = "pydeseq2"
        suggested_kwargs.setdefault("pseudobulk_de_backend", "pydeseq2")
    else:
        de_backend = "scanpy"
        suggested_kwargs.setdefault("de_method", "wilcoxon")

    warnings = list(diag.get("warnings", []))
    recommendations = list(diag.get("recommendations", []))
    power = diag.get("power_summary") or {}
    max_exact = power.get("max_exact_permutations_pseudobulk")

    if suggested_kwargs.get("use_permutation") and max_exact is not None:
        if max_exact < _PERM_FDR_MIN_SUCCESS:
            suggested_kwargs["use_permutation"] = False
            filter_preset = diag.get("suggested_preset") or "pseudobulk"
            warnings.append(
                f"Only ~{max_exact} exact pseudobulk permutation shuffles are possible "
                f"(<{_PERM_FDR_MIN_SUCCESS}); suggested_kwargs sets use_permutation=False. "
                "Use filter_active_genes(preset='pseudobulk') on all_results for exploratory "
                "gene lists, or DE p_adj from differential_expression()."
            )
        else:
            filter_preset = "significant"

    rec = {
        "workflow_preset": workflow_key,
        "preset_config": preset_config,
        "recommended_preset": filter_preset,
        "filter_preset": filter_preset,
        "de_backend": de_backend,
        "use_permutation": bool(suggested_kwargs.get("use_permutation", False)),
        "suggested_kwargs": suggested_kwargs,
        "warnings": warnings,
        "recommendations": recommendations,
        "power_summary": power,
        "full_diagnosis": diag,
    }

    logger.info(
        "Workflow recommendation: workflow_preset=%s, filter_preset=%s, de_backend=%s",
        rec["workflow_preset"],
        rec["filter_preset"],
        rec["de_backend"],
    )
    return rec
