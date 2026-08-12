"""
Synthetic toy datasets for a first, zero-download run of scATrans.

These are NOT real biology: real gene symbols are used (sampled from the
bundled gene-feature tables, so ``add_gene_features``/bias correction works
without warnings) but the per-gene mechanism assignment and all counts are
randomly simulated. Use ``load_toy`` to check that installation and the
default :func:`~scatrans.tl.partition_de_by_mechanism` call work end to end
before touching real data — not to draw any biological conclusion.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from .pp_bias import _open_package_data

__all__ = ["load_toy"]

_DEFAULT_FEATURE_FILE = {
    "mouse": "mouse_2020A_gene_features.parquet",
    "human": "human_GRCh38_2024A_gene_features.parquet",
}


def _sample_real_gene_names(organism: str, n_genes: int, rng: np.random.Generator) -> list[str]:
    """Sample real, unique protein-coding gene symbols from a bundled feature table.

    Falls back to placeholder names (``GeneNNNN``) if the bundled table can't
    be read for any reason — the toy dataset still works, just with a lower
    ``add_gene_features`` mapping rate.
    """
    organism_norm = str(organism).lower()
    filename = _DEFAULT_FEATURE_FILE.get(organism_norm)
    if filename is None:
        raise ValueError(f"organism must be 'mouse' or 'human', got {organism!r}")

    try:
        with _open_package_data(filename) as path:
            table = pd.read_parquet(path, columns=["gene_name", "gene_type"])
        pool = table.loc[table["gene_type"] == "protein_coding", "gene_name"].dropna().unique()
        if len(pool) >= n_genes:
            return list(rng.choice(pool, size=n_genes, replace=False))
    except Exception:
        pass
    return [f"Gene{i:04d}" for i in range(n_genes)]


def load_toy(
    n_cells: int = 240,
    n_genes: int = 200,
    n_samples: int = 4,
    organism: str = "mouse",
    seed: int = 0,
) -> AnnData:
    """Small synthetic AnnData with spliced/unspliced layers, ready for
    :func:`~scatrans.tl.partition_de_by_mechanism`.

    Two conditions (``"Control"`` / ``"Disease"``, half the cells each) and
    four gene categories are simulated: housekeeping (no change),
    transcription-driven up (spliced *and* unspliced both increase, unspliced
    proportionally more — a transient synthesis-rate increase), stabilization-
    driven up (spliced increases, unspliced unchanged — slower decay, same
    synthesis rate), and a generic down-regulated group. Ground truth is kept
    in ``adata.var["ground_truth_mechanism"]`` for comparison only; real data
    never has this column.

    Parameters
    ----------
    n_cells
        Total cells, split evenly Control/Disease.
    n_genes
        Total genes (~60% housekeeping, ~15% transcription-up, ~15%
        stabilization-up, ~10% down).
    n_samples
        Total biological replicates (donors), split evenly across the two
        conditions. Must be even and >= 2. Populates ``adata.obs["sample"]``
        for testing ``sample_col=``.
    organism
        ``"mouse"`` (default) or ``"human"`` — controls which bundled
        gene-symbol pool ``var_names`` are drawn from.
    seed
        Random seed; the same seed always returns the same object.

    Returns
    -------
    AnnData
        ``.X`` and ``.layers["spliced"]``: mature/spliced counts (float).
        ``.layers["unspliced"]``: nascent/unspliced counts (float).
        ``.obs["condition"]``: ``"Control"`` / ``"Disease"``.
        ``.obs["sample"]``: replicate id, nested within condition.
        ``.var["ground_truth_mechanism"]``: simulated label (teaching only).

    Examples
    --------
    >>> import scatrans as scat
    >>> adata = scat.datasets.load_toy()
    >>> result = scat.partition_de_by_mechanism(
    ...     adata, groupby="condition",
    ...     target_group="Disease", reference_group="Control",
    ... )
    >>> len(result.selected) > 0
    True
    """
    if n_samples < 2 or n_samples % 2 != 0:
        raise ValueError(f"n_samples must be even and >= 2, got {n_samples}")
    if n_cells < 2 * n_samples:
        raise ValueError("n_cells must be >= 2 * n_samples")

    rng = np.random.default_rng(seed)

    half = n_cells // 2
    n_control, n_disease = half, n_cells - half
    conditions = np.array(["Control"] * n_control + ["Disease"] * n_disease)

    samples = np.empty(n_cells, dtype=object)
    for label, mask_condition in (
        ("Control", conditions == "Control"),
        ("Disease", conditions == "Disease"),
    ):
        idx = np.flatnonzero(mask_condition)
        chunks = np.array_split(idx, n_samples // 2)
        for k, chunk in enumerate(chunks):
            samples[chunk] = f"{label.lower()}{k + 1}"

    obs = pd.DataFrame(
        {"condition": conditions, "sample": samples},
        index=[f"cell{i:04d}" for i in range(n_cells)],
    )

    n_txn = max(1, int(round(n_genes * 0.15)))
    n_stab = max(1, int(round(n_genes * 0.15)))
    n_down = max(1, int(round(n_genes * 0.10)))
    n_house = n_genes - n_txn - n_stab - n_down
    if n_house < 0:
        raise ValueError("n_genes too small for the default category split; use n_genes >= 20")

    truth = np.array(
        ["housekeeping"] * n_house
        + ["transcription_up"] * n_txn
        + ["stabilization_up"] * n_stab
        + ["down"] * n_down
    )
    rng.shuffle(truth)

    gene_names = _sample_real_gene_names(organism, n_genes, rng)

    mu_s_base = np.exp(rng.normal(np.log(10.0), 0.7, n_genes)).clip(2.0, 300.0)
    ratio = rng.uniform(0.15, 0.45, n_genes)
    mu_u_base = mu_s_base * ratio

    fold_s = np.ones(n_genes)
    fold_u = np.ones(n_genes)

    txn_mask = truth == "transcription_up"
    fold_s[txn_mask] = rng.uniform(1.8, 3.0, txn_mask.sum())
    fold_u[txn_mask] = fold_s[txn_mask] * rng.uniform(1.5, 2.5, txn_mask.sum())

    stab_mask = truth == "stabilization_up"
    fold_s[stab_mask] = rng.uniform(1.8, 3.0, stab_mask.sum())
    # unspliced (synthesis rate) unchanged: only decay slows down

    down_mask = truth == "down"
    down_fold = rng.uniform(0.3, 0.6, down_mask.sum())
    fold_s[down_mask] = down_fold
    fold_u[down_mask] = down_fold

    def _negbin(mu: np.ndarray, r: float = 4.0) -> np.ndarray:
        mu = np.clip(mu, 1e-3, None)
        p = r / (r + mu)
        return rng.negative_binomial(r, p)

    is_disease = conditions == "Disease"
    spliced = np.empty((n_cells, n_genes))
    unspliced = np.empty((n_cells, n_genes))
    spliced[~is_disease] = _negbin(np.tile(mu_s_base, (int((~is_disease).sum()), 1)))
    unspliced[~is_disease] = _negbin(np.tile(mu_u_base, (int((~is_disease).sum()), 1)))
    spliced[is_disease] = _negbin(np.tile(mu_s_base * fold_s, (int(is_disease.sum()), 1)))
    unspliced[is_disease] = _negbin(np.tile(mu_u_base * fold_u, (int(is_disease.sum()), 1)))

    adata = AnnData(
        X=spliced.astype(float),
        obs=obs,
        var=pd.DataFrame({"ground_truth_mechanism": truth}, index=gene_names),
    )
    adata.layers["spliced"] = spliced.astype(float)
    adata.layers["unspliced"] = unspliced.astype(float)
    adata.uns["scatrans_toy"] = {
        "seed": seed,
        "note": (
            "Synthetic demo dataset: gene symbols are real, but counts and "
            "mechanism assignment are simulated. Not biological data."
        ),
    }
    return adata
