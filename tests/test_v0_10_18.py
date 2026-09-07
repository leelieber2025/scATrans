"""Regression tests for cross-group comparison edge cases."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scatrans.tl.compare import compare_de_across_groups


def test_nullable_group_column(monkeypatch):
    import scatrans.tl.compare as module

    data = ad.AnnData(np.ones((4, 2)))
    data.obs["cell_type"] = pd.array(["T", "T", pd.NA, "B"], dtype="string")
    seen = []

    def de(subset, **kwargs):
        seen.append(subset.n_obs)
        return subset, pd.DataFrame(
            {"logFC": [2.0, -2.0], "p_adj": [0.001, 0.001]}, index=["up", "down"]
        )

    monkeypatch.setattr(module, "differential_expression", de)
    result = compare_de_across_groups(
        data, "cell_type", "condition", "A", "B", groups=["T", "B"], min_cells=1
    )
    assert seen == [2, 1]
    assert result.failed == {}
    assert result.summary["total"].tolist() == [2, 2]


@pytest.mark.parametrize("groups", [[1, "1"], ["T", "T"]])
def test_ambiguous_group_names_rejected(groups):
    data = ad.AnnData(np.ones((2, 1)))
    data.obs["cell_type"] = ["T", "B"]
    with pytest.raises(ValueError, match="unique string"):
        compare_de_across_groups(data, "cell_type", "condition", "A", "B", groups=groups)
