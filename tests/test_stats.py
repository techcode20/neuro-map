"""Unit tests for statistical testing (src/stats.py)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd
from stats import apply_fdr_correction, compute_empirical_pvalues


def test_apply_fdr_correction_known_values():
    """FDR correction should flag small p-values as significant."""
    pvals = pd.Series(
        [0.001, 0.01, 0.03, 0.2, 0.5, 0.8],
        index=["R1", "R2", "R3", "R4", "R5", "R6"],
    )
    result = apply_fdr_correction(pvals, alpha=0.05)
    assert "significant" in result.columns
    assert "fdr" in result.columns
    assert bool(result.loc["R1", "significant"]) is True
    assert bool(result.loc["R6", "significant"]) is False


def test_apply_fdr_correction_all_nan():
    """All-NaN p-values should not crash FDR correction."""
    pvals = pd.Series([np.nan, np.nan, np.nan], index=["R1", "R2", "R3"])
    result = apply_fdr_correction(pvals, alpha=0.05)
    assert result["significant"].sum() == 0
    assert result["fdr"].isna().all()


def test_compute_empirical_pvalues_bounds():
    """Empirical p-values should be in [0, 1]."""
    observed = pd.Series([1.0, 2.0, 3.0], index=["R1", "R2", "R3"])
    null = pd.DataFrame(
        {
            "set1": [0.5, 1.5, 2.5],
            "set2": [1.0, 2.0, 3.0],
            "set3": [1.5, 2.5, 3.5],
        },
        index=["R1", "R2", "R3"],
    )
    pvals = compute_empirical_pvalues(observed, null)
    assert (pvals >= 0).all()
    assert (pvals <= 1).all()
