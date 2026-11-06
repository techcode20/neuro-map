"""Unit tests for regional scoring (src/scoring.py)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
import pandas as pd
from scoring import compute_regional_scores, compute_donor_scores, compute_consistency_metric


def test_compute_regional_scores_returns_series():
    """compute_regional_scores should return a pd.Series indexed by region."""
    expr = pd.DataFrame(
        {
            "GENE_A": [1.0, 2.0, 3.0, 4.0],
            "GENE_B": [4.0, 3.0, 2.0, 1.0],
            "GENE_C": [2.0, 2.0, 2.0, 2.0],
        },
        index=["Region1", "Region2", "Region3", "Region4"],
    )
    genes = pd.DataFrame({"gene_symbol": ["GENE_A", "GENE_B"]})

    scores = compute_regional_scores(expr, genes)
    assert isinstance(scores, pd.Series)
    assert list(scores.index) == ["Region1", "Region2", "Region3", "Region4"]


def test_compute_regional_scores_empty_genes():
    """Empty gene list should return series of NaN/zeros, not crash."""
    expr = pd.DataFrame(
        {"GENE_A": [1.0, 2.0]}, index=["Region1", "Region2"]
    )
    genes = pd.DataFrame({"gene_symbol": []})

    scores = compute_regional_scores(expr, genes)
    assert isinstance(scores, pd.Series)
    assert len(scores) == 2


def test_compute_donor_scores_returns_dataframe():
    """compute_donor_scores should return a DataFrame indexed by region with donor columns."""
    expr = pd.DataFrame(
        {
            "GENE_A": [1.0, 2.0, 3.0, 4.0],
            "GENE_B": [4.0, 3.0, 2.0, 1.0],
        },
        index=["Region1", "Region2", "Region3", "Region4"],
    )
    donor_matrices = {
        "Donor1": expr.copy(),
        "Donor2": expr.copy() * 0.5,
    }
    genes = pd.DataFrame({"gene_symbol": ["GENE_A", "GENE_B"]})

    donor_scores = compute_donor_scores(donor_matrices, genes)
    assert isinstance(donor_scores, pd.DataFrame)
    assert "Donor1" in donor_scores.columns
    assert "Donor2" in donor_scores.columns
    assert list(donor_scores.index) == ["Region1", "Region2", "Region3", "Region4"]


def test_compute_consistency_metric_returns_series():
    """compute_consistency_metric should return a Series of positive donor counts."""
    donor_scores = pd.DataFrame(
        {
            "Donor1": [1.0, 2.0, 3.0],
            "Donor2": [1.1, 1.9, 3.1],
            "Donor3": [0.9, 2.1, 2.9],
        },
        index=["R1", "R2", "R3"],
    )
    consistency = compute_consistency_metric(donor_scores)
    assert isinstance(consistency, pd.Series)
    assert (consistency >= 0).all()
    assert (consistency <= 3).all()
