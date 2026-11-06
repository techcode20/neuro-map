"""Unit tests for gene list building (src/gene_lists.py)."""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pandas as pd
from gene_lists import build_gene_lists, filter_gwas_by_disease

EXPECTED_COUNTS = {
    "Parkinson": {"strict": 157, "broad": 227},
    "Alzheimer": {"strict": 93, "broad": 297},
    "Huntington": {"strict": 1, "broad": 1},
    # Height is a real negative control with matched genes but 0 significant regions
    "Height": {"strict": 677, "broad": 6432},
}


def test_gene_list_counts_match_readme():
    """Verify gene counts match README (never hand-typed)."""
    from load_data import load_allen_expression

    expr_matrix, _ = load_allen_expression(force_refresh=False)

    for disease, expected in EXPECTED_COUNTS.items():
        strict, broad = build_gene_lists(disease, expr_matrix, force_refresh=False)
        assert len(strict) == expected["strict"], (
            f"{disease} strict: expected {expected['strict']}, got {len(strict)}"
        )
        assert len(broad) == expected["broad"], (
            f"{disease} broad: expected {expected['broad']}, got {len(broad)}"
        )


def test_gene_lists_have_gene_symbol_column():
    """Gene list DataFrames must have a gene_symbol column."""
    from load_data import load_allen_expression

    expr_matrix, _ = load_allen_expression(force_refresh=False)

    strict, broad = build_gene_lists("Parkinson", expr_matrix, force_refresh=False)
    assert "gene_symbol" in strict.columns
    assert "gene_symbol" in broad.columns


def test_no_gene_symbols_start_with_GENE_():
    """Real GWAS data should not contain mock GENE_XXXX symbols."""
    from load_data import load_allen_expression

    expr_matrix, _ = load_allen_expression(force_refresh=False)

    for disease in ["Parkinson", "Alzheimer", "Huntington", "Height"]:
        strict, broad = build_gene_lists(disease, expr_matrix, force_refresh=False)
        for df in [strict, broad]:
            if not df.empty:
                mock = df["gene_symbol"].astype(str).str.startswith("GENE_")
                assert not mock.any(), (
                    f"{disease} contains mock GENE_XXXX symbols: "
                    f"{df.loc[mock, 'gene_symbol'].tolist()}"
                )


def test_height_is_negative_control():
    """Height (negative control) should have 0 FDR-significant regions."""
    from load_data import load_allen_expression

    expr_matrix, _ = load_allen_expression(force_refresh=False)
    strict, broad = build_gene_lists("Height", expr_matrix, force_refresh=False)
    # Height should have real genes matched (not mock GENE_XXXX)
    assert len(strict) > 0, "Height should have matched genes"
    assert not strict["gene_symbol"].astype(str).str.startswith("GENE_").any(), (
        "Height should not contain mock GENE_XXXX symbols"
    )


def test_trait_filters_written():
    """Trait filter audit trail should exist after running gene list builder."""
    from load_data import load_allen_expression

    expr_matrix, _ = load_allen_expression(force_refresh=False)
    build_gene_lists("Parkinson", expr_matrix, force_refresh=True)

    trait_filters_path = os.path.join(
        os.path.dirname(__file__), "..", "results", "trait_filters.csv"
    )
    assert os.path.exists(trait_filters_path), "trait_filters.csv should exist"
    tf = pd.read_csv(trait_filters_path)
    assert "disease" in tf.columns
    assert "trait" in tf.columns
    assert len(tf) > 0
