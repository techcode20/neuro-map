"""Integration tests for NeuroMap pipeline outputs.

Verifies that pipeline-generated artifacts exist and are well-formed.
These tests check the output of `python run_all.py` rather than re-running
the full pipeline (which takes several minutes).
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import pandas as pd

EXPECTED_DISEASES = ["Parkinson", "Alzheimer", "Huntington"]
EXPECTED_LIST_TYPES = ["strict", "broad"]
EXPECTED_GENE_COUNTS = {
    "Parkinson": {"strict": 157, "broad": 227},
    "Alzheimer": {"strict": 93, "broad": 297},
    "Huntington": {"strict": 1, "broad": 1},
}


def test_results_directory_exists():
    """results/ directory should exist after running run_all.py."""
    assert os.path.isdir("results"), "results/ directory should exist"


def test_stats_files_exist_for_all_diseases():
    """stats parquet files should exist for every disease/list combination
    that has >=10 genes (Huntington's with 1 gene is underpowered and skips stats)."""
    for disease in EXPECTED_DISEASES:
        for list_type in EXPECTED_LIST_TYPES:
            if disease == "Huntington":
                continue
            path = f"results/stats_{disease.lower()}_{list_type}.parquet"
            assert os.path.isfile(path), f"Missing stats file: {path}"


def test_gene_list_files_exist_for_all_diseases():
    """gene_list parquet files should exist for every disease/list combination."""
    for disease in EXPECTED_DISEASES:
        for list_type in EXPECTED_LIST_TYPES:
            path = f"data/processed/gene_list_{disease.lower()}_{list_type}.parquet"
            assert os.path.isfile(path), f"Missing gene list file: {path}"


def test_stats_files_have_required_columns():
    """stats parquet files should have required columns."""
    for disease in EXPECTED_DISEASES:
        for list_type in EXPECTED_LIST_TYPES:
            if disease == "Huntington":
                continue
            path = f"results/stats_{disease.lower()}_{list_type}.parquet"
            df = pd.read_parquet(path)
            for col in ["score", "pvalue", "fdr", "significant", "depleted"]:
                assert col in df.columns, f"Missing column {col} in {path}"


def test_gene_counts_match_readme():
    """Gene list counts should match README table."""
    for disease in EXPECTED_DISEASES:
        for list_type in EXPECTED_LIST_TYPES:
            path = f"data/processed/gene_list_{disease.lower()}_{list_type}.parquet"
            df = pd.read_parquet(path)
            expected = EXPECTED_GENE_COUNTS[disease][list_type]
            assert len(df) == expected, (
                f"{disease} {list_type}: expected {expected} genes, got {len(df)}"
            )


def test_huntington_stats_not_written_for_underpowered():
    """Huntington's (1 gene) should not have stats files written."""
    for list_type in EXPECTED_LIST_TYPES:
        path = f"results/stats_huntington_{list_type}.parquet"
        # Huntington's stats should NOT exist because <10 genes skips stats writing
        assert not os.path.isfile(path), (
            f"Huntington's {list_type} stats should not exist (underpowered), but {path} exists"
        )


def test_height_negative_control_zero_significant():
    """Height negative control should have 0 FDR-significant regions."""
    path = "results/stats_height_strict.parquet"
    if os.path.isfile(path):
        df = pd.read_parquet(path)
        n_sig = int(df["significant"].sum()) if "significant" in df.columns else 0
        assert n_sig == 0, f"Height negative control should have 0 significant regions, got {n_sig}"


def test_expression_matrix_exists():
    """Expression matrix parquet should exist after pipeline run."""
    path = "data/processed/allen_expression_desikan-killiany.parquet"
    assert os.path.isfile(path), f"Missing expression matrix: {path}"
    df = pd.read_parquet(path)
    assert df.shape[0] > 0, "Expression matrix should have rows"
    assert df.shape[1] > 10000, f"Expression matrix should have >10000 genes, got {df.shape[1]}"


def test_run_manifest_exists():
    """run_manifest.json should exist in results_summary/."""
    path = "results_summary/run_manifest.json"
    assert os.path.isfile(path), f"Missing run manifest: {path}"
    import json
    manifest = json.loads(open(path).read())
    assert "commit" in manifest
    assert "run_date" in manifest
    assert "data_type" in manifest
