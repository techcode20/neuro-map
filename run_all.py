#!/usr/bin/env python
"""
Run the full NeuroMap pipeline for all diseases.

This script rebuilds all results from raw data.
Usage: python run_all.py
"""

import pandas as pd
import numpy as np
from pathlib import Path

from src.load_data import load_allen_expression
from src.gene_lists import build_gene_lists
from src.scoring import compute_regional_scores, compute_donor_scores, compute_consistency_metric, save_scores
from src.stats import run_permutation_test, save_statistical_results
from src.validation import (
    validate_disease,
    load_negative_control_genes,
    run_negative_control_test,
)
from src.report import generate_html_report

# Configuration
DISEASES = ["Parkinson", "Alzheimer", "Huntington"]
LIST_TYPES = ["strict", "broad"]
N_RANDOM_SETS = 5000
N_NEGATIVE_CONTROL_SETS = 1000
RANDOM_STATE = 42
ATLAS = "desikan-killiany"


def run_pipeline_for_disease(
    disease: str,
    expression_matrix: pd.DataFrame,
    donor_matrices: dict,
    negative_control_results: pd.DataFrame = None,
) -> dict:
    """Run the complete pipeline for one disease."""
    print(f"\n{'='*60}")
    print(f"Processing {disease}...")
    print(f"{'='*60}")

    results = {}

    # Build gene lists
    print("\n1. Building gene lists...")
    strict_genes, broad_genes = build_gene_lists(
        disease, expression_matrix, force_refresh=False
    )

    for list_type, genes in [("strict", strict_genes), ("broad", broad_genes)]:
        print(f"\n  Processing {list_type} list ({len(genes)} genes)...")

        # Compute regional scores
        print("  Computing regional scores...")
        regional_scores = compute_regional_scores(expression_matrix, genes)

        # Donor-level scores for consistency
        print("  Computing donor-level scores...")
        donor_scores = compute_donor_scores(donor_matrices, genes)
        consistency = compute_consistency_metric(donor_scores)

        save_scores(regional_scores, disease, list_type, donor_scores, consistency)

        # Permutation test
        print("  Running permutation test...")
        observed, pvals, fdr_results, effects = run_permutation_test(
            expression_matrix,
            genes,
            n_sets=N_RANDOM_SETS,
            random_state=RANDOM_STATE,
        )

        save_statistical_results(
            disease, list_type, observed, pvals, fdr_results, effects
        )

        # Validation
        print("  Running validation...")
        stats_results = pd.DataFrame({
            "significant": fdr_results["significant"],
            "score": observed,
        })
        validation = validate_disease(
            disease=disease,
            stats_results=stats_results,
            donor_consistency=consistency,
            expression_matrix=expression_matrix,
            strict_genes=strict_genes,
            broad_genes=broad_genes,
            list_type=list_type,
            negative_control_results=negative_control_results,
        )

        results[f"{list_type}_genes"] = genes
        results[f"{list_type}_scores"] = regional_scores
        results[f"{list_type}_donor_scores"] = donor_scores
        results[f"{list_type}_consistency"] = consistency
        results[f"{list_type}_stats"] = pd.DataFrame({
            "score": observed,
            "pvalue": pvals,
            "fdr": fdr_results["fdr"],
            "significant": fdr_results["significant"],
            "cohens_d": effects["cohens_d"],
        })
        results[f"{list_type}_validation"] = validation

        # Generate report
        print("  Generating report...")
        generate_html_report(disease, list_type)

    return results


def main():
    print("NeuroMap Pipeline")
    print("=" * 60)

    # Step 1: Load Allen expression data
    print("\n1. Loading Allen Human Brain Atlas data...")
    expression_matrix, donor_matrices = load_allen_expression(atlas=ATLAS)

    print(f"\nExpression matrix shape: {expression_matrix.shape}")
    print(f"Number of donors: {len(donor_matrices)}")
    print(f"Regions: {expression_matrix.index.tolist()[:10]}...")

    # Step 2: Negative control (Height - no brain relevance)
    print("\n2. Running negative control (Height)...")
    control_genes = load_negative_control_genes(expression_matrix, "Height")
    negative_control_results = run_negative_control_test(
        expression_matrix, control_genes, n_sets=N_NEGATIVE_CONTROL_SETS
    )
    print(f"  Negative control significant regions: {negative_control_results['significant'].sum()}")

    # Step 3: Run pipeline for each disease
    all_results = {}
    for disease in DISEASES:
        all_results[disease] = run_pipeline_for_disease(
            disease, expression_matrix, donor_matrices, negative_control_results
        )

    # Step 4: Generate comparison reports (placeholder)
    print("\n4. Generating cross-disease comparison...")
    # TODO: Implement disease comparison

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print("=" * 60)
    print("\nResults saved to:")
    print("  - data/processed/  (expression matrices, gene lists)")
    print("  - results/         (scores, stats, validation, reports)")
    print("\nTo view dashboard: python -m streamlit run app.py")


if __name__ == "__main__":
    main()