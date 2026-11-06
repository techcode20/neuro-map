#!/usr/bin/env python
"""
Run the full NeuroMap pipeline for all diseases.

This script rebuilds all results from raw data.
Usage: python run_all.py
"""

import sys
import hashlib
import json
import argparse
import subprocess
import os
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

import pandas as pd
import numpy as np

from load_data import load_allen_expression, PROCESSED_DIR, PROCESSED_DIR
from gene_lists import build_gene_lists
from scoring import compute_regional_scores, compute_donor_scores, compute_consistency_metric, save_scores
from stats import run_permutation_test, save_statistical_results
from validation import (
    validate_disease,
    load_negative_control_genes,
    run_negative_control_test,
)
from report import generate_html_report

# Configuration
DISEASES = ["Parkinson", "Alzheimer", "Huntington"]
LIST_TYPES = ["strict", "broad"]
MIN_GENES_FOR_ANALYSIS = 10  # Minimum genes needed for meaningful regional analysis
N_RANDOM_SETS = 5000
N_NEGATIVE_CONTROL_SETS = 1000
RANDOM_STATE = 42
ATLAS = "desikan-killiany"

# Striatal positive control gene set (known to be enriched in caudate/putamen)
STRIATAL_GENES = ["DRD1", "DRD2", "PPP1R1B", "ADORA2A", "PENK", "TAC1", "GPR88", "RGS9"]

# Paths for artifacts
RESULTS_DIR = Path("results")
RESULTS_SUMMARY_DIR = Path("results_summary")
RESULTS_DIR.mkdir(exist_ok=True)
RESULTS_SUMMARY_DIR.mkdir(exist_ok=True)

# Striatal positive control gene set (known to be enriched in caudate/putamen)
STRIATAL_GENES = ["DRD1", "DRD2", "PPP1R1B", "ADORA2A", "PENK", "TAC1", "GPR88", "RGS9"]


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

        if len(genes) < MIN_GENES_FOR_ANALYSIS:
            print(
                f"  WARNING: Only {len(genes)} genes for {disease} ({list_type}); "
                f"minimum {MIN_GENES_FOR_ANALYSIS} required for meaningful analysis. "
                f"Skipping statistical analysis (results would be invalid)."
            )
            empty_scores = pd.Series(index=expression_matrix.index, dtype=float, name="score")
            empty_index = expression_matrix.index
            save_scores(empty_scores, disease, list_type)
            # Do NOT write statistical results for underpowered gene lists —
            # writing all-False/all-NaN stats would be misleading.
            continue

        if len(genes) == 0:
            print(
                f"  WARNING: No matched genes for {disease} ({list_type}); "
                f"skipping scoring/stats (real GWAS data not available)."
            )
            empty_scores = pd.Series(index=expression_matrix.index, dtype=float, name="score")
            empty_index = expression_matrix.index
            save_scores(empty_scores, disease, list_type)
            # Do NOT write statistical results when no genes match —
            # writing empty stats would be misleading.
            continue

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


def run_striatal_positive_control(
    expression_matrix: pd.DataFrame,
    donor_matrices: dict,
    n_sets: int = N_RANDOM_SETS,
    random_state: int = RANDOM_STATE,
) -> dict:
    """Run striatal positive control with known striatal gene set.
    
    The striatal gene set (DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9)
    is known to be enriched in caudate and putamen. This serves as a positive
    control for the regional pipeline.
    """
    print(f"\n{'='*60}")
    print("Running Striatal Positive Control...")
    print(f"{'='*60}")
    print(f"Striatal gene set ({len(STRIATAL_GENES)} genes): {STRIATAL_GENES}")

    # Create gene DataFrame in the expected format
    striatal_df = pd.DataFrame({"gene_symbol": STRIATAL_GENES})
    matched_genes = [g for g in STRIATAL_GENES if g in expression_matrix.columns]
    print(f"Matched {len(matched_genes)}/{len(STRIATAL_GENES)} genes to expression matrix")
    
    if len(matched_genes) == 0:
        print("WARNING: No striatal genes matched expression matrix")
        return {}

    matched_df = pd.DataFrame({"gene_symbol": matched_genes})

    # Compute regional scores
    print("  Computing regional scores...")
    regional_scores = compute_regional_scores(expression_matrix, matched_df)

    # Donor-level scores for consistency
    print("  Computing donor-level scores...")
    donor_scores = compute_donor_scores(donor_matrices, matched_df)
    consistency = compute_consistency_metric(donor_scores)

    save_scores(regional_scores, "StriatalControl", "strict", donor_scores, consistency)

    # Permutation test
    print("  Running permutation test...")
    observed, pvals, fdr_results, effects = run_permutation_test(
        expression_matrix,
        matched_df,
        n_sets=n_sets,
        random_state=random_state,
    )

    save_statistical_results(
        "StriatalControl", "strict", observed, pvals, fdr_results, effects
    )

    # Validation
    print("  Running validation...")
    stats_results = pd.DataFrame({
        "significant": fdr_results["significant"],
        "score": observed,
    })
    validation = validate_disease(
        disease="StriatalControl",
        stats_results=stats_results,
        donor_consistency=consistency,
        expression_matrix=expression_matrix,
        strict_genes=matched_df,
        broad_genes=matched_df,
        list_type="strict",
    )

    # Report key results for caudate/putamen
    caudate_putamen = [36, 37, 77, 78]  # L/R caudate/putamen region IDs
    cp_results = fdr_results[fdr_results.index.isin(caudate_putamen)].copy()
    cp_results["score"] = observed[caudate_putamen]
    cp_results["cohens_d"] = effects["cohens_d"][caudate_putamen]
    print(f"\n  Caudate/Putamen results:")
    for _, row in cp_results.iterrows():
        sig = "***" if row["significant"] else ""
        print(f"    Region {int(row.name)}: score={row['score']:.4f}, "
              f"p={row['pvalue']:.4g}, FDR={row['fdr']:.4g}, "
              f"cohens_d={row['cohens_d']:.3f}, significant={row['significant']} {sig}")

    return {
        "scores": regional_scores,
        "donor_scores": donor_scores,
        "consistency": consistency,
        "observed": observed,
        "pvals": pvals,
        "fdr_results": fdr_results,
        "effects": effects,
        "validation": validation,
    }


def main():
    parser = argparse.ArgumentParser(description="Run the full NeuroMap pipeline.")
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete data/processed/ and results/ before running.",
    )
    parser.add_argument(
        "--clean-all",
        action="store_true",
        help="Delete data/processed/, results/, and results_summary/ before running.",
    )
    args = parser.parse_args()

    if args.clean:
        import shutil

        for d in [PROCESSED_DIR, Path("results")]:
            if d.exists():
                shutil.rmtree(d)
                print(f"Cleaned {d}")
        RESULTS_DIR.mkdir(exist_ok=True)
        RESULTS_SUMMARY_DIR.mkdir(exist_ok=True)
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    if args.clean_all:
        import shutil

        for d in [PROCESSED_DIR, Path("results"), Path("results_summary")]:
            if d.exists():
                shutil.rmtree(d)
                print(f"Cleaned {d}")
        RESULTS_DIR.mkdir(exist_ok=True)
        RESULTS_SUMMARY_DIR.mkdir(exist_ok=True)
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print("NeuroMap Pipeline")
    print("=" * 60)

    # Step 1: Load Allen expression data
    print("\n1. Loading Allen Human Brain Atlas data...")
    expression_matrix, donor_matrices = load_allen_expression(atlas_name=ATLAS)

    # Mock-data guard: prevent writing results with synthetic data
    cols = [str(c) for c in expression_matrix.columns[:20]]
    is_real = expression_matrix.shape[1] > 10000 and not any(
        str(c).startswith("GENE_") for c in cols
    )
    if not is_real and not os.environ.get("ALLOW_MOCK"):
        sys.exit(
            "ABORT: expression matrix is mock data. "
            "Set REAL_DATA=1 and place real AHBA data in data/raw/allen/. "
            "Nothing written to results/."
        )

    # Detect mock data and warn
    is_mock = (
        expression_matrix.shape[1] < 10000
        or any(str(c).startswith("GENE_") for c in expression_matrix.columns[:20])
    )
    if is_mock:
        print()
        print("=" * 60)
        print("WARNING: MOCK DATA DETECTED")
        print("=" * 60)
        print("The expression matrix contains synthetic/mock gene symbols "
              f"({expression_matrix.shape[1]} genes).")
        print("No real Allen Human Brain Atlas data was loaded.")
        print()
        print("To use real AHBA data:")
        print("  1. Download microarray data from https://human.brain-map.org")
        print("  2. Place donor folders in data/raw/allen/")
        print("  3. Run: REAL_DATA=1 python run_all.py")
        print()
        print("NOTE: Running in mock mode. Results are NOT valid for scientific use.")
        print("      README.md and results_summary/ will NOT be updated.")
        print("=" * 60)
        print()

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

    # Save negative control stats for audit/reporting
    save_statistical_results(
        "Height",
        "strict",
        negative_control_results["score"],
        negative_control_results["pvalue"],
        negative_control_results[["fdr", "significant", "pvalue_two_sided", "fdr_two_sided", "depleted"]]
        if "pvalue_two_sided" in negative_control_results.columns
        else pd.DataFrame(
            index=negative_control_results.index,
            data={
                "fdr": negative_control_results["fdr"],
                "significant": negative_control_results["significant"],
                "pvalue_two_sided": pd.Series(dtype=float),
                "fdr_two_sided": pd.Series(dtype=float),
                "depleted": pd.Series(dtype=bool),
            },
        ),
        negative_control_results[["cohens_d", "null_mean", "null_std"]]
        if "cohens_d" in negative_control_results.columns
        else pd.DataFrame(
            index=negative_control_results.index,
            data={
                "cohens_d": pd.Series(dtype=float),
                "null_mean": pd.Series(dtype=float),
                "null_std": pd.Series(dtype=float),
            },
        ),
    )

    # Step 2b: Striatal positive control
    striatal_results = run_striatal_positive_control(
        expression_matrix, donor_matrices
    )

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

    # Write manifest to results_summary
    import hashlib
    import datetime

    manifest_dir = RESULTS_SUMMARY_DIR
    manifest_dir.mkdir(exist_ok=True)

    # Compute checksum of expression matrix
    expr_cache = PROCESSED_DIR / f"allen_expression_{ATLAS}.parquet"
    expr_checksum = None
    if expr_cache.exists():
        h = hashlib.sha256()
        with open(expr_cache, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                h.update(chunk)
        expr_checksum = h.hexdigest()

    manifest = {
        "commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "gwas_file": str(Path("data/raw/gwas_catalog_full.tsv").resolve()),
        "gwas_file_size": Path("data/raw/gwas_catalog_full.tsv").stat().st_size,
        "expression_matrix_checksum": expr_checksum,
        "expression_matrix_shape": list(expression_matrix.shape),
        "data_type": "mock (synthetic)" if is_mock else "real (Allen Human Brain Atlas)",
        "run_date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    manifest_path = manifest_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nManifest written to {manifest_path}")


if __name__ == "__main__":
    main()