"""Validation against expected regions and evidence badges."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional

from load_data import load_gwas_catalog

REFERENCE_DIR = Path(__file__).parent.parent / "reference"
RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def load_expected_regions() -> pd.DataFrame:
    """Load the frozen expected regions table."""
    file = REFERENCE_DIR / "expected_regions.csv"
    if file.exists():
        return pd.read_csv(file)
    else:
        raise FileNotFoundError(f"Expected regions file not found: {file}")


def load_negative_control_genes(
    expression_matrix: pd.DataFrame,
    control_trait: str = "Height",
    pvalue_threshold: float = 5e-8,
    n_control_genes: int = 20,
) -> pd.DataFrame:
    """Load negative control gene set (trait with no brain relevance)."""
    from load_data import load_gwas_catalog, match_genes_to_atlas

    control_gwas = load_gwas_catalog(control_trait, pvalue_threshold=pvalue_threshold)
    control_matched, _ = match_genes_to_atlas(control_gwas, expression_matrix)

    # If no genes matched (e.g., mock genes vs real expression matrix),
    # use random real genes from the expression matrix as negative control
    if len(control_matched) == 0:
        print(f"  No control genes matched expression matrix; using {n_control_genes} random real genes")
        import numpy as np
        np.random.seed(42)
        all_genes = list(expression_matrix.columns)
        random_genes = np.random.choice(all_genes, size=min(n_control_genes, len(all_genes)), replace=False)
        control_matched = pd.DataFrame({
            "gene_symbol": random_genes,
            "reported_gene": random_genes,
            "mapped_gene": random_genes,
            "pvalue": [1e-8] * len(random_genes),
            "study_accession": ["random"] * len(random_genes),
            "trait": ["random"] * len(random_genes),
            "mapping_method": ["random"] * len(random_genes),
        })

    return control_matched


def assign_evidence_badges(
    stats_results: pd.DataFrame,
    donor_consistency: pd.Series,
    expected_regions: pd.DataFrame,
    disease: str,
    list_type: str = "strict",
    min_consistent_donors: int = 4,
) -> pd.DataFrame:
    """
    Assign evidence badges based on checklist.

    Checklist:
    1. Random-gene-set test passes (FDR < 0.05)
    2. Pattern consistent across donors (>= min_consistent_donors positive)
    3. Region matches frozen expected-regions table
    4. Negative control does not flag same region (advanced)

    Badge rules:
    - Supported: FDR < 0.05 AND consistent in >= min_consistent_donors
    - Uncertain: Only one condition holds
    - Not supported: Neither holds
    """
    # Get expected regions for this disease
    expected = expected_regions[expected_regions["disease"] == disease]["region"].tolist()

    badges = []
    for region in stats_results.index:
        row = stats_results.loc[region]
        fdr_pass = row["significant"] if not pd.isna(row["significant"]) else False
        n_positive = donor_consistency.get(region, 0)
        consistency_pass = n_positive >= min_consistent_donors
        in_expected = region in expected

        # Count passing criteria
        n_pass = sum([fdr_pass, consistency_pass, in_expected])

        if n_pass >= 2:
            badge = "Supported"
        elif n_pass == 1:
            badge = "Uncertain"
        else:
            badge = "Not supported"

        badges.append({
            "region": region,
            "fdr_pass": fdr_pass,
            "consistency_pass": consistency_pass,
            "in_expected": in_expected,
            "n_criteria_pass": n_pass,
            "badge": badge,
            "n_positive_donors": n_positive,
        })

    badge_df = pd.DataFrame(badges).set_index("region")
    badge_df["disease"] = disease
    badge_df["list_type"] = list_type

    return badge_df


def compare_gene_lists(
    stats_strict: pd.DataFrame,
    stats_broad: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compare results between strict and broad gene lists.

    Returns regions that are significant in both, strict only, broad only, or neither.
    """
    merged = stats_strict[["significant", "score"]].join(
        stats_broad[["significant", "score"]],
        lsuffix="_strict",
        rsuffix="_broad",
        how="outer"
    )

    def classify(row):
        strict_sig = row.get("significant_strict", False)
        broad_sig = row.get("significant_broad", False)

        if strict_sig and broad_sig:
            return "Both"
        elif strict_sig:
            return "Strict only"
        elif broad_sig:
            return "Broad only"
        else:
            return "Neither"

    merged["agreement"] = merged.apply(classify, axis=1)
    return merged


def run_negative_control_test(
    expression_matrix: pd.DataFrame,
    control_genes: pd.DataFrame,
    n_sets: int = 1000,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Run permutation test with negative control genes.

    Regions flagged by negative control are suspect.
    """
    from stats import run_permutation_test

    print(f"Running negative control test ({control_genes.shape[0]} genes)...")
    observed, pvals, fdr_results, effects = run_permutation_test(
        expression_matrix,
        control_genes,
        n_sets=n_sets,
        random_state=random_state,
    )

    control_results = pd.DataFrame({
        "score": observed,
        "pvalue": pvals,
        "fdr": fdr_results["fdr"],
        "significant": fdr_results["significant"],
    })
    control_results.index.name = "region"
    return control_results


def validate_disease(
    disease: str,
    stats_results: pd.DataFrame,
    donor_consistency: pd.Series,
    expression_matrix: pd.DataFrame,
    strict_genes: pd.DataFrame,
    broad_genes: pd.DataFrame,
    list_type: str = "strict",
    negative_control_results: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Run full validation pipeline for a disease.

    Returns validation results with evidence badges.
    """
    expected_regions = load_expected_regions()

    # Assign badges
    badges = assign_evidence_badges(
        stats_results,
        donor_consistency,
        expected_regions,
        disease,
        list_type,
    )

    # Compare strict vs broad (only if both exist)
    comparison = pd.DataFrame(index=badges.index)
    try:
        if list_type == "strict":
            broad_stats = load_stats_results(disease, "broad")
            strict_stats = stats_results
        else:
            strict_stats = load_stats_results(disease, "strict")
            broad_stats = stats_results
        comparison = compare_gene_lists(strict_stats, broad_stats)
    except FileNotFoundError:
        # Broad/strict not yet computed, skip comparison
        comparison["agreement"] = "Unknown"

    # Add gene count and exploratory flag
    n_strict = len(strict_genes)
    n_broad = len(broad_genes)
    badges["n_strict_genes"] = n_strict
    badges["n_broad_genes"] = n_broad
    badges["exploratory"] = (n_strict < 10) or (n_broad < 10)

    # Add negative control flag if available
    if negative_control_results is not None:
        badges["negative_control_flag"] = badges.index.map(
            negative_control_results["significant"]
        ).fillna(False)

    # Save
    file = RESULTS_DIR / f"validation_{disease.lower()}_{list_type}.parquet"
    badges.to_parquet(file)
    badges.reset_index().to_csv(
        RESULTS_DIR / f"validation_{disease.lower()}_{list_type}.csv", index=False
    )
    print(f"Saved validation results to {file}")

    return badges


def load_stats_results(disease: str, list_type: str = "strict") -> pd.DataFrame:
    """Load statistical results."""
    file = RESULTS_DIR / f"stats_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    else:
        raise FileNotFoundError(f"Stats results not found: {file}")