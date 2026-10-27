"""Statistical testing: permutation tests, FDR correction, effect sizes."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, Optional, List
from scipy import stats
from statsmodels.stats.multitest import multipletests

from load_data import zscore_expression
from scoring import compute_regional_scores

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def generate_matched_random_gene_sets(
    expression_matrix: pd.DataFrame,
    disease_genes: pd.DataFrame,
    gene_column: str = "gene_symbol",
    n_sets: int = 5000,
    n_bins: int = 20,
    random_state: int = 42,
) -> List[List[str]]:
    """
    Generate random gene sets matched on mean expression level.

    Parameters
    ----------
    expression_matrix : pd.DataFrame
        Region x gene expression matrix.
    disease_genes : pd.DataFrame
        Disease gene list.
    gene_column : str
        Gene symbol column.
    n_sets : int
        Number of random gene sets to generate.
    n_bins : int
        Number of expression-level bins for matching.
    random_state : int
        Random seed for reproducibility.

    Returns
    -------
    list of list of str
        Each inner list is a random gene set of same size as disease_genes.
    """
    np.random.seed(random_state)

    disease_symbols = set(disease_genes[gene_column].tolist())
    all_genes = set(expression_matrix.columns)
    background_genes = list(all_genes - disease_symbols)

    # Compute mean expression per gene (across regions)
    mean_expr = expression_matrix[background_genes].mean(axis=0)

    # Bin background genes by mean expression
    bins = pd.qcut(mean_expr, q=n_bins, labels=False, duplicates="drop")
    bin_labels = bins.unique()

    # Get bin assignments for disease genes
    disease_mean_expr = expression_matrix[list(disease_symbols & all_genes)].mean(axis=0)
    disease_bins = pd.qcut(disease_mean_expr, q=n_bins, labels=False, duplicates="drop")

    # Count disease genes per bin
    bin_counts = disease_bins.value_counts().to_dict()

    # Generate random sets
    random_sets = []
    for _ in range(n_sets):
        random_set = []
        for bin_label, count in bin_counts.items():
            bin_genes = bins[bins == bin_label].index.tolist()
            if len(bin_genes) >= count:
                selected = np.random.choice(bin_genes, size=count, replace=False)
                random_set.extend(selected.tolist())
        random_sets.append(random_set)

    return random_sets


def compute_null_distribution(
    expression_matrix: pd.DataFrame,
    random_gene_sets: List[List[str]],
    zscore_axis: int = 0,
) -> pd.DataFrame:
    """
    Compute regional scores for all random gene sets.

    Parameters
    ----------
    expression_matrix : pd.DataFrame
        Region x gene expression matrix.
    random_gene_sets : list of list of str
        Random gene sets.
    zscore_axis : int
        Axis for z-scoring.

    Returns
    -------
    pd.DataFrame
        Regions x random_sets matrix of null scores.
    """
    # Z-score once
    zscored = zscore_expression(expression_matrix, axis=zscore_axis)

    # Compute scores for each random set
    null_scores = {}
    for i, gene_set in enumerate(random_gene_sets):
        available = [g for g in gene_set if g in zscored.columns]
        if len(available) > 0:
            scores = zscored[available].mean(axis=1)
            # Ensure it's a Series with region index
            null_scores[i] = scores

    null_df = pd.DataFrame(null_scores)
    null_df.index.name = "region"

    return null_df


def compute_empirical_pvalues(
    observed_scores: pd.Series,
    null_scores: pd.DataFrame,
    alternative: str = "greater",
) -> pd.Series:
    """
    Compute empirical p-values from null distribution.

    Parameters
    ----------
    observed_scores : pd.Series
        Observed regional scores.
    null_scores : pd.DataFrame
        Null distribution (regions x random_sets).
    alternative : str
        'greater' for one-sided test of enrichment.

    Returns
    -------
    pd.Series
        Empirical p-values per region.
    """
    n_sets = null_scores.shape[1]
    pvalues = {}

    for region in observed_scores.index:
        if region not in null_scores.index:
            pvalues[region] = np.nan
            continue

        observed = observed_scores[region]
        # Get null values for this region, ensure 1D array
        null_vals = null_scores.loc[region].values
        # Flatten in case it's 2D
        null_vals = np.asarray(null_vals).ravel()

        if alternative == "greater":
            # p = (1 + number of null >= observed) / (1 + n_sets)
            n_extreme = np.sum(null_vals >= observed)
        else:
            n_extreme = np.sum(null_vals <= observed)

        pvalues[region] = (1 + n_extreme) / (1 + n_sets)

    return pd.Series(pvalues, name="pvalue", index=observed_scores.index)


def apply_fdr_correction(
    pvalues: pd.Series,
    alpha: float = 0.05,
    method: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvalues : pd.Series
        Raw p-values.
    alpha : float
        FDR threshold.
    method : str
        Correction method ('fdr_bh' for Benjamini-Hochberg).

    Returns
    -------
    pd.DataFrame
        DataFrame with pvalue, fdr, and reject columns.
    """
    valid = pvalues.dropna()
    reject, fdr_corrected, _, _ = multipletests(
        valid.values, alpha=alpha, method=method
    )

    result = pd.DataFrame({
        "pvalue": pvalues,
        "fdr": np.nan,
        "significant": False,
    })
    result.loc[valid.index, "fdr"] = fdr_corrected
    result.loc[valid.index, "significant"] = reject

    return result


def compute_effect_sizes(
    observed_scores: pd.Series,
    null_scores: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute effect sizes relative to null distribution.

    Returns observed score, null mean, null std, and Cohen's d.
    """
    effect_data = []

    for region in observed_scores.index:
        if region not in null_scores.index:
            continue

        observed = observed_scores[region]
        null_vals = null_scores.loc[region].values

        null_mean = np.mean(null_vals)
        null_std = np.std(null_vals, ddof=1)

        if null_std > 0:
            cohens_d = (observed - null_mean) / null_std
        else:
            cohens_d = np.nan

        effect_data.append({
            "region": region,
            "observed_score": observed,
            "null_mean": null_mean,
            "null_std": null_std,
            "cohens_d": cohens_d,
        })

    return pd.DataFrame(effect_data).set_index("region")


def run_permutation_test(
    expression_matrix: pd.DataFrame,
    disease_genes: pd.DataFrame,
    gene_column: str = "gene_symbol",
    n_sets: int = 5000,
    n_bins: int = 20,
    random_state: int = 42,
    zscore_axis: int = 0,
) -> Tuple[pd.Series, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run full permutation test pipeline.

    Returns
    -------
    observed_scores : pd.Series
    pvalues : pd.Series
    fdr_results : pd.DataFrame
    effect_sizes : pd.DataFrame
    """
    print("Computing observed scores...")
    observed = compute_regional_scores(
        expression_matrix, disease_genes, gene_column, zscore_axis
    )

    print(f"Generating {n_sets} matched random gene sets...")
    random_sets = generate_matched_random_gene_sets(
        expression_matrix, disease_genes, gene_column, n_sets, n_bins, random_state
    )

    print("Computing null distribution...")
    null_dist = compute_null_distribution(expression_matrix, random_sets, zscore_axis)

    print("Computing empirical p-values...")
    pvals = compute_empirical_pvalues(observed, null_dist)

    print("Applying FDR correction...")
    fdr_results = apply_fdr_correction(pvals)

    print("Computing effect sizes...")
    effects = compute_effect_sizes(observed, null_dist)

    return observed, pvals, fdr_results, effects


def save_statistical_results(
    disease: str,
    list_type: str,
    observed_scores: pd.Series,
    pvalues: pd.Series,
    fdr_results: pd.DataFrame,
    effect_sizes: pd.DataFrame,
) -> None:
    """Save all statistical results."""
    # Combined results
    results = pd.DataFrame({
        "score": observed_scores,
        "pvalue": pvalues,
        "fdr": fdr_results["fdr"],
        "significant": fdr_results["significant"],
        "cohens_d": effect_sizes["cohens_d"],
        "null_mean": effect_sizes["null_mean"],
        "null_std": effect_sizes["null_std"],
    })
    results.index.name = "region"
    results["disease"] = disease
    results["list_type"] = list_type

    file = RESULTS_DIR / f"stats_{disease.lower()}_{list_type}.parquet"
    results.to_parquet(file)
    print(f"Saved statistical results to {file}")