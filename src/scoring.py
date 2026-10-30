"""Regional scoring for disease-gene expression enrichment."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict

from load_data import zscore_expression

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def compute_regional_scores(
    expression_matrix: pd.DataFrame,
    disease_genes: pd.DataFrame,
    gene_column: str = "gene_symbol",
    zscore_axis: int = 0,
) -> pd.Series:
    """
    Compute regional disease-gene expression scores.

    Score for region r = mean of z-scored expression values
    across the selected disease genes.

    Parameters
    ----------
    expression_matrix : pd.DataFrame
        Region x gene expression matrix.
    disease_genes : pd.DataFrame
        DataFrame containing disease gene identifiers.
    gene_column : str
        Column in disease_genes with gene symbols matching expression_matrix columns.
    zscore_axis : int
        Axis along which to z-score. 0 = across regions (per gene), 1 = across genes (per region).

    Returns
    -------
    pd.Series
        Regional scores indexed by region name.
    """
    # Z-score expression matrix across regions (per gene)
    zscored = zscore_expression(expression_matrix, axis=zscore_axis)

    # Get disease gene symbols that exist in expression matrix
    gene_symbols = disease_genes[gene_column].tolist()
    available_genes = [g for g in gene_symbols if g in zscored.columns]

    if len(available_genes) == 0:
        raise ValueError("No disease genes found in expression matrix")

    # Compute mean z-score across disease genes for each region
    regional_scores = zscored[available_genes].mean(axis=1, skipna=True)
    regional_scores.name = "score"

    # Drop regions with all-NaN expression (no data for any disease gene)
    nan_mask = regional_scores.isna()
    if nan_mask.any():
        dropped = regional_scores.index[nan_mask].tolist()
        print(f"  Dropped {len(dropped)} regions with no expression data: {dropped}")
        regional_scores = regional_scores.dropna()

    return regional_scores

    return regional_scores


def compute_donor_scores(
    donor_matrices: Dict[str, pd.DataFrame],
    disease_genes: pd.DataFrame,
    gene_column: str = "gene_symbol",
) -> pd.DataFrame:
    """
    Compute regional scores per donor for consistency checking.

    Parameters
    ----------
    donor_matrices : dict
        Dictionary mapping donor ID to region x gene matrix.
    disease_genes : pd.DataFrame
        Disease gene list.
    gene_column : str
        Gene symbol column.

    Returns
    -------
    pd.DataFrame
        Regions x donors matrix of scores.
    """
    donor_scores = {}
    for donor_id, matrix in donor_matrices.items():
        scores = compute_regional_scores(matrix, disease_genes, gene_column)
        donor_scores[donor_id] = scores

    return pd.DataFrame(donor_scores)


def compute_consistency_metric(
    donor_scores: pd.DataFrame,
    min_donors: int = 4,
) -> pd.Series:
    """
    Compute donor consistency metric.

    Returns the number of donors where the region's score is positive
    (or above median), indicating consistent directionality.

    Parameters
    ----------
    donor_scores : pd.DataFrame
        Regions x donors score matrix.
    min_donors : int
        Minimum number of donors for "consistent" call.

    Returns
    -------
    pd.Series
        Consistency count per region.
    """
    # Count donors with positive score
    positive_count = (donor_scores > 0).sum(axis=1)
    positive_count.name = "n_positive_donors"
    return positive_count


def save_scores(
    scores: pd.Series,
    disease: str,
    list_type: str,
    donor_scores: Optional[pd.DataFrame] = None,
    consistency: Optional[pd.Series] = None,
) -> None:
    """Save regional scores to results directory."""
    # Main scores
    scores_df = scores.to_frame("score")
    scores_df.index.name = "region"
    scores_df["disease"] = disease
    scores_df["list_type"] = list_type
    scores_file = RESULTS_DIR / f"scores_{disease.lower()}_{list_type}.parquet"
    scores_df.to_parquet(scores_file)

    # Donor-level scores
    if donor_scores is not None:
        donor_file = RESULTS_DIR / f"donor_scores_{disease.lower()}_{list_type}.parquet"
        donor_scores.to_parquet(donor_file)

    # Consistency
    if consistency is not None:
        consistency_df = consistency.to_frame("n_positive_donors")
        consistency_df.index.name = "region"
        consistency_file = RESULTS_DIR / f"consistency_{disease.lower()}_{list_type}.parquet"
        consistency_df.to_parquet(consistency_file)

    print(f"Saved scores for {disease} ({list_type}) to {scores_file}")


def load_scores(disease: str, list_type: str = "strict") -> pd.DataFrame:
    """Load saved regional scores."""
    file = RESULTS_DIR / f"scores_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    else:
        raise FileNotFoundError(f"Scores not found: {file}")