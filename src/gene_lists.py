"""Gene list construction for NeuroMap: strict and broad versions."""

import pandas as pd
from pathlib import Path
from typing import Tuple

from load_data import load_gwas_catalog, match_genes_to_atlas

PROCESSED_DIR = Path(__file__).parent.parent / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def build_gene_lists(
    disease: str,
    expression_matrix: pd.DataFrame,
    pvalue_threshold_strict: float = 5e-8,
    pvalue_threshold_broad: float = 1e-5,
    force_refresh: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build strict and broad gene lists for a disease.

    Parameters
    ----------
    disease : str
        Disease name.
    expression_matrix : pd.DataFrame
        Expression matrix to match genes against.
    pvalue_threshold_strict : float
        P-value threshold for strict list (genome-wide significant).
    pvalue_threshold_broad : float
        P-value threshold for broad list (suggestive).
    force_refresh : bool
        If True, re-fetch from GWAS Catalog.

    Returns
    -------
    strict_genes : pd.DataFrame
        Genes passing strict threshold, matched to expression matrix.
    broad_genes : pd.DataFrame
        Genes passing broad threshold, matched to expression matrix.
    """
    # Load GWAS data for both thresholds
    strict_gwas = load_gwas_catalog(
        disease, pvalue_threshold=pvalue_threshold_strict, cache=not force_refresh
    )
    broad_gwas = load_gwas_catalog(
        disease, pvalue_threshold=pvalue_threshold_broad, cache=not force_refresh
    )

    # Match to expression matrix
    strict_matched, strict_missing = match_genes_to_atlas(strict_gwas, expression_matrix)
    broad_matched, broad_missing = match_genes_to_atlas(broad_gwas, expression_matrix)

    # Add list type column
    strict_matched["list_type"] = "strict"
    broad_matched["list_type"] = "broad"

    # Save
    strict_file = PROCESSED_DIR / f"gene_list_{disease.lower()}_strict.parquet"
    broad_file = PROCESSED_DIR / f"gene_list_{disease.lower()}_broad.parquet"
    strict_matched.to_parquet(strict_file)
    broad_matched.to_parquet(broad_file)

    print(f"{disease}: Strict={len(strict_matched)} genes, Broad={len(broad_matched)} genes")

    return strict_matched, broad_matched


def load_gene_list(disease: str, list_type: str = "strict") -> pd.DataFrame:
    """Load saved gene list."""
    file = PROCESSED_DIR / f"gene_list_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    else:
        raise FileNotFoundError(f"Gene list not found: {file}")


def get_disease_genes(
    disease: str,
    expression_matrix: pd.DataFrame,
    list_type: str = "strict",
    force_refresh: bool = False,
) -> pd.DataFrame:
    """Convenience function to get disease gene list."""
    if list_type == "strict":
        strict, _ = build_gene_lists(disease, expression_matrix, force_refresh=force_refresh)
        return strict
    elif list_type == "broad":
        _, broad = build_gene_lists(disease, expression_matrix, force_refresh=force_refresh)
        return broad
    else:
        raise ValueError(f"Unknown list_type: {list_type}")