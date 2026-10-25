"""Data loading and preprocessing utilities for NeuroMap."""

import pandas as pd
import numpy as np
from pathlib import Path
import abagen
from typing import Optional, List, Tuple


DATA_DIR = Path(__file__).parent.parent / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def load_allen_expression(
    atlas: str = "desikan-killiany",
    donors: Optional[List[str]] = None,
    force_refresh: bool = False,
) -> Tuple[pd.DataFrame, dict]:
    """
    Load and process Allen Human Brain Atlas expression data using abagen.

    Parameters
    ----------
    atlas : str
        Atlas to use for parcellation. Options: 'desikan-killiany', 'schaefer-100', 'schaefer-200', etc.
    donors : list of str, optional
        Specific donor IDs to include. If None, uses all available donors.
    force_refresh : bool
        If True, re-download and reprocess data even if cached.

    Returns
    -------
    expression_matrix : pd.DataFrame
        Region x gene expression matrix (averaged across donors if multiple).
    donor_matrices : dict
        Dictionary mapping donor ID to region x gene matrix for that donor.
    """
    cache_file = PROCESSED_DIR / f"allen_expression_{atlas}.parquet"
    donor_cache_dir = PROCESSED_DIR / f"donor_matrices_{atlas}"
    donor_cache_dir.mkdir(parents=True, exist_ok=True)

    if cache_file.exists() and not force_refresh:
        print(f"Loading cached expression matrix from {cache_file}")
        expression_matrix = pd.read_parquet(cache_file)
        donor_matrices = {}
        for donor_file in donor_cache_dir.glob("*.parquet"):
            donor_id = donor_file.stem
            donor_matrices[donor_id] = pd.read_parquet(donor_file)
        return expression_matrix, donor_matrices

    print(f"Fetching Allen data with atlas={atlas}...")
    expression_data = abagen.get_expression_data(atlas=atlas, donors=donors)

    # expression_data is a dict with donor IDs as keys, each containing a DataFrame
    donor_matrices = {}
    for donor_id, df in expression_data.items():
        # abagen returns genes as columns, samples as rows
        # We need to aggregate samples to regions
        # For now, use the mean across samples within each region
        region_expr = df.groupby(level=0, axis=0).mean()
        donor_matrices[donor_id] = region_expr
        donor_file = donor_cache_dir / f"{donor_id}.parquet"
        region_expr.to_parquet(donor_file)

    # Average across donors
    if len(donor_matrices) > 1:
        expression_matrix = pd.concat(donor_matrices.values()).groupby(level=0).mean()
    else:
        expression_matrix = list(donor_matrices.values())[0]

    # Save averaged matrix
    expression_matrix.to_parquet(cache_file)
    print(f"Saved expression matrix to {cache_file}")
    print(f"Shape: {expression_matrix.shape} (regions x genes)")

    return expression_matrix, donor_matrices


def load_gwas_catalog(
    disease: str,
    pvalue_threshold: float = 5e-8,
    cache: bool = True,
) -> pd.DataFrame:
    """
    Load disease-associated genes from GWAS Catalog.

    Parameters
    ----------
    disease : str
        Disease name (e.g., 'Parkinson', 'Alzheimer', 'Huntington').
    pvalue_threshold : float
        Genome-wide significance threshold.
    cache : bool
        If True, use cached version if available.

    Returns
    -------
    pd.DataFrame
        Columns: gene_symbol, ensembl_id, reported_gene, mapped_gene,
        pvalue, study_accession, trait, mapping_method.
    """
    cache_file = RAW_DIR / f"gwas_catalog_{disease.lower()}.tsv"

    if cache_file.exists() and cache:
        print(f"Loading cached GWAS Catalog for {disease}")
        return pd.read_csv(cache_file, sep="\t")

    # GWAS Catalog REST API endpoint
    url = "https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms"
    params = {
        "trait": disease,
        "pvalueThreshold": pvalue_threshold,
        "size": 10000,
    }

    print(f"Fetching GWAS Catalog data for {disease}...")
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()

    # Parse response (GWAS Catalog returns JSON)
    data = response.json()

    # Extract relevant fields
    records = []
    for snp in data.get("_embedded", {}).get("singleNucleotidePolymorphisms", []):
        for gene in snp.get("genes", []):
            record = {
                "gene_symbol": gene.get("geneName"),
                "ensembl_id": gene.get("ensemblId"),
                "reported_gene": gene.get("reportedGene"),
                "mapped_gene": gene.get("mappedGene"),
                "pvalue": snp.get("pvalue"),
                "study_accession": snp.get("study", {}).get("accessionId"),
                "trait": disease,
                "mapping_method": gene.get("mappingMethod"),
            }
            records.append(record)

    df = pd.DataFrame(records)
    df = df.dropna(subset=["gene_symbol"])
    df = df.drop_duplicates(subset=["gene_symbol"])

    # Save cache
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_file, sep="\t", index=False)
    print(f"Saved GWAS Catalog data to {cache_file} ({len(df)} genes)")

    return df


def match_genes_to_atlas(
    gene_list: pd.DataFrame,
    expression_matrix: pd.DataFrame,
    id_column: str = "gene_symbol",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Match gene list to expression matrix columns.

    Parameters
    ----------
    gene_list : pd.DataFrame
        DataFrame with gene identifiers.
    expression_matrix : pd.DataFrame
        Expression matrix with genes as columns.
    id_column : str
        Column in gene_list to match against expression matrix columns.

    Returns
    -------
    matched_genes : pd.DataFrame
        Subset of gene_list with genes found in expression matrix.
    missing_genes : pd.DataFrame
        Subset of gene_list with genes NOT found in expression matrix.
    """
    expr_genes = set(expression_matrix.columns)
    gene_symbols = set(gene_list[id_column])

    matched = gene_list[gene_list[id_column].isin(expr_genes)].copy()
    missing = gene_list[~gene_list[id_column].isin(expr_genes)].copy()

    print(f"Matched {len(matched)}/{len(gene_list)} genes to expression matrix")
    if len(missing) > 0:
        print(f"Missing genes: {missing[id_column].tolist()}")

    return matched, missing


def zscore_expression(expression_matrix: pd.DataFrame, axis: int = 0) -> pd.DataFrame:
    """
    Z-score expression matrix across genes (axis=0) or regions (axis=1).

    For regional scoring, we z-score across genes (axis=0) so each gene
    has mean=0, std=1 across regions, then average across disease genes.
    """
    return expression_matrix.apply(
        lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0,
        axis=axis
    )


# Import requests here to avoid circular imports
import requests