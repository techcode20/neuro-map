"""Gene list builder from GWAS Catalog TSV download."""

import re
import pandas as pd
from pathlib import Path

PROCESSED_DIR = Path(__file__).parent.parent / "data" / "processed"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

P_THRESHOLD = 5e-8  # fixed before looking at results
DISEASE_TRAIT_MAP = {
    "Parkinson": "Parkinson",
    "Alzheimer": "Alzheimer",
    "Huntington": "Huntington",
}

def split_genes(value):
    if pd.isna(value):
        return []
    parts = re.split(r"[,;]|\s+-\s+", str(value))
    bad = {"", "intergenic", "nr", "nan"}
    return [p.strip() for p in parts if p.strip().lower() not in bad]

def _load_gwas_table(disease):
    """Load and filter GWAS catalog for a disease."""
    path = "data/raw/gwas_parkinsons.tsv"
    df = pd.read_csv(path, sep="\t", low_memory=False)
    
    # Filter for disease
    trait = DISEASE_TRAIT_MAP.get(disease, disease)
    df = df[df["DISEASE/TRAIT"].str.contains(trait, case=False, na=False)]
    
    df = df[pd.to_numeric(df["P-VALUE"], errors="coerce") < P_THRESHOLD]
    return df

def _build_gene_table(df):
    """Build gene table from filtered GWAS dataframe."""
    study_col = "STUDY ACCESSION" if "STUDY ACCESSION" in df.columns else "STUDY"
    
    rows = []
    for _, r in df.iterrows():
        for g in split_genes(r["REPORTED GENE(S)"]):
            rows.append((g, "reported", r[study_col], r["SNPS"], r["P-VALUE"]))
        for g in split_genes(r["MAPPED_GENE"]):
            rows.append((g, "mapped", r[study_col], r["SNPS"], r["P-VALUE"]))
    return pd.DataFrame(rows, columns=["gene", "mapping_type", "study", "snp", "p_value"])

def _match_genes_to_expression(genes, expression_matrix):
    """Match gene list to expression matrix columns."""
    expr_genes = set(expression_matrix.columns)
    matched = [g for g in genes if g in expr_genes]
    return matched

def build_gene_lists(disease, expression_matrix, force_refresh=False):
    """
    Build strict and broad gene lists for a disease, matched to expression matrix.
    
    Expected API for run_all.py pipeline.
    """
    cache_strict = PROCESSED_DIR / f"gene_list_{disease.lower()}_strict.parquet"
    cache_broad = PROCESSED_DIR / f"gene_list_{disease.lower()}_broad.parquet"
    
    if cache_strict.exists() and cache_broad.exists() and not force_refresh:
        strict = pd.read_parquet(cache_strict)
        broad = pd.read_parquet(cache_broad)
        return strict, broad
    
    # Load and filter GWAS data
    df = _load_gwas_table(disease)
    table = _build_gene_table(df)
    
    # Strict = reported genes, Broad = all genes
    strict_genes = table[table.mapping_type == "reported"].gene.drop_duplicates().tolist()
    broad_genes = table.gene.drop_duplicates().tolist()
    
# Match to expression matrix
    strict_matched = _match_genes_to_expression(strict_genes, expression_matrix)
    broad_matched = _match_genes_to_expression(broad_genes, expression_matrix)

    # Save as DataFrames with gene_symbol column (matching scoring.py API)
    strict_df = pd.DataFrame({"gene_symbol": strict_matched})
    broad_df = pd.DataFrame({"gene_symbol": broad_matched})
    
    strict_df.to_parquet(cache_strict)
    broad_df.to_parquet(cache_broad)
    
    print(f"{disease}: Strict={len(strict_matched)} genes, Broad={len(broad_matched)} genes")
    
    return strict_df, broad_df

def build_gene_table(path):
    """Standalone function for direct GWAS file processing (Parkinson's only)."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[df["DISEASE/TRAIT"].str.contains("Parkinson", case=False, na=False)]
    df = df[pd.to_numeric(df["P-VALUE"], errors="coerce") < P_THRESHOLD]
    return _build_gene_table(df)

if __name__ == "__main__":
    table = build_gene_table("data/raw/gwas_parkinsons.tsv")
    table.to_csv("data/processed/pd_genes_all.csv", index=False)
    
    strict = table[table.mapping_type == "reported"].gene.drop_duplicates()
    broad = table.gene.drop_duplicates()
    strict.to_csv("data/processed/pd_genes_strict.csv", index=False)
    broad.to_csv("data/processed/pd_genes_broad.csv", index=False)
    print("Strict genes:", len(strict), "| Broad genes:", len(broad))