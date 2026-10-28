"""Gene list builder from GWAS Catalog TSV download."""

import re
import pandas as pd

P_THRESHOLD = 5e-8  # fixed before looking at results

def split_genes(value):
    if pd.isna(value):
        return []
    parts = re.split(r"[,;]|\s+-\s+", str(value))
    bad = {"", "intergenic", "nr", "nan"}
    return [p.strip() for p in parts if p.strip().lower() not in bad]

def build_gene_table(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df[pd.to_numeric(df["P-VALUE"], errors="coerce") < P_THRESHOLD]

    rows = []
    for _, r in df.iterrows():
        for g in split_genes(r["REPORTED GENE(S)"]):
            rows.append((g, "reported", r["STUDY ACCESSION"], r["SNPS"], r["P-VALUE"]))
        for g in split_genes(r["MAPPED_GENE"]):
            rows.append((g, "mapped", r["STUDY ACCESSION"], r["SNPS"], r["P-VALUE"]))
    return pd.DataFrame(rows, columns=["gene", "mapping_type", "study", "snp", "p_value"])

if __name__ == "__main__":
    table = build_gene_table("data/raw/gwas_parkinsons.tsv")
    table.to_csv("data/processed/pd_genes_all.csv", index=False)

    strict = table[table.mapping_type == "reported"].gene.drop_duplicates()
    broad = table.gene.drop_duplicates()
    strict.to_csv("data/processed/pd_genes_strict.csv", index=False)
    broad.to_csv("data/processed/pd_genes_broad.csv", index=False)
    print("Strict genes:", len(strict), "| Broad genes:", len(broad))