import pandas as pd

# Load expression matrix
expr = pd.read_parquet('data/processed/allen_expression_desikan-killiany.parquet')
expr_genes = set(expr.columns)

# Load Parkinson's strict gene list
pd_strict = pd.read_csv('data/processed/pd_genes_strict.csv')
pd_genes = set(pd_strict['gene'])

# Match genes
matched = pd_genes & expr_genes
missing = pd_genes - expr_genes

print(f"Expression matrix shape: {expr.shape}")
print(f"Number of genes in atlas: {len(expr_genes)}")
print(f"Parkinson's strict genes: {len(pd_genes)}")
print(f"Matched genes: {len(matched)}")
print(f"Missing genes: {len(missing)}")
print(f"Missing: {sorted(missing)}")

# Save match counts for README
match_summary = pd.DataFrame({
    'metric': ['total_atlas_genes', 'total_pd_genes', 'matched', 'missing'],
    'count': [len(expr_genes), len(pd_genes), len(matched), len(missing)]
})
match_summary.to_csv('data/processed/gene_match_counts.csv', index=False)
print("Saved to data/processed/gene_match_counts.csv")
print(match_summary)