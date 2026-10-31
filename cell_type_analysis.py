"""Cell-type enrichment analysis for NeuroMap.

Tests whether disease GWAS genes are enriched for specific brain cell-type
marker genes using Fisher's exact test, with FDR correction across all
disease x cell-type tests.

Positive control: microglia markers should be enriched among Alzheimer's
genes (known from GWAS literature).
"""

import pandas as pd
import numpy as np
import sys
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path

sys.path.insert(0, "src")
from load_data import load_allen_expression

RESULTS_DIR = Path("results")
REFERENCE_DIR = Path("reference")


def load_marker_genes(path="reference/cell_markers.csv"):
    """Load cell-type marker genes."""
    df = pd.read_csv(path)
    markers = {}
    for cell_type in df["cell_type"].unique():
        genes = df[df["cell_type"] == cell_type]["gene_symbol"].tolist()
        markers[cell_type] = set(genes)
    return markers


def fisher_exact_test(disease_genes, marker_genes, background_genes):
    """Fisher's exact test for enrichment of marker genes among disease genes."""
    d_set = set(disease_genes)
    m_set = set(marker_genes)
    bg_set = set(background_genes)

    # Count: disease genes that are markers vs not markers
    n_disease_marker = len(d_set & m_set)
    n_disease_nonmarker = len(d_set - m_set)
    n_nondisease_marker = len((bg_set - d_set) & m_set)
    n_nondisease_nonmarker = len((bg_set - d_set) - m_set)

    table = [[n_disease_marker, n_disease_nonmarker],
             [n_nondisease_marker, n_nondisease_nonmarker]]

    odds_ratio, p_value = stats.fisher_exact(table, alternative="greater")

    return {
        "n_disease_genes": len(d_set),
        "n_marker_genes": len(m_set),
        "n_background_genes": len(bg_set),
        "n_overlap": n_disease_marker,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
    }


def main():
    # Load expression matrix (real data)
    expr, donors = load_allen_expression(atlas_name="desikan-killiany", force_refresh=False)
    background_genes = set(expr.columns)
    print(f"Background genes: {len(background_genes)}")

    # Load marker genes
    markers = load_marker_genes()
    print(f"Cell types: {list(markers.keys())}")
    for ct, genes in markers.items():
        print(f"  {ct}: {len(genes)} genes")

    # Load disease gene lists
    diseases = ["Parkinson", "Alzheimer", "Huntington"]
    disease_genes = {}
    for disease in diseases:
        for lt in ["strict", "broad"]:
            try:
                df = pd.read_parquet(f"data/processed/gene_list_{disease.lower()}_{lt}.parquet")
                key = f"{disease}_{lt}"
                disease_genes[key] = df["gene_symbol"].tolist()
                print(f"{key}: {len(disease_genes[key])} genes")
            except FileNotFoundError:
                print(f"{disease}_{lt}: not found, skipping")

    # Run Fisher's exact test for all disease x cell_type combinations
    results = []
    for dkey, dgenes in disease_genes.items():
        disease = dkey.rsplit("_", 1)[0]  # Strip _strict/_broad
        list_type = dkey.rsplit("_", 1)[1]
        d_set = set(dgenes) & background_genes

        for ctype, mgenes in markers.items():
            m_set = mgenes & background_genes
            result = fisher_exact_test(d_set, m_set, background_genes)
            result["disease"] = disease
            result["list_type"] = list_type
            result["cell_type"] = ctype
            results.append(result)

    results_df = pd.DataFrame(results)

    # FDR correction across all tests
    valid = results_df.dropna(subset=["p_value"])
    reject, fdr_corrected, _, _ = multipletests(
        valid["p_value"].values, alpha=0.05, method="fdr_bh"
    )
    results_df.loc[valid.index, "fdr"] = fdr_corrected
    results_df.loc[valid.index, "significant"] = reject

    # Sort by FDR
    results_df = results_df.sort_values("fdr")

    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(RESULTS_DIR / "cell_type_enrichment.csv", index=False)

    # Print results
    print("\n" + "=" * 80)
    print("CELL-TYPE ENRICHMENT RESULTS (Fisher's exact test, FDR-corrected)")
    print("=" * 80)
    print()
    print(f"{'Disease':<15} {'List':<8} {'Cell Type':<25} {'OR':>8} {'p':>10} {'FDR':>10} {'Sig':>6} {'Overlap':>8}")
    print("-" * 80)
    for _, r in results_df.iterrows():
        sig = "**" if r.get("significant", False) else ""
        print(f"{r['disease']:<15} {r['list_type']:<8} {r['cell_type']:<25} "
              f"{r['odds_ratio']:>8.2f} {r['p_value']:>10.4f} {r.get('fdr', float('nan')):>10.4f} "
              f"{sig:>6} {r['n_overlap']:>8}")

    # Positive control check: microglia in Alzheimer's
    print()
    print("POSITIVE CONTROL: Microglia in Alzheimer's")
    alz_microglia = results_df[
        (results_df["disease"] == "Alzheimer") & (results_df["cell_type"] == "Microglia")
    ]
    for _, r in alz_microglia.iterrows():
        print(f"  {r['list_type']}: OR={r['odds_ratio']:.2f}, p={r['p_value']:.4f}, FDR={r.get('fdr', float('nan')):.4f}, sig={r.get('significant', False)}")


if __name__ == "__main__":
    main()