"""Cell-type enrichment analysis for NeuroMap.

Tests whether disease GWAS genes are enriched for specific brain cell-type
marker genes using Fisher's exact test, with FDR correction across all
disease x cell-type tests.

Positive control: microglia markers should be enriched among Alzheimer's
genes (known from GWAS literature).

Diseases with fewer than 10 genes are reported as "not analyzable"
(underpowered) and skipped from the enrichment test.
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Resolve paths relative to repo root (this file lives in src/)
REPO_ROOT = Path(__file__).parent.parent
RESULTS_DIR = REPO_ROOT / "results"
REFERENCE_DIR = REPO_ROOT / "reference"
PROCESSED_DIR = REPO_ROOT / "data" / "processed"

sys.path.insert(0, str(REPO_ROOT / "src"))
from load_data import load_allen_expression

MIN_GENES_FOR_ANALYSIS = 10


def load_marker_genes(path="cell_markers.csv"):
    """Load cell-type marker genes."""
    df = pd.read_csv(REFERENCE_DIR / path)
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
            f = PROCESSED_DIR / f"gene_list_{disease.lower()}_{lt}.parquet"
            if f.exists():
                df = pd.read_parquet(f)
                key = f"{disease}_{lt}"
                n_genes = len(df)
                disease_genes[key] = df["gene_symbol"].tolist()
                print(f"{key}: {n_genes} genes")
            else:
                print(f"{disease}_{lt}: not found, skipping")

    # Run Fisher's exact test for all disease x cell_type combinations
    results = []
    for dkey, dgenes in disease_genes.items():
        disease = dkey.rsplit("_", 1)[0]
        list_type = dkey.rsplit("_", 1)[1]

        # Skip diseases with too few genes (underpowered)
        if len(dgenes) < MIN_GENES_FOR_ANALYSIS:
            print(f"\n{dkey}: {len(dgenes)} genes < {MIN_GENES_FOR_ANALYSIS} — "
                  f"not analyzable, skipping enrichment test.")
            for ctype in markers:
                results.append({
                    "disease": disease,
                    "list_type": list_type,
                    "cell_type": ctype,
                    "n_disease_genes": len(dgenes),
                    "n_marker_genes": len(markers[ctype]),
                    "n_background_genes": len(background_genes),
                    "n_overlap": 0,
                    "odds_ratio": float("nan"),
                    "p_value": float("nan"),
                    "fdr": float("nan"),
                    "significant": False,
                    "analyzable": False,
                })
            continue

        d_set = set(dgenes) & background_genes

        for ctype, mgenes in markers.items():
            m_set = mgenes & background_genes
            result = fisher_exact_test(d_set, m_set, background_genes)
            result["disease"] = disease
            result["list_type"] = list_type
            result["cell_type"] = ctype
            result["analyzable"] = True
            results.append(result)

    results_df = pd.DataFrame(results)

    # FDR correction across all analyzable tests
    valid = results_df[results_df["analyzable"] & results_df["p_value"].notna()]
    if len(valid) > 0:
        reject, fdr_corrected, _, _ = multipletests(
            valid["p_value"].values, alpha=0.05, method="fdr_bh"
        )
        results_df.loc[valid.index, "fdr"] = fdr_corrected
        results_df.loc[valid.index, "significant"] = reject

    # Sort by FDR (analyzable first, then by FDR)
    results_df = results_df.sort_values(
        ["analyzable", "fdr"], ascending=[False, True],
        na_position="last"
    )

    # Save results
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(RESULTS_DIR / "cell_type_enrichment.csv", index=False)

    # Print results
    print("\n" + "=" * 80)
    print("CELL-TYPE ENRICHMENT RESULTS (Fisher's exact test, FDR-corrected)")
    print("=" * 80)
    print()
    print(f"{'Disease':<15} {'List':<8} {'Cell Type':<25} {'OR':>8} {'p':>10} {'FDR':>10} {'Sig':>6} {'Overlap':>8} {'Analyzable':>10}")
    print("-" * 100)
    for _, r in results_df.iterrows():
        sig = "**" if r.get("significant", False) else ""
        ana = "Yes" if r.get("analyzable", False) else "No"
        print(f"{r['disease']:<15} {r['list_type']:<8} {r['cell_type']:<25} "
              f"{r['odds_ratio']:>8.2f} {r['p_value']:>10.4g} {r.get('fdr', float('nan')):>10.4g} "
              f"{sig:>6} {r['n_overlap']:>8} {ana:>10}")

    # Positive control check: microglia in Alzheimer's
    print()
    print("POSITIVE CONTROL: Microglia in Alzheimer's")
    alz_microglia = results_df[
        (results_df["disease"] == "Alzheimer") &
        (results_df["cell_type"] == "Microglia") &
        (results_df["analyzable"] == True)
    ]
    for _, r in alz_microglia.iterrows():
        print(f"  {r['list_type']}: OR={r['odds_ratio']:.2f}, p={r['p_value']:.4g}, "
              f"FDR={r.get('fdr', float('nan')):.4g}, sig={r.get('significant', False)}")


if __name__ == "__main__":
    main()
