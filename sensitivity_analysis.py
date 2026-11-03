"""Sensitivity analyses for NeuroMap.

Recomputed after the gene-list rebuild (Step 3 of the audit).
All p-values are unadjusted, exploratory.

Outputs:
  * results/overlap_loci.csv — overlap genes with locus assignments
  * Prints leave-one-locus-out table and OR without MHC locus
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

RESULTS_DIR = Path("results")
PROCESSED_DIR = Path("data/processed")
REFERENCE_DIR = Path("reference")

P_THRESHOLD = 5e-8
LOCUS_WINDOW = 500_000  # 500kb
MHC_CHR = 6
MHC_START = 25_000_000
MHC_END = 35_000_000


def load_microglia_markers(path="reference/cell_markers.csv"):
    df = pd.read_csv(path)
    return set(df[df["cell_type"] == "Microglia"]["gene_symbol"].tolist())


def load_strict_genes(disease):
    f = PROCESSED_DIR / f"gene_list_{disease.lower()}_strict.parquet"
    if not f.exists():
        return set()
    return set(pd.read_parquet(f)["gene_symbol"].tolist())


def load_provenance(disease, list_type="strict"):
    f = PROCESSED_DIR / f"gwas_genes_{disease.lower()}_{list_type}.parquet"
    if not f.exists():
        return pd.DataFrame()
    return pd.read_parquet(f)


def assign_loci(df):
    """Assign locus IDs based on 500kb windows around lead-SNP positions.

    Genes on the same chromosome within LOCUS_WINDOW bp of each other
    are grouped into one locus. Genes with no position info get their
    own 'no_position' locus.
    """
    df = df.copy()
    df["locus"] = "no_position"

    has_pos = df["chr_id"].notna() & df["chr_pos"].notna()
    if not has_pos.any():
        return df

    sub = df[has_pos].copy()
    sub = sub.sort_values(["chr_id", "chr_pos"])

    locus_id = 0
    current_chr = None
    current_end = None

    for idx, row in sub.iterrows():
        cid = row["chr_id"]
        cpos = row["chr_pos"]
        if current_chr != cid or (cpos - current_end) > LOCUS_WINDOW:
            locus_id += 1
            current_chr = cid
            current_end = cpos
        else:
            current_end = max(current_end, cpos)
        df.at[idx, "locus"] = f"chr{cid}_locus{locus_id}"

    return df


def is_mhc(row):
    """True if the gene's position falls in the MHC region (chr6:25-35Mb)."""
    if pd.isna(row["chr_id"]) or pd.isna(row["chr_pos"]):
        return False
    return int(row["chr_id"]) == MHC_CHR and MHC_START <= row["chr_pos"] <= MHC_END


def fisher_or(a_set, b_set, background):
    """Fisher's exact test odds ratio and p-value for overlap of two sets."""
    d_set = set(a_set) & background
    m_set = set(b_set) & background
    n_dm = len(d_set & m_set)
    n_d = len(d_set)
    n_m = len(m_set)
    n_bg = len(background)
    table = [
        [n_dm, n_d - n_dm],
        [n_m - n_dm, n_bg - n_d - n_m + n_dm],
    ]
    or_, p = stats.fisher_exact(table, alternative="greater")
    return or_, p, n_dm


def main():
    microglia = load_microglia_markers()
    background = set()
    expr = pd.read_parquet(PROCESSED_DIR / "allen_expression_desikan-killiany.parquet")
    background = set(expr.columns)

    print("=== Overlap genes between each strict list and microglia markers ===\n")

    all_overlap_rows = []

    for disease in ["Parkinson", "Alzheimer", "Huntington"]:
        dgenes = load_strict_genes(disease)
        overlap = dgenes & microglia
        prov = load_provenance(disease, "strict")
        overlap_prov = prov[prov["gene"].isin(overlap)].copy()

        # Add locus assignments and MHC flag
        overlap_prov = assign_loci(overlap_prov)
        overlap_prov["in_mhc"] = overlap_prov.apply(is_mhc, axis=1)
        overlap_prov["disease"] = disease

        print(f"--- {disease} strict ---")
        print(f"Disease genes: {len(dgenes)}")
        print(f"Microglia overlap genes ({len(overlap)}): {sorted(overlap)}")
        print(f"In MHC locus (chr6:25-35Mb): {overlap_prov['in_mhc'].sum()}")
        print()

        for _, r in overlap_prov.iterrows():
            all_overlap_rows.append(
                {
                    "disease": disease,
                    "gene": r["gene"],
                    "chr_id": r["chr_id"],
                    "chr_pos": r["chr_pos"],
                    "locus": r["locus"],
                    "in_mhc": bool(r["in_mhc"]),
                    "p_value": r["p_value"],
                    "study": r["study"],
                }
            )

    # Save overlap loci
    loci_df = pd.DataFrame(all_overlap_rows)
    loci_df.to_csv(RESULTS_DIR / "overlap_loci.csv", index=False)
    print(f"Saved {len(loci_df)} overlap gene rows to results/overlap_loci.csv\n")

    # Sensitivity: OR without MHC locus, OR with largest locus removed, LOO
    print("=== Sensitivity analyses (unadjusted, exploratory) ===\n")

    for disease in ["Parkinson", "Alzheimer"]:
        dgenes = load_strict_genes(disease)
        overlap = sorted(dgenes & microglia)

        # Full OR
        or_full, p_full, n_full = fisher_or(dgenes, microglia, background)
        print(f"--- {disease} strict × Microglia ---")
        print(f"Full: OR={or_full:.2f}, p={p_full:.4g}, overlap={n_full} genes: {overlap}")

        # OR without MHC locus
        overlap_prov = loci_df[loci_df["disease"] == disease]
        non_mhc_overlap = overlap_prov[~overlap_prov["in_mhc"]]["gene"].tolist()
        dgenes_no_mhc = set(g for g in dgenes if g not in set(overlap_prov[overlap_prov["in_mhc"]]["gene"]))
        or_nomhc, p_nomhc, n_nomhc = fisher_or(dgenes_no_mhc, microglia, background)
        print(
            f"Without MHC locus: OR={or_nomhc:.2f}, p={p_nomhc:.4g}, "
            f"overlap={n_nomhc} genes: {sorted(set(overlap) & set(dgenes_no_mhc))}"
        )

        # OR with largest locus removed
        locus_counts = overlap_prov["locus"].value_counts()
        if not locus_counts.empty:
            largest_locus = locus_counts.index[0]
            largest_genes = set(
                overlap_prov[overlap_prov["locus"] == largest_locus]["gene"]
            )
            dgenes_no_largest = set(g for g in dgenes if g not in largest_genes)
            or_nolarge, p_nolarge, n_nolarge = fisher_or(
                dgenes_no_largest, microglia, background
            )
            print(
                f"Without largest locus ({largest_locus}): OR={or_nolarge:.2f}, "
                f"p={p_nolarge:.4g}, overlap={n_nolarge} genes: "
                f"{sorted(set(overlap) & dgenes_no_largest)}"
            )

        # Leave-one-locus-out
        print("Leave-one-locus-out:")
        for locus in sorted(overlap_prov["locus"].unique()):
            locus_genes = set(
                overlap_prov[overlap_prov["locus"] == locus]["gene"]
            )
            dgenes_loo = set(g for g in dgenes if g not in locus_genes)
            or_loo, p_loo, n_loo = fisher_or(dgenes_loo, microglia, background)
            remaining = sorted(set(overlap) & dgenes_loo)
            print(
                f"  Without {locus}: OR={or_loo:.2f}, p={p_loo:.4g}, "
                f"overlap={n_loo} genes: {remaining}"
            )
        print()


if __name__ == "__main__":
    main()
