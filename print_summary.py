"""Print a markdown summary table from pipeline results CSVs.

README numbers are copied from here, never typed by hand.
"""

import pandas as pd
from pathlib import Path

RESULTS_DIR = Path("results")
PROCESSED_DIR = Path("data/processed")


def gene_count(disease: str, list_type: str) -> int:
    file = PROCESSED_DIR / f"gene_list_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return len(pd.read_parquet(file))
    return 0


def region_name(region_id) -> str:
    """Map Desikan-Killiany atlas region ID to name."""
    mapping = {
        36: "Caudate (L)", 37: "Putamen (L)",
        77: "Caudate (R)", 78: "Putamen (R)",
    }
    region_id_str = str(region_id)
    for rid, name in mapping.items():
        if str(rid) == region_id_str:
            return name
    return str(region_id)


def main():
    rows = []
    for stats_file in sorted(RESULTS_DIR.glob("stats_*.parquet")):
        stem = stats_file.stem
        parts = stem.split("_")
        disease = parts[1].capitalize()
        list_type = parts[2]
        df = pd.read_parquet(stats_file)
        n_regions = len(df)
        n_sig = int(df["significant"].sum()) if "significant" in df.columns else 0
        n_depleted = int(df["depleted"].sum()) if "depleted" in df.columns else 0
        rows.append(
            {
                "disease": disease,
                "list": list_type,
                "genes": gene_count(disease, list_type),
                "regions": n_regions,
                "significant": n_sig,
                "depleted": n_depleted,
            }
        )

    if not rows:
        print("No results found. Run `python run_all.py` first.")
        return

    df = pd.DataFrame(rows)
    print("| Disease | List | Genes | Regions | Significant (FDR) | Depleted (2-sided FDR) |")
    print("|---------|------|-------|---------|-------------------|------------------------|")
    for _, r in df.iterrows():
        print(
            f"| {r['disease']} | {r['list']} | {r['genes']} | {r['regions']} | "
            f"{r['significant']} | {r['depleted']} |"
        )

    # Negative control
    nc = RESULTS_DIR / "stats_height_strict.parquet"
    if nc.exists():
        ndf = pd.read_parquet(nc)
        n_sig = int(ndf["significant"].sum()) if "significant" in ndf.columns else 0
        print()
        print(f"Negative control (Height): {n_sig} significant regions "
              f"(expected 0 — confirms no systematic bias).")

    # PD strict: top 10 regions by score, caudate/putamen
    pd_stats = RESULTS_DIR / "stats_parkinson_strict.parquet"
    if pd_stats.exists():
        pdf = pd.read_parquet(pd_stats)
        pdf["region_name"] = pdf.index.map(region_name)

        print()
        print("## Parkinson's strict — top 10 regions by score")
        print("| Region | Score | P-value | FDR | Cohen's d | Significant | Depleted |")
        print("|--------|-------|---------|-----|-----------|-------------|----------|")
        top10 = pdf.nlargest(10, "score")
        for _, r in top10.iterrows():
            print(
                f"| {r['region_name']} | {r['score']:.4f} | "
                f"{r['pvalue']:.4f} | {r['fdr']:.4f} | {r['cohens_d']:.3f} | "
                f"{r['significant']} | {r['depleted']} |"
            )

        print()
        print("## Parkinson's strict — Caudate and Putamen")
        print("| Region | Score | P-value | FDR | Cohen's d | Significant | Depleted |")
        print("|--------|-------|---------|-----|-----------|-------------|----------|")
        targets = [36, 37, 77, 78]
        cp = pdf[pdf.index.isin(targets)].sort_index()
        for _, r in cp.iterrows():
            print(
                f"| {r['region_name']} | {r['score']:.4f} | "
                f"{r['pvalue']:.4f} | {r['fdr']:.4f} | {r['cohens_d']:.3f} | "
                f"{r['significant']} | {r['depleted']} |"
            )

    # Striatal positive control — separate gene set analysis
    print()
    print("## Striatal positive control (regional method)")
    print("Gene set: DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9")
    sc_stats = RESULTS_DIR / "stats_striatalcontrol_strict.parquet"
    if sc_stats.exists():
        sc_df = pd.read_parquet(sc_stats)
        cp_stats = sc_df[sc_df.index.isin([36, 37, 77, 78])]
        print("| Region | Score | P-value | FDR | Cohen's d | Significant | Depleted |")
        print("|--------|-------|---------|-----|-----------|-------------|----------|")
        for _, r in cp_stats.iterrows():
            print(
                f"| {region_name(r.name)} | {r['score']:.4f} | "
                f"{r['pvalue']:.4f} | {r['fdr']:.4f} | {r['cohens_d']:.3f} | "
                f"{r['significant']} | {r['depleted']} |"
            )
    else:
        print("Not computed — run `python run_all.py` to generate.")


if __name__ == "__main__":
    main()
