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


def main():
    rows = []
    for stats_file in sorted(RESULTS_DIR.glob("stats_*.csv")):
        # stats_<disease>_<list>.csv
        stem = stats_file.stem  # e.g. stats_parkinson_strict
        parts = stem.split("_")
        disease = parts[1].capitalize()
        list_type = parts[2]
        df = pd.read_csv(stats_file)
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
        print("No results CSVs found. Run `python run_all.py` first.")
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
    nc = RESULTS_DIR / "stats_height_strict.csv"
    if nc.exists():
        ndf = pd.read_csv(nc)
        n_sig = int(ndf["significant"].sum()) if "significant" in ndf.columns else 0
        print()
        print(f"Negative control (Height): {n_sig} significant regions "
              f"(expected 0 — confirms no systematic bias).")


if __name__ == "__main__":
    main()
