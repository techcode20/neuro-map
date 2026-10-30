"""Null-model calibration test for NeuroMap statistics engine.

Independent calibration: does NOT use stats.py. Acts as a reference implementation
to verify that the null model and FDR correction behave correctly on random data.

Healthy result: ~5 or fewer runs out of 100 have any FDR-significant region,
and the p-value histogram is roughly flat.
"""

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def run_once(expr, n_genes=155, n_null=500, rng=None):
    """Run one null calibration: random genes, matched null, FDR correction."""
    if rng is None:
        rng = np.random.default_rng(42)

    z = ((expr - expr.mean()) / expr.std(ddof=0)).values  # regions x genes
    n_genes_total = z.shape[1]

    # Bin ALL genes by mean expression (single shared call)
    mean_per_gene = expr.mean().values
    bins = pd.qcut(mean_per_gene, 10, labels=False, duplicates="drop")
    bins_arr = np.asarray(bins).astype(int)

    # Pick random "disease" genes
    idx = rng.choice(n_genes_total, n_genes, replace=False)
    obs = z[:, idx].mean(axis=1)

    # Build null pool (exclude disease genes)
    pool = np.setdiff1d(np.arange(n_genes_total), idx)
    by_bin = {b: pool[bins_arr[pool] == b] for b in np.unique(bins_arr)}
    counts = pd.Series(bins_arr[idx]).value_counts()

    null = np.empty((n_null, z.shape[0]))
    for k in range(n_null):
        pick = np.concatenate([
            rng.choice(by_bin[b], c, replace=False)
            for b, c in counts.items()
            if b in by_bin and len(by_bin[b]) >= c
        ])
        null[k] = z[:, pick].mean(axis=1)

    pvals = (1 + (null >= obs).sum(axis=0)) / (1 + n_null)
    return pvals


def main():
    expr = pd.read_parquet("data/processed/allen_expression_desikan-killiany.parquet")
    expr = expr.dropna(axis=0, how="all").dropna(axis=1)  # drop empty regions / NaN genes
    print(f"Expression matrix: {expr.shape[0]} regions x {expr.shape[1]} genes")

    rng = np.random.default_rng(42)
    hits = 0
    allp = []
    n_runs = 100

    for _ in range(n_runs):
        p = run_once(expr, n_genes=155, n_null=500, rng=rng)
        allp.append(p)
        sig_any = multipletests(p, method="fdr_bh")[0].any()
        hits += sig_any

    hist = np.histogram(np.concatenate(allp), bins=10, range=(0, 1))[0]

    print(f"Runs with any FDR-significant region: {hits} / {n_runs}")
    print(f"P-value histogram (10 bins): {hist.tolist()}")
    print()
    if hits <= 5:
        print("CALIBRATION PASSED: null model and FDR behave correctly.")
    else:
        print("CALIBRATION FAILED: too many false positives. Null model is broken.")


if __name__ == "__main__":
    main()