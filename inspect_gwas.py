import pandas as pd

path = "data/raw/gwas_parkinsons.tsv"
cols = pd.read_csv(path, sep="\t", nrows=0).columns
print("Columns:", list(cols))
trait = next(c for c in cols if c.upper().startswith("DISEASE"))

counts, total = {}, 0
for ch in pd.read_csv(path, sep="\t", usecols=[trait], chunksize=200_000, dtype=str):
    total += len(ch)
    s = ch[trait].fillna("")
    for key in ["arkinson", "lzheimer", "untington"]:
        for v, n in s[s.str.contains(key, case=False)].value_counts().items():
            counts[(key, v)] = counts.get((key, v), 0) + n

print("Total rows:", total)
for key in ["arkinson", "lzheimer", "untington"]:
    rows = sorted([(n, v) for (k, v), n in counts.items() if k == key], reverse=True)[:10]
    print(f"\n{key}: {sum(n for n, _ in rows)} rows in top 10 traits")
    for n, v in rows:
        print(f"  {n:5d}  {v}")
