"""Gene list builder from a single GWAS Catalog full file.

One parser serves every disease. Trait inclusion / exclusion rules are decided
*before* looking at results and written to ``results/trait_filters.csv`` so the
filter set is auditable and justifiable.

Trait inclusion rules (fixed before running):
  * Include a trait if it contains the disease keyword (case-insensitive).
  * Exclude traits containing any of:
      progression, polygenic, risk score, family history, resilience,
      motor subtype, pleiotropy, atrophy, proxy, MTAG, dyskinesia, adjusted
  * Include age-at-onset traits (still disease-associated).
  * Keep only P-VALUE < 5e-8 (genome-wide significance).

Strict list  = reported-gene mappings only (REPORTED GENE(S)).
Broad list   = reported + mapped genes (REPORTED GENE(S) + MAPPED_GENE).

Outputs:
  * data/processed/gwas_filtered_<disease>.parquet   (trait-filtered associations)
  * data/processed/gwas_genes_<disease>_<list>.parquet (full provenance per gene)
  * data/processed/gene_list_<disease>_<list>.parquet   (unique gene_symbols)
  * results/trait_filters.csv                          (matched traits per disease)
"""

import re
import pandas as pd
from pathlib import Path

RAW_DIR = Path(__file__).parent.parent / "data" / "raw"
PROCESSED_DIR = Path(__file__).parent.parent / "data" / "processed"
RESULTS_DIR = Path(__file__).parent.parent / "results"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

P_THRESHOLD = 5e-8  # genome-wide significance, fixed before looking at results
CHUNKSIZE = 200_000

FULL_CATALOG = RAW_DIR / "gwas_catalog_full.tsv"

# Exclusion keywords (decided before running). A trait containing any of these
# is dropped even if it also contains the disease keyword.
EXCLUDE_KEYWORDS = [
    "progression", "polygenic", "risk score", "family history", "resilience",
    "motor subtype", "pleiotropy", "atrophy", "proxy", "mtag", "dyskinesia",
    "adjusted",
]

DISEASE_KEYWORD = {
    "Parkinson": "Parkinson",
    "Alzheimer": "Alzheimer",
    "Huntington": "Huntington",
    "Height": "Height",
}

STRICT_KEYWORD = "reported"
BROAD_KEYWORDS = ("reported", "mapped")


def _trait_passes(trait: str, disease: str) -> bool:
    """Apply the fixed trait inclusion / exclusion rules."""
    t = str(trait).lower()
    kw = DISEASE_KEYWORD[disease].lower()
    if kw not in t:
        return False
    for excl in EXCLUDE_KEYWORDS:
        if excl in t:
            return False
    return True


def _split_genes(value):
    """Split a GWAS gene cell into a list of cleaned gene symbols."""
    if pd.isna(value):
        return []
    parts = re.split(r"[,;]|\s+-\s+", str(value))
    bad = {"", "intergenic", "nr", "nan", "na", "none"}
    return [p.strip() for p in parts if p.strip().lower() not in bad]


def filter_gwas_by_disease(disease: str, force_refresh: bool = False) -> pd.DataFrame:
    """Read the full catalog in chunks and keep only trait-filtered rows for a disease.

    Returns a DataFrame with columns:
        study, trait, chr_id, chr_pos, p_value, reported_genes, mapped_genes
    """
    cache = PROCESSED_DIR / f"gwas_filtered_{disease.lower()}.parquet"
    if cache.exists() and not force_refresh:
        return pd.read_parquet(cache)

    if not FULL_CATALOG.exists():
        raise FileNotFoundError(f"Full GWAS catalog not found: {FULL_CATALOG}")

    trait_col = None
    # Detect trait column name from the header.
    header = pd.read_csv(FULL_CATALOG, sep="\t", nrows=0)
    for c in header.columns:
        if c.upper().startswith("DISEASE"):
            trait_col = c
            break
    if trait_col is None:
        raise ValueError("No DISEASE/TRAIT column found in full catalog")

    matched_traits = set()
    rows = []

    for ch in pd.read_csv(
        FULL_CATALOG,
        sep="\t",
        chunksize=CHUNKSIZE,
        dtype=str,
        usecols=lambda c: c in {
            "STUDY", trait_col, "CHR_ID", "CHR_POS",
            "REPORTED GENE(S)", "MAPPED_GENE", "P-VALUE",
        },
    ):
        # Filter by trait rules.
        keep = ch[trait_col].apply(lambda t: _trait_passes(t, disease))
        ch = ch[keep]
        if ch.empty:
            continue
        # Record matched traits (for the audit trail).
        for t in ch[trait_col].dropna().unique():
            matched_traits.add((disease, t))
        # Filter by genome-wide significance.
        ch = ch[pd.to_numeric(ch["P-VALUE"], errors="coerce") < P_THRESHOLD]
        if ch.empty:
            continue
        rows.append(
            ch[
                [
                    "STUDY",
                    trait_col,
                    "CHR_ID",
                    "CHR_POS",
                    "REPORTED GENE(S)",
                    "MAPPED_GENE",
                    "P-VALUE",
                ]
            ].copy()
        )

    if rows:
        filtered = pd.concat(rows, ignore_index=True)
    else:
        filtered = pd.DataFrame(
            columns=[
                "STUDY", trait_col, "CHR_ID", "CHR_POS",
                "REPORTED GENE(S)", "MAPPED_GENE", "P-VALUE",
            ]
        )

    # Normalise column names.
    filtered = filtered.rename(
        columns={
            "STUDY": "study",
            trait_col: "trait",
            "CHR_ID": "chr_id",
            "CHR_POS": "chr_pos",
            "REPORTED GENE(S)": "reported_genes",
            "MAPPED_GENE": "mapped_genes",
            "P-VALUE": "p_value",
        }
    )
    filtered["p_value"] = pd.to_numeric(filtered["p_value"], errors="coerce")
    filtered["chr_pos"] = pd.to_numeric(filtered["chr_pos"], errors="coerce")
    filtered["chr_id"] = pd.to_numeric(filtered["chr_id"], errors="coerce")

    filtered.to_parquet(cache, index=False)

    # Write / append the trait filter audit trail.
    _write_trait_filters(matched_traits)

    print(
        f"{disease}: {len(filtered)} trait-filtered associations "
        f"({len(matched_traits)} distinct traits kept)"
    )
    return filtered


def _write_trait_filters(new_traits):
    """Append matched (disease, trait) pairs to results/trait_filters.csv."""
    trail_file = RESULTS_DIR / "trait_filters.csv"
    new_rows = pd.DataFrame(
        [(d, t) for d, t in new_traits], columns=["disease", "trait"]
    )
    if trail_file.exists():
        existing = pd.read_csv(trail_file)
        combined = pd.concat([existing, new_rows], ignore_index=True).drop_duplicates()
    else:
        combined = new_rows
    combined = combined.sort_values(["disease", "trait"]).reset_index(drop=True)
    combined.to_csv(trail_file, index=False)


def _expand_gene_rows(filtered: pd.DataFrame, mapping_types) -> pd.DataFrame:
    """Expand trait-filtered associations into one row per (gene, mapping type)."""
    rows = []
    for _, r in filtered.iterrows():
        for mapping in mapping_types:
            if mapping == "reported":
                genes = _split_genes(r["reported_genes"])
            else:
                genes = _split_genes(r["mapped_genes"])
            for g in genes:
                rows.append(
                    {
                        "gene": g,
                        "mapping_type": mapping,
                        "study": r["study"],
                        "snp": r.get("snp"),
                        "chr_id": r["chr_id"],
                        "chr_pos": r["chr_pos"],
                        "p_value": r["p_value"],
                        "trait": r["trait"],
                    }
                )
    return pd.DataFrame(
        rows,
        columns=[
            "gene", "mapping_type", "study", "snp", "chr_id",
            "chr_pos", "p_value", "trait",
        ],
    )


def build_gene_lists(disease: str, expression_matrix: pd.DataFrame,
                     force_refresh: bool = False):
    """Build strict and broad gene lists for a disease, matched to expression matrix.

    Strict = reported-gene mappings only; Broad = reported + mapped genes.
    Returns (strict_df, broad_df) DataFrames with a ``gene_symbol`` column.
    Provenance tables (with study/SNP/chr/p-value) are also saved.
    """
    cache_strict = PROCESSED_DIR / f"gene_list_{disease.lower()}_strict.parquet"
    cache_broad = PROCESSED_DIR / f"gene_list_{disease.lower()}_broad.parquet"

    if cache_strict.exists() and cache_broad.exists() and not force_refresh:
        strict = pd.read_parquet(cache_strict)
        broad = pd.read_parquet(cache_broad)
        return strict, broad

    # Trait-filter the full catalog for this disease (writes trait_filters.csv).
    filtered = filter_gwas_by_disease(disease, force_refresh=force_refresh)

    expr_genes = set(expression_matrix.columns)

    list_defs = {
        "strict": ({"reported"}, cache_strict,
                   PROCESSED_DIR / f"gwas_genes_{disease.lower()}_strict.parquet"),
        "broad": (set(BROAD_KEYWORDS), cache_broad,
                  PROCESSED_DIR / f"gwas_genes_{disease.lower()}_broad.parquet"),
    }

    result_dfs = {}
    for list_type, (mapping_types, list_cache, prov_cache) in list_defs.items():
        expanded = _expand_gene_rows(filtered, mapping_types)
        if not expanded.empty:
            # Deduplicate by gene, keeping the most significant p-value.
            expanded = expanded.sort_values("p_value")
            expanded = expanded.drop_duplicates(subset=["gene"], keep="first")
            # Provenance table: every gene that passed the filter.
            expanded.to_parquet(prov_cache, index=False)
            # Match to expression matrix columns.
            matched = expanded[expanded["gene"].isin(expr_genes)].copy()
            matched = matched.rename(columns={"gene": "gene_symbol"})
        else:
            matched = pd.DataFrame(
                columns=["gene_symbol", "mapping_type", "study", "snp",
                         "chr_id", "chr_pos", "p_value", "trait"]
            )
            expanded = matched.copy()

        matched.to_parquet(list_cache, index=False)
        result_dfs[list_type] = matched
        print(
            f"{disease} {list_type}: {len(matched)} genes matched to expression "
            f"matrix ({len(expanded)} passed filter before matching)"
        )

    return result_dfs["strict"], result_dfs["broad"]


def build_gene_table(path: str):
    """Standalone helper: trait-filter Parkinson's from a real GWAS file path.

    Kept for backward compatibility with the old __main__ block.
    """
    # This path-based entry is superseded by filter_gwas_by_disease / the full
    # catalog. Retained so the module still imports cleanly.
    raise NotImplementedError(
        "build_gene_table(path) is deprecated; use filter_gwas_by_disease() "
        "with the single full catalog instead."
    )


if __name__ == "__main__":
    # Smoke test: filter each disease and report counts (no expression match).
    for disease in ["Parkinson", "Alzheimer", "Huntington", "Height"]:
        filter_gwas_by_disease(disease, force_refresh=True)
