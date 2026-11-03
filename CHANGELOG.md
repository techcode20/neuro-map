# Changelog

## v1.1.0 (unreleased) — Data-pipeline audit

### Fixed
- `src/gene_lists.py`: was hardcoded to read `data/raw/gwas_parkinsons.tsv` for every disease, so Alzheimer's and Huntington's could not produce gene lists. Now reads the single full catalog `data/raw/gwas_catalog_full.tsv` in chunks and filters per disease by trait name.
- Mock/real gene-symbol mismatch: in mock mode the pipeline previously matched real GWAS symbols against synthetic `GENE_XXXX` expression genes (0 matches). Now prefers normalized catalog files (`GENE_XXXX`) when the expression matrix is synthetic, and real GWAS Catalog files when it is real.
- `src/scoring.py`: `compute_regional_scores` raised `ValueError` when no disease genes matched; now returns an empty series with a warning so the pipeline does not crash.
- `src/stats.py`: `compute_effect_sizes` crashed on empty results; now returns an empty DataFrame.
- `run_all.py`: skips scoring/stats for empty gene lists with a clear warning instead of crashing.
- `app.py`: removed dead duplicate `main()` (5-tab version shadowed by 7-tab version); removed duplicated "Genes Driving Top Regions" block; fixed `st.sidebar.radio()` called inside a tab (Streamlit error).
- `data/raw/gwas_parkinsons.tsv` renamed to `data/raw/gwas_catalog_full.tsv` (it was the full GWAS catalog, 1.2M rows, misnamed).

### Changed
- Gene lists are now trait-filtered from the full catalog with rules decided before running and written to `results/trait_filters.csv`:
  - include trait containing the disease keyword (case-insensitive)
  - exclude traits containing: progression, polygenic, risk score, family history, resilience, motor subtype, pleiotropy, atrophy, proxy, mtag, dyskinesia, adjusted
  - include age-at-onset traits
  - keep P-VALUE < 5e-8
- Strict list = reported-gene mappings only; Broad = reported + mapped.
- Gene tables saved with provenance (gene, mapping_type, study, snp, chr_id, chr_pos, p_value, trait) to `data/processed/`.
- All pipeline outputs written to `results/*.csv` in addition to Parquet (app.py still reads Parquet).
- README gene counts regenerated from `print_summary.py` (never hand-typed):
  - Parkinson's: strict 157 / broad 227 (was 162 / 257)
  - Alzheimer's: strict 93 / broad 297 (was 166 / 1166)
  - Huntington's: strict 1 / broad 1 (was 3 / 3)
  - 0 regions pass FDR correction for all diseases
- Cell-type enrichment recomputed (`results/cell_type_enrichment.csv`):
  - Alzheimer's × Microglia: OR=17.96 strict (FDR=3.5e-04), OR=10.62 broad (FDR=3.0e-05)
  - Parkinson's × Microglia: OR=10.36 strict (FDR=2.9e-03), OR=7.06 broad (FDR=1.1e-02)
- Removed 20 unused scratch scripts/outputs (`fix_readme*.py`, `check_readme*.py`, `read_*.py`, `cell_type_analysis.py` was restored from git, `sensitivity_*.py`, `verify_final.py`, `test_setup.py`, `*.txt` outputs).

### Added
- `results/trait_filters.csv` — audit trail of matched GWAS traits per disease.
- `print_summary.py` — prints the README numbers table from `results/*.csv`.
- `CHANGELOG.md` — this file.
- Honest banner on README: results are preliminary pending re-verification.

### Known limitations / pending
- Leave-one-locus-out and locus-collapse ORs not yet recomputed with the new gene lists.
- Huntington's is underpowered (1 gene passes the GWAS filter).
- Negative control (Height) uses mock `GENE_XXXX` genes; real negative-control GWAS TSV not yet obtained.
