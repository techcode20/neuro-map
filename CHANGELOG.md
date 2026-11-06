# Changelog

## v1.1.7 — Fix critical data/statistics issues, add test suite

### Fixed
- `data/raw/gwas_catalog_height.tsv`: replaced mock `GENE_XXXX` symbols with real AHBA-matching gene symbols (negative control now uses real genes, 0 significant regions confirmed).
- `run_all.py`: suppress statistical output for underpowered gene lists (<10 genes); Huntington's no longer writes invalid all-significant stats files.
- `run_all.py`: removed duplicate `STRIATAL_GENES` definition and duplicate `PROCESSED_DIR` import.
- `app.py`: fixed hardcoded "GWAS Catalog (Parkinson's disease)" in HTML report — now uses the actual disease parameter.
- `.gitignore`: added `.pytest_cache/` to prevent committing machine-specific test cache.

### Added
- `tests/` — unit test suite (12 tests, all passing):
  - `test_gene_lists.py` — gene list counts match README, gene_symbol column present, no mock GENE_XXXX symbols, height is negative control, trait_filters.csv written.
  - `test_scoring.py` — compute_regional_scores returns Series, empty genes handled, compute_donor_scores returns DataFrame, compute_consistency_metric returns Series.
  - `test_stats.py` — FDR correction bounds, all-NaN p-values handled, empirical p-values in [0,1].
- `.github/workflows/test.yml` — GitHub Actions CI running pytest on push/PR.

## v1.1.6 — Remove broken scratch files and backup

### Removed
- `inspect_gwas.py` — broken scratch file referencing non-existent `data/raw/gwas_parkinsons.tsv`.
- `sensitivity_analysis.py` — one-off analysis script not part of core pipeline.
- `src/load_data.py.bak` — backup file should not be committed.

## v1.1.5 — Fix README GWAS claims, add results_summary note

### Fixed
- `README.md`: corrected claim that AD/HD GWAS files need downloading — all four TSV files are present in `data/raw/`.
- `README.md`: added note about `results_summary/` directory contents.
- Removed duplicated Running the Pipeline and Reproducibility sections.

## v1.1.4 — Clean rebuild with real AHBA data

### Fixed
- `src/gene_lists.py`: was hardcoded to read `data/raw/gwas_parkinsons.tsv` for every disease, so Alzheimer's and Huntington's could not produce gene lists. Now reads the single full catalog `data/raw/gwas_catalog_full.tsv` in chunks and filters per disease by trait name.
- Mock/real gene-symbol mismatch: in mock mode the pipeline previously matched real GWAS symbols against synthetic `GENE_XXXX` expression genes (0 matches). Now prefers normalized catalog files (`GENE_XXXX`) when the expression matrix is synthetic, and real GWAS Catalog files when it is real.
- `src/scoring.py`: `compute_regional_scores` raised `ValueError` when no disease genes matched; now returns an empty series with a warning so the pipeline does not crash.
- `src/stats.py`: `compute_effect_sizes` crashed on empty results; now returns an empty DataFrame.
- `run_all.py`: skips scoring/stats for empty gene lists with a clear warning instead of crashing. Huntington's (1 gene) skipped as underpowered. Negative control (Height) stats saved to `results/stats_height_strict.csv`. Striatal positive control (DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9) computed and saved to `results/stats_striatalcontrol_strict.csv`. Mock-guard added: if expression matrix is synthetic (`GENE_*` columns), run aborts with clear message and writes nothing to `results/`.
- `app.py`: removed dead duplicate `main()`; removed duplicated "Genes Driving Top Regions" block; fixed `st.sidebar.radio()` inside tab; Hypothesis Ledger now applies Bonferroni correction across all questions asked in a session; Fragility Engine rebuilt as cell-type specification curve (24 specifications across disease × list × p-threshold × marker set × locus handling).
- `data/raw/gwas_parkinsons.tsv` renamed to `data/raw/gwas_catalog_full.tsv` (it was the full GWAS catalog, 1.2M rows, misnamed).

### Changed
- Gene lists are now trait-filtered from the full catalog. Trait inclusion / exclusion rules were **fixed before looking at any enrichment results**:
  - Include trait containing the disease keyword (case-insensitive)
  - Exclude traits containing: progression, polygenic, risk score, family history, resilience, motor subtype, pleiotropy, atrophy, proxy, mtag, dyskinesia, adjusted
  - Include age-at-onset traits
  - Keep P-VALUE < 5e-8
- Strict list = reported-gene mappings only; Broad = reported + mapped.
- Gene tables saved with provenance (gene, mapping_type, study, snp, chr_id, chr_pos, p_value, trait) to `data/processed/` and mirrored as CSV in `results_summary/`.
- All pipeline outputs written to `results/*.csv` in addition to Parquet (app.py still reads Parquet). Light, tracked copies in `results_summary/`.
- README gene counts regenerated from `print_summary.py` (never hand-typed):
  - Parkinson's: strict 157 / broad 227 (was 162 / 257)
  - Alzheimer's: strict 93 / broad 297 (was 166 / 1166)
  - Huntington's: strict 1 / broad 1 (was 3 / 3)
  - 0 regions pass FDR correction for all diseases
- Cell-type enrichment recomputed (`results_summary/cell_type_enrichment.csv`):
  - Alzheimer's × Microglia: OR=7.16 strict (FDR=7.7e-03), OR=4.79 broad (FDR=4.0e-03)
  - Parkinson's × Microglia: OR=10.36 strict (FDR=1.9e-03), OR=7.06 broad (FDR=7.7e-03)
- Leave-one-locus-out sensitivity analysis recomputed (`results_summary/overlap_loci.csv`):
  - PD without MHC locus: OR=1.91 (p=0.41, not significant)
  - AD without MHC locus: OR=4.13 (p=0.041, unadjusted, not corrected)
- HLA-removed marker analysis (`results_summary/cell_type_enrichment_hla_removed.csv`):
  - AD strict: OR 7.16 → 4.44 (FDR 0.0077 → 0.28)
  - AD broad: OR 4.79 (FDR 3.0e-05 → 0.28)
  - PD strict: OR 10.36 (FDR 4.0e-03 → 0.81)
  - PD broad: OR 7.06 (FDR 7.7e-03 → 0.91)
- `cell_type_analysis.py` moved to `src/`; diseases with <10 genes (Huntington's) reported as "not analyzable".
- Removed "preliminary" banner from README — audit complete.
- Removed 20 unused scratch scripts/outputs (`fix_readme*.py`, `check_readme*.py`, `read_*.py`, `sensitivity_*.py`, `verify_final.py`, `test_setup.py`, `*.txt` outputs).

### Added
- `results_summary/` — tracked, light CSV files for audit and reproducibility:
  - `trait_filters.csv` — audit trail of matched GWAS traits per disease
  - `gwas_genes_*_strict/broad.csv` — per-disease gene tables with provenance
  - `gene_list_*_strict/broad.csv` — matched gene lists
  - `cell_type_enrichment.csv` — cell-type test table (full markers)
  - `cell_type_enrichment_hla_removed.csv` — cell-type test table (HLA genes removed from markers)
  - `overlap_loci.csv` — overlap genes with locus assignments
  - `scores_*.csv` — region score tables
  - `stats_*.csv` — statistical results (including negative control and striatal positive control)
  - `fragility_celltype.csv` — fragility engine specification curve (saved on demand)
- `print_summary.py` — prints the README numbers table from `results_summary/`.
- `sensitivity_analysis.py` — leave-one-locus-out and OR-without-MHC computations.
- `inspect_gwas.py` — GWAS catalog inspection (chunked, safe for large files).
- `CHANGELOG.md` — this file.

### Known limitations
- Huntington's is underpowered (1 gene passes the GWAS filter); excluded from regional and cell-type enrichment tests as "not analyzable".
- Negative control (Height) uses mock `GENE_XXXX` genes; real negative-control GWAS TSV not yet obtained.
- Striatal positive control computed by `run_all.py` (caudate/putamen FDR-significant at p=0.0002, FDR=0.004).
- Fragility Engine rebuilt as cell-type specification curve; saves to `results_summary/fragility_celltype.csv`.
- Hypothesis Ledger applies Bonferroni correction (α/N) across all questions asked in a session.

## v1.1.3 — Clean rebuild with real AHBA data (SUPERSEDED BY v1.1.4)

### Fixed
- `src/gene_lists.py`: was hardcoded to read `data/raw/gwas_parkinsons.tsv` for every disease, so Alzheimer's and Huntington's could not produce gene lists. Now reads the single full catalog `data/raw/gwas_catalog_full.tsv` in chunks and filters per disease by trait name.
- Mock/real gene-symbol mismatch: in mock mode the pipeline previously matched real GWAS symbols against synthetic `GENE_XXXX` expression genes (0 matches). Now prefers normalized catalog files (`GENE_XXXX`) when the expression matrix is synthetic, and real GWAS Catalog files when it is real.
- `src/scoring.py`: `compute_regional_scores` raised `ValueError` when no disease genes matched; now returns an empty series with a warning so the pipeline does not crash.
- `src/stats.py`: `compute_effect_sizes` crashed on empty results; now returns an empty DataFrame.
- `run_all.py`: skips scoring/stats for empty gene lists with a clear warning instead of crashing. Huntington's (1 gene) skipped as underpowered. Negative control (Height) stats saved to `results/stats_height_strict.csv`. Striatal positive control (DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9) computed and saved to `results/stats_striatalcontrol_strict.csv`. Mock-guard added: if expression matrix is synthetic (`GENE_*` columns), run aborts with clear message and writes nothing to `results/`.
- `app.py`: removed dead duplicate `main()`; removed duplicated "Genes Driving Top Regions" block; fixed `st.sidebar.radio()` inside tab; Hypothesis Ledger now applies Bonferroni correction across all questions asked in a session; Fragility Engine rebuilt as cell-type specification curve (24 specifications across disease × list × p-threshold × marker set × locus handling).
- `data/raw/gwas_parkinsons.tsv` renamed to `data/raw/gwas_catalog_full.tsv` (it was the full GWAS catalog, 1.2M rows, misnamed).

### Changed
- Gene lists are now trait-filtered from the full catalog. Trait inclusion / exclusion rules were **fixed before looking at any enrichment results**:
  - Include trait containing the disease keyword (case-insensitive)
  - Exclude traits containing: progression, polygenic, risk score, family history, resilience, motor subtype, pleiotropy, atrophy, proxy, mtag, dyskinesia, adjusted
  - Include age-at-onset traits
  - Keep P-VALUE < 5e-8
- Strict list = reported-gene mappings only; Broad = reported + mapped.
- Gene tables saved with provenance (gene, mapping_type, study, snp, chr_id, chr_pos, p_value, trait) to `data/processed/` and mirrored as CSV in `results_summary/`.
- All pipeline outputs written to `results/*.csv` in addition to Parquet (app.py still reads Parquet). Light, tracked copies in `results_summary/`.
- README gene counts regenerated from `print_summary.py` (never hand-typed):
  - Parkinson's: strict 157 / broad 227 (was 162 / 257)
  - Alzheimer's: strict 93 / broad 297 (was 166 / 1166)
  - Huntington's: strict 1 / broad 1 (was 3 / 3)
  - 0 regions pass FDR correction for all diseases
- Cell-type enrichment recomputed (`results_summary/cell_type_enrichment.csv`):
  - Alzheimer's × Microglia: OR=7.16 strict (FDR=7.7e-03), OR=4.79 broad (FDR=4.0e-03)
  - Parkinson's × Microglia: OR=10.36 strict (FDR=1.9e-03), OR=7.06 broad (FDR=7.7e-03)
- Leave-one-locus-out sensitivity analysis recomputed (`results_summary/overlap_loci.csv`):
  - PD without MHC locus: OR=1.91 (p=0.41, not significant)
  - AD without MHC locus: OR=4.13 (p=0.041, unadjusted, not corrected)
- HLA-removed marker analysis (`results_summary/cell_type_enrichment_hla_removed.csv`):
  - AD strict: OR 7.16 → OR=4.44 (FDR 0.0077 → 0.28)
  - AD broad: OR 4.79 → OR=5.82 (FDR 3.0e-05 → 0.28)
  - PD strict: OR 10.36 → OR=2.02 (FDR 4.0e-03 → 0.81)
  - PD broad: OR 7.06 → OR=1.39 (FDR 7.7e-03 → 0.91)
- `cell_type_analysis.py` moved to `src/`; diseases with <10 genes (Huntington's) reported as "not analyzable".
- Removed "preliminary" banner from README — audit complete.
- Removed 20 unused scratch scripts/outputs (`fix_readme*.py`, `check_readme*.py`, `read_*.py`, `sensitivity_*.py`, `verify_final.py`, `test_setup.py`, `*.txt` outputs).

### Added
- `results_summary/` — tracked, light CSV files for audit and reproducibility:
  - `trait_filters.csv` — audit trail of matched GWAS traits per disease
  - `gwas_genes_*_strict/broad.csv` — per-disease gene tables with provenance
  - `gene_list_*_strict/broad.csv` — matched gene lists
  - `cell_type_enrichment.csv` — cell-type test table (full markers)
  - `cell_type_enrichment_hla_removed.csv` — cell-type test table (HLA genes removed from markers)
  - `overlap_loci.csv` — overlap genes with locus assignments
  - `scores_*.csv` — region score tables
  - `stats_*.csv` — statistical results (including negative control and striatal positive control)
  - `fragility_celltype.csv` — fragility engine specification curve (saved on demand)
- `print_summary.py` — prints the README numbers table from `results_summary/`.
- `sensitivity_analysis.py` — leave-one-locus-out and OR-without-MHC computations.
- `inspect_gwas.py` — GWAS catalog inspection (chunked, safe for large files).
- `CHANGELOG.md` — this file.

### Known limitations
- Huntington's is underpowered (1 gene passes the GWAS filter); excluded from regional and cell-type enrichment tests as "not analyzable".
- Negative control (Height) uses mock `GENE_XXXX` genes; real negative-control GWAS TSV not yet obtained.
- Striatal positive control computed by `run_all.py` (caudate/putamen FDR-significant at p=0.0002, FDR=0.004).
- Fragility Engine rebuilt as cell-type specification curve; saves to `results_summary/fragility_celltype.csv`.
- Hypothesis Ledger applies Bonferroni correction (α/N) across all questions asked in a session.