# NeuroMap Reproduction Guide

This document describes the exact steps to reproduce the real-data (non-mock) results from a fresh clone.

## Prerequisites

- Python 3.12+
- pip
- A web browser with access to https://human.brain-map.org (login required for downloads)
- ~3 GB free disk space for AHBA microarray data
- ~4 GB free RAM for processing

## Step 0: Clone and install

```bash
git clone <repo-url> NeuroMap
cd NeuroMap
python -m venv .venv
# Windows:
.venv\Scripts\Activate.ps1
# macOS/Linux:
source .venv/bin/activate
python -m pip install -r requirements.txt
```

## Step 1: Download Allen Human Brain Atlas microarray data (MANUAL - requires login)

The Allen Brain Map requires a logged-in browser session. Downloads are large (~400 MB per donor, ~2.5 GB total).

1. Open https://human.brain-map.org in your browser and **sign in** (create an account if needed).
2. Download each donor ZIP file by pasting the URL into your browser address bar:

| Donor | File ID | Download URL | Save as |
|-------|---------|--------------|---------|
| 9861 | 178238387 | https://human.brain-map.org/api/v2/well_known_file_download/178238387 | `data/raw/allen/zips/donor9861.zip` |
| 10021 | 178238373 | https://human.brain-map.org/api/v2/well_known_file_download/178238373 | `data/raw/allen/zips/donor10021.zip` |
| 12876 | 178238359 | https://human.brain-map.org/api/v2/well_known_file_download/178238359 | `data/raw/allen/zips/donor12876.zip` |
| 14380 | 178238316 | https://human.brain-map.org/api/v2/well_known_file_download/178238316 | `data/raw/allen/zips/donor14380.zip` |
| 15496 | 178238266 | https://human.brain-map.org/api/v2/well_known_file_download/178238266 | `data/raw/allen/zips/donor15496.zip` |
| 15697 | 178236545 | https://human.brain-map.org/api/v2/well_known_file_download/178236545 | `data/raw/allen/zips/donor15697.zip` |

The `data/raw/allen/zips/` directory already exists in the repo.

3. Extract each ZIP into its corresponding donor folder:

```powershell
# Windows PowerShell (run from repo root)
foreach ($d in @("9861","10021","12876","14380","15496","15697")) {
    $target = "data/raw/allen/microarray/normalized_microarray_donor$d"
    if (-not (Test-Path $target)) { New-Item -ItemType Directory -Path $target -Force | Out-Null }
    Copy-Item "data/raw/allen/zips/donor$d.zip" "$target\" -Force
    Expand-Archive "data/raw/allen/zips/donor$d.zip" -DestinationPath "$target" -Force
}

# Verify: each folder should contain 5 CSV files
dir data/raw/allen/microarray/normalized_microarray_donor9861
```

Expected files per donor: `MicroarrayExpression.csv` (~960 MB), `Ontology.csv`, `PACall.csv`, `Probes.csv`, `SampleAnnot.csv`

4. Verify zip integrity (optional):

```powershell
python -c "import zipfile; print(zipfile.ZipFile('data/raw/allen/zips/donor9861.zip').testzip())"
# Should print: None
```

## Step 2: Download GWAS Catalog Parkinson's associations (MANUAL)

1. Go to https://www.ebi.ac.uk/gwas/efotraits/EFO_0000435 (Parkinson's disease trait page).
2. Click **Download associations** (the full tab-separated file).
3. Save as `data/raw/gwas_parkinsons.tsv` (this replaces any mock file).
4. Note the download date for your records.

The file should have columns including: `STUDY ACCESSION`, `SNPS`, `REPORTED GENE(S)`, `MAPPED_GENE`, `P-VALUE`, `DISEASE/TRAIT`.

## Step 3: Run the pipeline with real data

```powershell
# Windows PowerShell
$env:REAL_DATA=1
python run_all.py
```

This will:
1. Load the 6-donor AHBA expression data via abagen (local files, no download)
2. Build gene lists from the downloaded GWAS Catalog TSV
3. Run negative control, permutation tests, FDR correction, validation
4. Generate HTML reports for each disease/list combination

Expected runtime: ~10-20 minutes for 6 donors (mostly abagen processing).

## Step 4: Verify results match committed numbers

After running, verify:

```powershell
# Check expression matrix shape (should be 83 regions x ~15,633 genes)
python -c "import pandas as pd; d=pd.read_parquet('data/processed/allen_expression_desikan-killiany.parquet'); print(d.shape)"

# Check Parkinson's strict gene count (should be 162)
python -c "import pandas as pd; d=pd.read_parquet('data/processed/gene_list_parkinson_strict.parquet'); print(len(d))"

# Check negative control shows 0 significant regions
python -c "import pandas as pd; d=pd.read_parquet('results/stats_parkinson_strict.parquet'); print(d[d['significant']==True].shape)"
```

## Step 5: Launch dashboard

```powershell
python -m streamlit run app.py
```

Then open http://localhost:8501 in your browser.

## Expected Results Summary

| Metric | Value |
|--------|-------|
| Expression matrix | 83 regions × 15,633 genes |
| Parkinson's strict genes | 162 |
| Parkinson's broad genes | 257 |
| Alzheimer's strict genes | 166 |
| Alzheimer's broad genes | 1,166 |
| Huntington's strict genes | 3 |
| Huntington's broad genes | 3 |
| Parkinson's significant regions (FDR) | 0 |
| Negative control significant regions | 0 |
| Caudate/putamen cohens_d | ≈ -2.4 (descriptive, not FDR-significant) |
| Microglia enrichment (AD strict) | OR=16.97, FDR<0.0001 (circular - driven by 8 shared genes) |
| Microglia enrichment (PD strict) | OR=10.02, FDR=0.005 (circular - driven by 5 shared genes) |
| Dopaminergic neuron enrichment (PD) | OR=0, zero overlap, low power |

## Notes

- `data/` and `results/` are gitignored - they are not in the repository. You must generate them by running the pipeline.
- The `data/raw/allen/zips/` directory exists but is empty - you must download the ZIP files manually (Step 1).
- The GWAS Catalog REST API is deprecated (returns 410 Gone). The pipeline falls back to reading the downloaded TSV file directly.
- Real AHBA processing takes ~10-20 minutes. Mock data (REAL_DATA not set) runs in seconds but uses synthetic gene symbols.