# NeuroMap

An interactive tool that maps where disease-associated genes are expressed across the human brain, tests whether the pattern is statistically credible, and explains the result.

## Aim

Build a Python application that takes a neurological disease (Parkinson's first, then Alzheimer's and Huntington's) and answers three questions:

1. Which brain regions show unusually strong expression of the disease's associated genes?
2. Is that pattern statistically meaningful, or could it appear by chance?
3. Does the result agree with published knowledge about which regions the disease affects?

The output is an interactive dashboard with a 3D brain map, ranked regions, evidence badges, gene-level explanations, and exportable results.

## The Problem

Disease-associated genes exist throughout the body, yet diseases damage specific brain regions. For example, Parkinson's prominently affects the substantia nigra, Alzheimer's the entorhinal cortex and connected memory regions, and Huntington's the striatum. These are simplified descriptions, since each disease affects multiple regions and cell types.

NeuroMap explores one angle on this: whether disease genes are unusually expressed in the regions known to be affected.

## Scope and Honest Limits

**State these in the README and in the dashboard.**

- The data come from healthy donor brains, not patients. The tool shows where disease genes are expressed, not where damage occurs.
- Expression enrichment is not the same as vulnerability. The top-scoring region is not automatically the most damaged one.
- GWAS-associated genes are not all proven disease-causing genes. Evidence type and mapping method are kept for every gene.
- Random gene sets are not extra donors. Permutation tests do not add biological samples.
- It is a research exploration tool and is **not for clinical use**.
- No machine-learning model is added without a proper training and evaluation dataset.

## Data Sources

| Source | What it provides | How it is used |
|--------|------------------|----------------|
| Allen Human Brain Atlas | Gene expression from six adult control brains, with anatomical annotations | Region-by-gene expression matrices |
| GWAS Catalog | Published disease-variant associations, reported and mapped genes, study accessions | Disease gene lists with evidence provenance |
| abagen (Python toolbox) | Processing of Allen data: probe-to-gene collapsing, donor normalization, region alignment | Builds the expression matrix reproducibly |
| Allen single-nucleus data (advanced) | Cell-type expression | Cell-type enrichment layer |
| Published literature | Which regions each disease affects | Frozen, cited reference table used for validation |

## System Architecture

The project is a pipeline of six stages feeding one dashboard:

1. **Ingest**: abagen fetches and processes Allen data. The GWAS Catalog file or REST API supplies disease associations.
2. **Clean and align**: match gene identifiers (symbol versus Ensembl ID), map samples to atlas regions, check missing values, keep donor-level matrices.
3. **Score**: compute a regional disease-gene score for each region.
4. **Test**: matched random gene sets, empirical p-values, FDR correction, effect sizes, donor consistency, negative control.
5. **Validate**: compare against the frozen expected-regions table and assign evidence badges.
6. **Present**: Streamlit dashboard, 3D viewer, exports and the automatic report.

Data flow: raw data goes into `data/raw`, processed matrices into `data/processed`, results tables into `results`, and the dashboard reads only from processed and results files.

## Technology Stack

| Tool | Purpose |
|------|---------|
| Python | Main language |
| Pandas, NumPy | Data handling |
| abagen | Allen atlas processing |
| SciPy | Statistics and permutations |
| Statsmodels | Benjamini-Hochberg FDR correction |
| Matplotlib, Seaborn, Plotly | Charts and interactive plots |
| nilearn or Plotly mesh | 3D brain surface viewer |
| NetworkX (optional) | Gene co-expression network view |
| Streamlit | Dashboard |
| Jupyter or Google Colab | Exploration and validation notebooks |
| Git and GitHub | Version control and portfolio |

## Project Structure

```
NeuroMap/
├── data/raw/            original downloads, never edited
├── data/processed/      clean matrices and gene lists
├── results/             scores, p-values, badges
├── notebooks/           01_explore, 02_scoring, 03_validation
├── src/                 load_data.py, gene_lists.py, scoring.py,
│                        stats.py, validation.py, report.py
├── reference/           expected_regions.csv (cited, frozen)
├── app.py               Streamlit dashboard
├── run_all.py           one command rebuilds all results
├── requirements.txt     pinned versions
└── README.md
```

## Key Decisions

| Decision | Choice |
|----------|--------|
| Brain atlas | Desikan-Killiany first. It includes caudate and putamen but not the substantia nigra, so note the gap or build a custom atlas with the substantia nigra from the Allen structure ontology |
| Gene rule | Genome-wide significant GWAS hits (p < 5e-8) with a documented mapping rule |
| Gene lists | A strict version and a broad version, for sensitivity checks |
| Expected regions | Cited table (disease, region, source), written before running anything |
| Score | Mean z-scored expression across disease genes per region |
| Null model | 5,000 random gene sets matched on expression level |
| Correction | FDR across all regions tested |
| Random seed | Fixed (42), so results are reproducible |

## Statistical Method

**Step A: Observed score.** For each region r, the score is the mean of z-scored expression values over the G selected disease genes. The exact standardization must be justified for the dataset, because raw expression values from different genes and probes are not on comparable scales.

**Step B: Null distribution.** Generate many random gene sets of the same size G from an eligible background, matched on expression level where feasible. Compute each set's regional score with the identical procedure.

**Step C: Empirical p-value.** For a one-sided test of unusually high expression: p = (1 + number of null scores at least as large as the observed score) / (1 + number of random sets).

**Step D: Multiple testing.** Apply Benjamini-Hochberg FDR across the prespecified family of regional tests.

**Step E: Effect size.** Report the observed score relative to the null mean and standard deviation, so results are not just red or green.

**Step F: Spatial effects (advanced).** The random-gene-set test addresses gene-set specificity. It does not address spatial autocorrelation between neighbouring regions. A spatial null model answers that separate question, and the dashboard shows both p-values side by side.

## Evidence Engine and Badges

Every region receives a checklist and a badge. Define the rules in the README before seeing results.

**Checklist per region**

1. Random-gene-set test passes (FDR below 0.05)
2. Pattern is consistent across donors
3. Pattern holds for both strict and broad gene lists
4. Negative-control trait does not flag the same region
5. Region matches the frozen expected-regions table
6. Spatial null result (advanced, when available)

**Badge rule (example)**

- **Supported**: FDR below 0.05 and consistent in at least 4 of 6 donors
- **Uncertain**: only one of those conditions holds
- **Not supported**: neither holds

Document mismatches with published evidence as findings. Never change the method to force agreement.

## Features

### Core

- **Regional Analysis**: ranked regions, scores, adjusted p-values, effect sizes, badge
- **Gene Evidence**: which genes drive a region's score, with study accession and reported-versus-mapped type
- **Export**: CSV of results and an automatic HTML report

### Standout

- Evidence checklist per region
- Interactive 3D brain viewer with click-to-inspect regions
- Robustness sandbox: sliders for GWAS threshold, strict or broad list, number of random sets and donor subset, with live updates
- Disease Comparison: shared versus disease-specific regions and similarity between disease profiles

### Advanced (choose one or two)

- Cell-type layer: test enrichment in cell types such as dopaminergic neurons or medium spiny neurons
- Spatially aware null models
- Gene co-expression network view

## Installation

```bash
git clone <repo-url>
cd NeuroMap
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## Running the Pipeline

```bash
# Run the full pipeline
python run_all.py

# Launch the dashboard
python -m streamlit run app.py
```

## Reproducibility

- All random seeds are fixed (42)
- Package versions are pinned in `requirements.txt`
- The expected regions table is frozen in `reference/expected_regions.csv`
- `run_all.py` rebuilds all results from raw data

## License

MIT License