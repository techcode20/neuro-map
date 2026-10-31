# NeuroMap

An interactive tool that maps where disease-associated genes are expressed across the human brain, tests whether the pattern is statistically credible, and explains the result.

## Aim

Build a Python application that takes a neurological disease (Parkinson's first, then Alzheimer's and Huntington's) and answers three questions:

1. Which brain regions show unusually strong expression of the disease's associated genes?
2. Is that pattern statistically meaningful, or could it appear by chance?
3. Does the result agree with published knowledge about which regions the disease affects?

The output is an interactive dashboard with a 3D brain map, ranked regions, evidence badges, gene-level explanations, and exportable results.

## Key Findings

- **Bulk regional expression:** No significant enrichment of Parkinson's GWAS genes in any brain region (0 regions pass FDR). Caudate/putamen show depletion (cohens_d ≈ -2.4) but it does not survive two-sided FDR — descriptive only.
- **Cell-type layer:** Microglia marker overlap present for both diseases but rests on about three loci, mostly in the MHC; not strong evidence for microglial involvement. AD-microglia OR=16.97 → OR=10.85 without MHC locus (still significant); PD-microglia OR=10.02 → OR=2.01 without MHC locus (not significant).
- **Validation:** Negative control returns 0 significant regions; calibration test returns 0/100 false positives; striatal gene set positive control returns FDR-significant caudate/putamen enrichment.

## Validation Summary

| Test | Result | Interpretation |
|------|--------|----------------|
| Negative control (Height GWAS) | 0 significant regions | Method calibrated, no false positives |
| Calibration test | 0/100 runs with any FDR hit | Null model correct |
| Positive control (striatal genes) | Caudate/putamen FDR-significant | Pipeline detects real signal |
| PD regional enrichment | 0 regions significant | Honest null result |
| AD–microglia enrichment | OR=16.97 → OR=10.85 without MHC | Remains significant; driven by non-MHC genes |
| PD–microglia enrichment | OR=10.02 → OR=2.01 without MHC | Not significant; driven by MHC locus |
| PD–dopaminergic neuron enrichment | Zero overlap, CI not estimable | Low power — cannot conclude absence |

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

- **Cell-type layer**: test enrichment in cell types such as dopaminergic neurons or medium spiny neurons
- Spatially aware null models
- Gene co-expression network view

## Cell-Type Enrichment Results

Fisher's exact test (FDR-corrected across all disease × cell-type tests) of GWAS genes against published cell-type marker lists:

| Disease | List | Cell Type | Odds Ratio | p-value | FDR | Significant |
|---------|------|-----------|------------|---------|-----|-------------|
| Alzheimer's | strict | Microglia | 16.97 | <0.0001 | <0.0001 | Yes |
| Parkinson's | strict | Microglia | 10.02 | 0.0002 | 0.0050 | Yes |
| Alzheimer's | broad | Microglia | 3.57 | 0.0005 | 0.0071 | Yes |
| Parkinson's | broad | Microglia | 6.21 | 0.0019 | 0.0199 | Yes |
| Alzheimer's | strict | Oligodendrocytes | 9.47 | 0.0051 | 0.0427 | Yes |

### Sensitivity Analyses

**Leave-one-locus-out (microglia enrichment):**
- **AD strict (8 overlap genes):** OR=16.97 with all overlap → OR=10.85 without MHC locus (HLA-DRA, HLA-DRB1, HLA-DRB5, HLA-DQB1). Remaining overlap genes: CD33, IL34, NR1H3, SPI1, TREM2 (5 genes across chr19, chrX, chr12, chr11, chr6p21 — independent chromosomes). Enrichment remains significant (p=0.00017), not driven by MHC.
- **PD strict (5 overlap genes):** OR=10.02 with all overlap → OR=2.01 without MHC locus (p=0.399, not significant). Only FCGR2A (chr1) remains as overlap. PD-microglia enrichment is driven by MHC locus genes.

**Independent loci among overlap genes (500kb windows):**
- **AD overlap genes (8):** HLA-DRA, HLA-DRB1, HLA-DRB5 on chr6 (1 MHC locus), HLA-DQB1 on chr6 (1 MHC locus), FCGR2A on chr1 (1 locus). CD33, IL34, NR1H3, SPI1, TREM2 not in significant GWAS hits. → 3 independent loci total.
- **PD overlap genes (5):** HLA-DRA, HLA-DRB1, HLA-DRB5, HLA-DQB1 on chr6 (2 MHC loci), FCGR2A on chr1 (1 locus). → 3 independent loci total.
- Most overlap genes sit in the MHC region (chr6), where genes are in linkage disequilibrium — this is a hypothesis, not something directly tested.

**Dopaminergic neuron CI:**
- OR=0 with zero overlap (23 markers vs 162 PD genes). CI not estimable ([0, inf]). Low power — cannot conclude absence.

**Locus collapse (Parkinson microglia):**
- OR 10.02 → 10.12 when collapsing to one gene per 500kb locus — enrichment remains robust.

### What Did Not Hold Up

| Claim | Result | Why |
|-------|--------|-----|
| AD–microglia enrichment | OR=16.97 → OR=10.85 without MHC locus (still significant) | Remaining overlap genes (CD33, IL34, NR1H3, SPI1, TREM2) across 4 chromosomes — not driven by MHC alone |
| PD–microglia enrichment | OR=10.02 → OR=2.01 without MHC locus (not significant) | Enrichment driven by MHC locus genes on chr6; only FCGR2A (chr1) remains |
| PD–dopaminergic neuron enrichment | Zero overlap, CI not estimable | Low power (23 markers vs 162 genes) — cannot conclude absence |
| Caudate/putamen depletion for PD | Does not survive two-sided FDR | Descriptive only, not a significant finding |
| Huntington's cell-type analysis | Only 3 genes | Underpowered |

### Positive Control (regional method)

The real positive control for the regional method is the striatal gene set (DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9). Caudate and putamen are FDR-significant (p=0.0002, FDR=0.004), confirming the regional pipeline can detect a true signal when one exists.

## What the Story Becomes

"Disease genes are not enriched in the expected brain regions by bulk expression, a validated and calibrated method confirms that, and the cell-type layer shows which cell types the genes actually belong to - but with important caveats about overlap-gene independence and low power."

This is a stronger, more honest narrative than a forced positive result:
- Bulk regional expression shows **no significant enrichment** for Parkinson's (0 regions pass FDR)
- Caudate/putamen show **depletion** (not enrichment), consistent across all 6 donors (cohens_d ≈ -2.4), but this **does not survive two-sided FDR correction** - treat as descriptive only
- The **cell-type layer** reveals that PD and AD GWAS signals track microglia (neuroinflammation), not cell-type-specific vulnerability - but:
  - **AD–microglia**: OR=16.97 → OR=10.85 without MHC locus (still significant). Driven by non-MHC genes (CD33, IL34, NR1H3, SPI1, TREM2) across 4 different chromosomes. Likely real.
  - **PD–microglia**: OR=10.02 → OR=2.01 without MHC locus (not significant). Driven by MHC locus genes on chr6. Likely reflects LD in the MHC region rather than independent microglial signals.
  - **HLA caveat**: HLA genes are expressed in several cell types besides microglia (e.g., dendritic cells, B cells), so microglia-enrichment claims based on HLA overlap may partly reflect antigen-presentation genes in general.

**Microglia marker overlap is present for both diseases but rests on about three loci, mostly in the MHC; it is not strong evidence for microglial involvement.**

## Installation

```bash
git clone <repo-url>
cd NeuroMap
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

For real-data reproduction steps, see [REPRODUCE.md](REPRODUCE.md).

## Running the Pipeline

```bash
# Run the full pipeline
python run_all.py

# Launch the dashboard
python -m streamlit run app.py
```

### Dashboard Robustness Sandbox

The dashboard sidebar includes sliders for GWAS p-value threshold, number of random gene sets, and donor subset. These are **UI controls only** in this release — they set session state but do not re-run the pipeline with different parameters. To change these parameters, edit `src/gene_lists.py` (P_THRESHOLD) and `src/stats.py` (n_sets) and re-run `run_all.py`. A fully parameterized dashboard is planned for a future release.

## Reproducibility

- All random seeds are fixed (42)
- Package versions are pinned in `requirements.txt`
- The expected regions table is frozen in `reference/expected_regions.csv`
- `run_all.py` rebuilds all results from raw data

## License

MIT License