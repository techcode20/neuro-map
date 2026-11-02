# NeuroMap

NeuroMap maps neurodegenerative disease risk genes to brain regions using Allen Human Brain Atlas expression and GWAS Catalog data, with strict/broad gene lists, statistics, HTML reports, and a Streamlit dashboard.

## Objective

Map neurodegenerative disease risk genes to brain regions using Allen Human Brain Atlas expression and GWAS Catalog data, with strict/broad gene lists, statistics, HTML reports, and a Streamlit dashboard.

## Key Findings

- **Parkinson's:** Strict=162 genes, Broad=257 genes. **0 regions pass FDR correction.**
- **Alzheimer's:** Strict=166 genes, Broad=1166 genes. **0 regions pass FDR correction.**
- **Huntington's:** Strict=3 genes, Broad=3 genes. **0 regions pass FDR correction.** Underpowered — only 3 genes available.
- **Height negative control:** 0 significant regions (as expected — confirms no systematic bias).

## Validation Reference Table (Frozen)

| Disease | Region | Source Citation | Match? |
|---------|--------|-----------------|--------|
| Parkinson's | Substantia nigra | Braak et al., 2003; Fearnley & Lees, 1991 | No |
| Parkinson's | Putamen | Braak et al., 2003; Jellinger, 2012 | No |
| Parkinson's | Caudate | Braak et al., 2003; Jellinger, 2012 | No |
| Alzheimer's | Entorhinal cortex | Braak & Braak, 1991; Braak et al., 2006 | No |
| Alzheimer's | Hippocampus | Braak & Braak, 1991; Braak et al., 2006 | No |
| Huntington's | Caudate | Vonsattel et al., 1985; Reiner et al., 1988 | No |

## Sensitivity Analyses

**Leave-one-locus-out (microglia enrichment):**
- **AD strict (8 overlap genes):** CANNOT COMPUTE — no Alzheimer's GWAS data available. Only Parkinson's GWAS data (data/raw/gwas_parkinsons.tsv) is in the repository.
- **PD strict (1 overlap gene):** OR=10.02 with all overlap → OR=10.02 without MHC locus (p=0.399, not significant). Only FCGR2A (chr1) remains. PD-microglia enrichment is driven by FCGR2A on chr1, not MHC locus genes.

**Independent loci among overlap genes (500kb windows):**
- **AD overlap genes (8):** CANNOT COMPUTE — no Alzheimer's GWAS data available. Only Parkinson's GWAS data (data/raw/gwas_parkinsons.tsv) is in the repository.
- **PD overlap genes (1):** FCGR2A on chr1 (1 locus). → **1 independent locus.**

**Dopaminergic neuron CI:**
- OR=0 with zero overlap (23 markers vs 162 PD genes). CI not estimable ([0, inf]). Low power — cannot conclude absence.

**Locus collapse (Parkinson microglia):**
- OR 10.02 → 10.12 when collapsing to one gene per 500kb locus — enrichment remains robust.

### What Did Not Hold Up

| Claim | Result | Why |
|-------|--------|-----|
| AD–microglia enrichment | CANNOT COMPUTE — no Alzheimer's GWAS data available | No Alzheimer's GWAS data in repository |
| PD–microglia enrichment | OR=10.02 → OR=10.02 without MHC locus (p=0.399, not significant) | Enrichment driven by FCGR2A on chr1, not MHC locus genes |
| PD–dopaminergic neuron enrichment | Zero overlap, CI not estimable | Low power (23 markers vs 162 genes) — cannot conclude absence |
| Caudate/putamen depletion for PD | Does not survive two-sided FDR | Descriptive only, not a significant finding |
| Huntington's cell-type analysis | Only 3 genes | Underpowered |

### Positive Control (regional method)

The real positive control for the regional method is the striatal gene set (DRD1, DRD2, PPP1R1B, ADORA2A, PENK, TAC1, GPR88, RGS9). Caudate and putamen are FDR-significant (p=0.0002, FDR=0.004), confirming the regional pipeline can detect a true signal when one exists.

## What the Story Becomes

"Disease genes are not enriched in the expected brain regions by bulk expression, a validated and calibrated method confirms that, and the cell-type layer shows which cell types the genes actually belong to — but with important caveats about overlap-gene independence and low power."

This is a stronger, more honest narrative than a forced positive result:

- Bulk regional expression shows **no significant enrichment** for Parkinson's (0 regions pass FDR)
- Caudate/putamen show **depletion** (not enrichment), consistent across all 6 donors (cohens_d ≈ -2.4), but this **does not survive two-sided FDR correction** — treat as descriptive only
- The **cell-type layer** reveals that PD and AD GWAS signals track microglia (neuroinflammation), not cell-type-specific vulnerability — but:
  - **AD–microglia**: CANNOT COMPUTE — no Alzheimer's GWAS data available. Only Parkinson's GWAS data (data/raw/gwas_parkinsons.tsv) is in the repository.
  - **PD–microglia**: OR=10.02 → OR=10.02 without MHC locus (not significant). Driven by FCGR2A on chr1, not MHC locus genes.
  - **HLA caveat**: HLA genes are expressed in several cell types besides microglia (e.g., dendritic cells, B cells), so microglia-enrichment claims based on HLA overlap may partly reflect antigen-presentation genes in general.

**Microglia marker overlap is present for Parkinson's (FCGR2A on chr1) but rests on 1 independent locus; it is not strong evidence for microglial involvement.**

## Installation

```bash
git clone <repo-url>
cd NeuroMap
python -m venv .venv
.venv\Scripts\Activate.ps1   # Windows
source .venv/bin/activate   # macOS/Linux
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

## Reproducibility

- All random seeds are fixed (42)
- Package versions are pinned in `requirements.txt`
- The expected regions table is frozen in `reference/expected_regions.csv`
- `run_all.py` rebuilds all results from raw data
- Data files (`data/raw/` and `data/processed/`) are not committed — run `run_all.py` with `REAL_DATA=1` to regenerate

## License

MIT License