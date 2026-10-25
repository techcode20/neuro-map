"""Automatic HTML report generation."""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional
from datetime import datetime

RESULTS_DIR = Path(__file__).parent.parent / "results"
REFERENCE_DIR = Path(__file__).parent.parent / "reference"


def generate_html_report(
    disease: str,
    list_type: str = "strict",
    output_dir: Optional[Path] = None,
) -> Path:
    """
    Generate an HTML report for a disease analysis.

    Parameters
    ----------
    disease : str
        Disease name.
    list_type : str
        'strict' or 'broad'.
    output_dir : Path, optional
        Output directory. Defaults to results/.

    Returns
    -------
    Path
        Path to generated HTML file.
    """
    if output_dir is None:
        output_dir = RESULTS_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load all results
    stats_file = RESULTS_DIR / f"stats_{disease.lower()}_{list_type}.parquet"
    validation_file = RESULTS_DIR / f"validation_{disease.lower()}_{list_type}.parquet"
    scores_file = RESULTS_DIR / f"scores_{disease.lower()}_{list_type}.parquet"
    donor_file = RESULTS_DIR / f"donor_scores_{disease.lower()}_{list_type}.parquet"
    gene_list_file = Path(__file__).parent.parent / "data" / "processed" / f"gene_list_{disease.lower()}_{list_type}.parquet"

    stats = pd.read_parquet(stats_file) if stats_file.exists() else pd.DataFrame()
    validation = pd.read_parquet(validation_file) if validation_file.exists() else pd.DataFrame()
    scores = pd.read_parquet(scores_file) if scores_file.exists() else pd.DataFrame()
    donor_scores = pd.read_parquet(donor_file) if donor_file.exists() else pd.DataFrame()
    gene_list = pd.read_parquet(gene_list_file) if gene_list_file.exists() else pd.DataFrame()

    # Load expected regions
    expected_file = REFERENCE_DIR / "expected_regions.csv"
    expected = pd.read_csv(expected_file) if expected_file.exists() else pd.DataFrame()
    expected_disease = expected[expected["disease"] == disease] if not expected.empty else pd.DataFrame()

    # Generate HTML
    html = _build_html_report(
        disease=disease,
        list_type=list_type,
        stats=stats,
        validation=validation,
        scores=scores,
        donor_scores=donor_scores,
        gene_list=gene_list,
        expected_regions=expected_disease,
    )

    output_file = output_dir / f"report_{disease.lower()}_{list_type}.html"
    with open(output_file, "w") as f:
        f.write(html)

    print(f"Generated report: {output_file}")
    return output_file


def _build_html_report(
    disease: str,
    list_type: str,
    stats: pd.DataFrame,
    validation: pd.DataFrame,
    scores: pd.DataFrame,
    donor_scores: pd.DataFrame,
    gene_list: pd.DataFrame,
    expected_regions: pd.DataFrame,
) -> str:
    """Build the HTML content."""

    # Top regions by score
    if not stats.empty:
        top_regions = stats.nlargest(10, "score")[["score", "pvalue", "fdr", "cohens_d", "significant"]]
    else:
        top_regions = pd.DataFrame()

    # Evidence badge summary
    if not validation.empty:
        badge_counts = validation["badge"].value_counts()
    else:
        badge_counts = pd.Series()

    # Gene list summary
    n_genes = len(gene_list) if not gene_list.empty else 0
    gene_summary = ""
    if not gene_list.empty:
        gene_summary = "<ul>"
        for _, gene in gene_list.head(20).iterrows():
            study = gene.get("study_accession", "N/A")
            mapping = gene.get("mapping_method", "N/A")
            gene_summary += f"<li>{gene['gene_symbol']} (study: {study}, mapping: {mapping})</li>"
        if n_genes > 20:
            gene_summary += f"<li>... and {n_genes - 20} more</li>"
        gene_summary += "</ul>"

    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>NeuroMap Report: {disease} ({list_type})</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
        h1, h2, h3 {{ color: #2c3e50; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #34495e; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .badge-supported {{ background-color: #27ae60; color: white; padding: 4px 8px; border-radius: 4px; }}
        .badge-uncertain {{ background-color: #f39c12; color: white; padding: 4px 8px; border-radius: 4px; }}
        .badge-not-supported {{ background-color: #e74c3c; color: white; padding: 4px 8px; border-radius: 4px; }}
        .metric {{ display: inline-block; margin: 10px 20px 10px 0; padding: 10px; background: #ecf0f1; border-radius: 4px; }}
        .limitations {{ background: #fff3cd; border: 1px solid #ffc107; padding: 15px; border-radius: 4px; margin: 20px 0; }}
    </style>
</head>
<body>
    <h1>NeuroMap Analysis Report</h1>
    <p><strong>Disease:</strong> {disease} | <strong>Gene List:</strong> {list_type} | <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

    <div class="limitations">
        <h3>⚠️ Limitations</h3>
        <ul>
            <li>Data from healthy donor brains, not patients. Shows where disease genes are expressed, not where damage occurs.</li>
            <li>Expression enrichment ≠ vulnerability. Top-scoring region is not automatically the most damaged one.</li>
            <li>GWAS-associated genes are not all proven disease-causing genes. Evidence type and mapping method kept for every gene.</li>
            <li>Random gene sets are not extra donors. Permutation tests do not add biological samples.</li>
            <li>This is a research exploration tool and is <strong>not for clinical use</strong>.</li>
        </ul>
    </div>

    <h2>Summary</h2>
    <div class="metric"><strong>Genes in list:</strong> {n_genes}</div>
    <div class="metric"><strong>Regions tested:</strong> {len(stats) if not stats.empty else 0}</div>
    <div class="metric"><strong>Significant (FDR < 0.05):</strong> {int(stats['significant'].sum()) if not stats.empty and 'significant' in stats.columns else 0}</div>
    <div class="metric"><strong>Supported regions:</strong> {badge_counts.get('Supported', 0)}</div>
    <div class="metric"><strong>Uncertain regions:</strong> {badge_counts.get('Uncertain', 0)}</div>
    <div class="metric"><strong>Not supported regions:</strong> {badge_counts.get('Not supported', 0)}</div>

    <h2>Top 10 Regions by Score</h2>
    {top_regions.to_html(float_format=lambda x: f'{x:.4f}') if not top_regions.empty else '<p>No data available</p>'}

    <h2>Evidence Badges</h2>
    {validation[['badge', 'fdr_pass', 'consistency_pass', 'in_expected', 'n_positive_donors']].to_html(
        formatters={'badge': lambda x: f'<span class="badge-{x.lower().replace(" ", "-")}">{x}</span>'}
    ) if not validation.empty else '<p>No validation data available</p>'}

    <h2>Expected Regions (from literature)</h2>
    {expected_regions[['region', 'source', 'confidence']].to_html(index=False) if not expected_regions.empty else '<p>No expected regions defined</p>'}

    <h2>Gene List</h2>
    {gene_summary}

    <h2>Donor Consistency</h2>
    {donor_scores.to_html(float_format=lambda x: f'{x:.4f}') if not donor_scores.empty else '<p>No donor-level data available</p>'}

    <hr>
    <p><em>Generated by NeuroMap pipeline. For research use only.</em></p>
</body>
</html>
"""
    return html


def generate_all_reports() -> list:
    """Generate reports for all diseases and list types."""
    diseases = ["Parkinson", "Alzheimer", "Huntington"]
    list_types = ["strict", "broad"]
    reports = []

    for disease in diseases:
        for list_type in list_types:
            stats_file = RESULTS_DIR / f"stats_{disease.lower()}_{list_type}.parquet"
            if stats_file.exists():
                report_file = generate_html_report(disease, list_type)
                reports.append(report_file)

    return reports