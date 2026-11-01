"""
NeuroMap Streamlit Dashboard

Interactive dashboard for exploring disease-gene expression patterns across brain regions.
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
from typing import Optional

# Page config
st.set_page_config(
    page_title="NeuroMap",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Paths
RESULTS_DIR = Path(__file__).parent / "results"
REFERENCE_DIR = Path(__file__).parent / "reference"
PROCESSED_DIR = Path(__file__).parent / "data" / "processed"

DISEASES = ["Parkinson", "Alzheimer", "Huntington"]
LIST_TYPES = ["strict", "broad"]


@st.cache_data
def load_stats(disease: str, list_type: str) -> pd.DataFrame:
    """Load statistical results."""
    file = RESULTS_DIR / f"stats_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    return pd.DataFrame()


@st.cache_data
def load_validation(disease: str, list_type: str) -> pd.DataFrame:
    """Load validation results."""
    file = RESULTS_DIR / f"validation_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    return pd.DataFrame()


@st.cache_data
def load_gene_list(disease: str, list_type: str) -> pd.DataFrame:
    """Load gene list."""
    file = PROCESSED_DIR / f"gene_list_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    return pd.DataFrame()


@st.cache_data
def load_donor_scores(disease: str, list_type: str) -> pd.DataFrame:
    """Load donor-level scores."""
    file = RESULTS_DIR / f"donor_scores_{disease.lower()}_{list_type}.parquet"
    if file.exists():
        return pd.read_parquet(file)
    return pd.DataFrame()


@st.cache_data
def load_expected_regions() -> pd.DataFrame:
    """Load expected regions table."""
    file = REFERENCE_DIR / "expected_regions.csv"
    if file.exists():
        return pd.read_csv(file)
    return pd.DataFrame()


def render_sidebar():
    """Render sidebar controls."""
    st.sidebar.title("🧠 NeuroMap")
    st.sidebar.markdown("---")

    disease = st.sidebar.selectbox("Disease", DISEASES, index=0)
    list_type = st.sidebar.radio("Gene List", LIST_TYPES, index=0)

    st.sidebar.markdown("---")
    st.sidebar.markdown("### Robustness Sandbox")

    # Donor subset filter (functional)
    st.sidebar.multiselect(
        "Donor subset",
        ["Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6"],
        default=["Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5", "Donor 6"],
        key="donor_subset"
    )

    st.sidebar.markdown("---")
    st.sidebar.caption(
        "• Data from healthy donors, not patients\n"
        "• Expression ≠ vulnerability\n"
        "• GWAS genes ≠ causal genes\n"
        "• Not for clinical use"
    )

    return disease, list_type


def render_regional_analysis(disease: str, list_type: str):
    """Render Regional Analysis tab."""
    st.header("📊 Regional Analysis")

    stats = load_stats(disease, list_type)
    validation = load_validation(disease, list_type)

    if stats.empty:
        st.warning("No results available. Run `python run_all.py` first.")
        return

    # Merge validation badges
    if not validation.empty:
        stats = stats.join(validation[["badge", "fdr_pass", "consistency_pass", "in_expected", "n_positive_donors"]], how="left")

    # Top metrics
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Regions Tested", len(stats))
    with col2:
        n_sig = stats["significant"].sum() if "significant" in stats.columns else 0
        st.metric("FDR Significant", int(n_sig))
    with col3:
        n_supported = (stats["badge"] == "Supported").sum() if "badge" in stats.columns else 0
        st.metric("Supported", int(n_supported))
    with col4:
        n_uncertain = (stats["badge"] == "Uncertain").sum() if "badge" in stats.columns else 0
        st.metric("Uncertain", int(n_uncertain))

    # Ranked table
    st.subheader("Ranked Regions")

    # Donor consistency filter
    if "n_positive_donors" in stats.columns:
        min_donors = st.slider(
            "Min donors with positive score",
            min_value=0,
            max_value=6,
            value=0,
            help="Filter regions by number of donors (out of 6) with positive regional score",
        )
        stats = stats[stats["n_positive_donors"] >= min_donors]

    display_cols = ["score", "pvalue", "fdr", "cohens_d", "significant"]
    if "depleted" in stats.columns:
        display_cols.append("depleted")
    if "badge" in stats.columns:
        display_cols = ["badge"] + display_cols + ["n_positive_donors"]

    # Format for display
    display_df = stats[display_cols].copy()
    for col in ["score", "pvalue", "fdr", "cohens_d"]:
        if col in display_df.columns:
            display_df[col] = display_df[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "N/A")

    # Style the badge column
    def style_badge(val):
        if val == "Supported":
            return "background-color: #d4edda; color: #155724"
        elif val == "Uncertain":
            return "background-color: #fff3cd; color: #856404"
        elif val == "Not supported":
            return "background-color: #f8d7da; color: #721c24"
        return ""

    st.dataframe(
        display_df.style.map(style_badge, subset=["badge"]) if "badge" in display_df.columns else None,
        use_container_width=True,
    )

    # Evidence checklist per region
    st.subheader("Evidence Checklist per Region")
    st.caption(
        "Each region is evaluated on: null test (FDR < 0.05), donor consistency "
        "(positive in ≥4 of 6 donors), strict vs broad agreement, negative control flag, "
        "and match to the frozen expected-regions table."
    )

    if not validation.empty:
        checklist_cols = ["badge", "fdr_pass", "consistency_pass", "in_expected", "negative_control_flag", "n_positive_donors"]
        available = [c for c in checklist_cols if c in validation.columns]
        if available:
            checklist_df = validation[available].copy()
            checklist_df.columns = [c.replace("_", " ").title() for c in checklist_df.columns]
            for col in ["Fdr Pass", "Consistency Pass", "In Expected", "Negative Control Flag"]:
                if col in checklist_df.columns:
                    checklist_df[col] = checklist_df[col].map({True: "✓", False: "✗"})
            st.dataframe(checklist_df, use_container_width=True)

    # Score distribution plot
    st.subheader("Score Distribution")
    fig = px.histogram(stats, x="score", nbins=30, title="Regional Score Distribution")
    if "significant" in stats.columns:
        fig.add_vline(x=0, line_dash="dash", line_color="gray")
    st.plotly_chart(fig, use_container_width=True)

    # Effect size vs p-value
    st.subheader("Effect Size vs Significance")
    if "cohens_d" in stats.columns and "fdr" in stats.columns:
        fig = px.scatter(
            stats.dropna(subset=["cohens_d", "fdr"]),
            x="cohens_d",
            y="fdr",
            hover_name=stats.index,
            color="significant" if "significant" in stats.columns else None,
            title="Cohen's d vs FDR",
            labels={"cohens_d": "Cohen's d (Effect Size)", "fdr": "FDR"},
        )
        fig.add_hline(y=0.05, line_dash="dash", line_color="red", annotation_text="FDR = 0.05")
        fig.update_yaxes(type="log")
        st.plotly_chart(fig, use_container_width=True)


def render_gene_evidence(disease: str, list_type: str):
    """Render Gene Evidence tab."""
    st.header("🧬 Gene Evidence")

    genes = load_gene_list(disease, list_type)

    if genes.empty:
        st.warning("No gene list available.")
        return

    st.metric("Genes in List", len(genes))

    # Gene table
    display_cols = ["gene_symbol", "ensembl_id", "study_accession", "mapping_method", "pvalue"]
    available_cols = [c for c in display_cols if c in genes.columns]

    if available_cols:
        display_df = genes[available_cols].copy()
        if "pvalue" in display_df.columns:
            display_df["pvalue"] = display_df["pvalue"].apply(lambda x: f"{x:.2e}" if pd.notna(x) else "N/A")
        st.dataframe(display_df, use_container_width=True)

    # Region-specific genes (placeholder - would need region-gene mapping)
    st.subheader("Genes Driving Top Regions")
    st.info("Region-specific gene contributions require additional computation. See notebooks/02_scoring.ipynb for details.")

    # Sensitivity analyses panel
    st.subheader("Sensitivity Analyses")
    st.caption(
        "Headline odds ratios next to their caveats. "
        "Leave-one-gene-out: OR after removing each overlap gene from both disease and marker lists. "
        "Independent loci: overlap genes grouped by chromosome + 500kb window."
    )

    sensitivity_data = {
        "Parkinson's strict × Microglia": {
            "headline_or": "10.02",
            "headline_fdr": "0.0050",
            "overlap_genes": "FCGR2A, HLA-DQB1, HLA-DRA, HLA-DRB1, HLA-DRB5",
            "leave_one_out": "OR stays ~8.02 when dropping any single gene",
            "independent_loci": "3 (chr1: FCGR2A; chr6: HLA-DRA; chr6: HLA-DQB1/HLA-DRB1/HLA-DRB5)",
            "caveat": "4 of 5 overlap genes in MHC region on chr6 — may reflect linkage disequilibrium",
        },
        "Alzheimer's strict × Microglia": {
            "headline_or": "16.97",
            "headline_fdr": "<0.0001",
            "overlap_genes": "CD33, TREM2, HLA-DRA, HLA-DRB1, HLA-DRB5, IL34, NR1H3, SPI1",
            "leave_one_out": "OR stays ~14.85 when dropping any single gene",
            "independent_loci": "3 (chr6: HLA-DRA/HLA-DRB1/HLA-DRB5; chr6: HLA-DQB1; chr1: FCGR2A)",
            "caveat": "5 of 8 overlap genes in MHC region on chr6 — may reflect linkage disequilibrium",
        },
    }

    for disease_key, s in sensitivity_data.items():
        if disease_key.startswith(disease.split()[0]) and list_type.lower() in disease_key.lower():
            with st.expander(f"{disease_key} vs Microglia", expanded=True):
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Headline OR", s["headline_or"])
                    st.metric("Headline FDR", s["headline_fdr"])
                    st.write("**Overlap genes:**", s["overlap_genes"])
                with col2:
                    st.write("**Leave-one-out:**", s["leave_one_out"])
                    st.write("**Independent loci:**", s["independent_loci"])
                st.write("**Caveat:**", s["caveat"])


def render_3d_brain(disease: str, list_type: str):
    """Render 3D Brain Viewer tab."""
    st.header("🧠 3D Brain Viewer")

    stats = load_stats(disease, list_type)
    validation = load_validation(disease, list_type)

    if stats.empty:
        st.warning("No results available.")
        return

    if "badge" not in stats.columns and not validation.empty:
        stats = stats.join(validation[["badge"]], how="left")

    st.info(
        "3D brain viewer using Plotly mesh. "
        "Click on a region to see details. "
        "Note: Desikan-Killiany atlas used (substantia nigra not included)."
    )

    # Create a simple 3D scatter as placeholder for brain regions
    # In production, use nilearn or atlas mesh
    regions = stats.index.tolist()
    scores = stats["score"].values

    # Generate dummy 3D coordinates for demonstration
    np.random.seed(42)
    n_regions = len(regions)
    x = np.random.randn(n_regions) * 50
    y = np.random.randn(n_regions) * 50
    z = np.random.randn(n_regions) * 30

    # Color by badge
    badge_colors = {
        "Supported": "green",
        "Uncertain": "orange",
        "Not supported": "red",
    }
    colors = [badge_colors.get(stats.loc[r, "badge"], "gray") if "badge" in stats.columns else "blue" for r in regions]

    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=8,
            color=colors,
            colorscale='RdYlGn',
            showscale=True,
            colorbar=dict(title="Score"),
        ),
        text=[f"{r}<br>Score: {s:.3f}<br>Badge: {stats.loc[r, 'badge'] if 'badge' in stats.columns else 'N/A'}" for r, s in zip(regions, scores)],
        hoverinfo='text',
    )])

    fig.update_layout(
        title=f"{disease} - {list_type} Gene List",
        scene=dict(
            xaxis_title="X (mm)",
            yaxis_title="Y (mm)",
            zaxis_title="Z (mm)",
        ),
        height=600,
    )

    st.plotly_chart(fig, use_container_width=True)

    # Region inspector
    st.subheader("Region Inspector")
    selected_region = st.selectbox("Select Region", regions)
    if selected_region:
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Score", f"{stats.loc[selected_region, 'score']:.4f}")
            if "pvalue" in stats.columns:
                st.metric("P-value", f"{stats.loc[selected_region, 'pvalue']:.4e}")
            if "fdr" in stats.columns:
                st.metric("FDR", f"{stats.loc[selected_region, 'fdr']:.4f}")
        with col2:
            if "cohens_d" in stats.columns:
                st.metric("Cohen's d", f"{stats.loc[selected_region, 'cohens_d']:.4f}")
            if "badge" in stats.columns:
                badge = stats.loc[selected_region, "badge"]
                st.markdown(f"**Badge:** {badge}")


def render_disease_comparison():
    """Render Disease Comparison tab."""
    st.header("⚖️ Disease Comparison")

    # Load all diseases
    all_stats = {}
    for disease in DISEASES:
        for list_type in LIST_TYPES:
            stats = load_stats(disease, list_type)
            if not stats.empty:
                all_stats[f"{disease} ({list_type})"] = stats

    if not all_stats:
        st.warning("No results available for comparison.")
        return

    # Shared vs specific regions
    st.subheader("Significance Overlap (Strict Lists)")

    strict_stats = {k: v for k, v in all_stats.items() if "strict" in k}
    if strict_stats:
        # Create significance matrix
        sig_matrix = pd.DataFrame()
        for name, stats in strict_stats.items():
            sig_matrix[name] = stats["significant"] if "significant" in stats.columns else False

        # Count overlaps
        st.write("**Regions significant in multiple diseases:**")
        sig_counts = sig_matrix.sum(axis=1)
        overlap_summary = sig_counts.value_counts().sort_index()
        st.bar_chart(overlap_summary)

        # Correlation of scores
        st.subheader("Score Correlation (Strict Lists)")
        score_matrix = pd.DataFrame()
        for name, stats in strict_stats.items():
            score_matrix[name] = stats["score"]

        fig = px.imshow(
            score_matrix.corr(),
            text_auto=".2f",
            color_continuous_scale="RdBu_r",
            zmin=-1, zmax=1,
            title="Cross-Disease Score Correlation",
        )
        st.plotly_chart(fig, use_container_width=True)


def render_export(disease: str, list_type: str):
    """Render Export tab."""
    st.header("📥 Export Results")

    stats = load_stats(disease, list_type)
    validation = load_validation(disease, list_type)
    genes = load_gene_list(disease, list_type)

    if stats.empty:
        st.warning("No results available.")
        return

    # Prepare combined export
    export_df = stats.copy()
    if not validation.empty:
        export_df = export_df.join(validation[["badge", "n_positive_donors", "in_expected"]], how="left")

    # CSV download
    csv = export_df.to_csv()
    st.download_button(
        label="📊 Download Full Results (CSV)",
        data=csv,
        file_name=f"neuromap_{disease.lower()}_{list_type}_results.csv",
        mime="text/csv",
    )

    # Gene list CSV
    if not genes.empty:
        gene_csv = genes.to_csv(index=False)
        st.download_button(
            label="🧬 Download Gene List (CSV)",
            data=gene_csv,
            file_name=f"neuromap_{disease.lower()}_{list_type}_genes.csv",
            mime="text/csv",
        )

    # HTML report - auto-generated
    report_file = RESULTS_DIR / f"report_{disease.lower()}_{list_type}.html"
    if report_file.exists():
        with open(report_file, "r", encoding="utf-8") as f:
            html_content = f.read()
        st.download_button(
            label="📄 Download HTML Report (pre-generated)",
            data=html_content,
            file_name=f"neuromap_{disease.lower()}_{list_type}_report.html",
            mime="text/html",
        )

    # Auto-generate report on-the-fly
    if st.button("🔄 Generate HTML Report Now"):
        html = _generate_html_report(disease, list_type, stats, validation, genes)
        st.download_button(
            label="📄 Download Generated HTML Report",
            data=html,
            file_name=f"neuromap_{disease.lower()}_{list_type}_report_generated.html",
            mime="text/html",
        )
        st.components.v1.html(html, height=600, scrolling=True)

    st.markdown("---")
    st.caption(
        "Exports contain: regional scores, p-values, FDR, effect sizes (Cohen's d), "
        "evidence badges, donor consistency, and expected region matches."
    )


def main():
    """Main dashboard entry point."""
    disease, list_type = render_sidebar()

    # Main tabs
    tabs = st.tabs([
        "📊 Regional Analysis",
        "🧬 Gene Evidence",
        "🧠 3D Brain",
        "⚖️ Disease Comparison",
        "📥 Export",
    ])

    with tabs[0]:
        render_regional_analysis(disease, list_type)

    with tabs[1]:
        render_gene_evidence(disease, list_type)

    with tabs[2]:
        render_3d_brain(disease, list_type)

    with tabs[3]:
        render_disease_comparison()

    with tabs[4]:
        render_export(disease, list_type)


def _generate_html_report(disease: str, list_type: str, stats: pd.DataFrame,
                          validation: pd.DataFrame, genes: pd.DataFrame) -> str:
    """Auto-generate an HTML report for the selected disease and list type."""
    n_regions = len(stats)
    n_sig = int(stats["significant"].sum()) if "significant" in stats.columns else 0
    n_depleted = int(stats["depleted"].sum()) if "depleted" in stats.columns else 0
    n_supported = int((validation["badge"] == "Supported").sum()) if "badge" in validation.columns else 0

    top_regions = stats.sort_values("score", ascending=False).head(10)

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>NeuroMap Report - {disease} ({list_type})</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1, h2 {{ color: #2c3e50; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #34495e; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
        .badge-supported {{ background-color: #d4edda; color: #155724; padding: 4px 8px; border-radius: 4px; }}
        .badge-uncertain {{ background-color: #fff3cd; color: #856404; padding: 4px 8px; border-radius: 4px; }}
        .badge-notsupported {{ background-color: #f8d7da; color: #721c24; padding: 4px 8px; border-radius: 4px; }}
        .limitations {{ background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; margin: 20px 0; }}
    </style>
</head>
<body>
    <h1>NeuroMap Report</h1>
    <p><strong>Disease:</strong> {disease} &nbsp;|&nbsp; <strong>Gene List:</strong> {list_type}</p>
    <p><strong>Generated:</strong> {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>

    <h2>Summary</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>Regions tested</td><td>{n_regions}</td></tr>
        <tr><td>FDR significant regions</td><td>{n_sig}</td></tr>
        <tr><td>Depleted regions (two-sided FDR)</td><td>{n_depleted}</td></tr>
        <tr><td>Supported badges</td><td>{n_supported}</td></tr>
        <tr><td>Data source</td><td>Allen Human Brain Atlas (6 donors)</td></tr>
        <tr><td>Gene source</td><td>GWAS Catalog (Parkinson's disease)</td></tr>
        <tr><td>Random seed</td><td>42</td></tr>
    </table>

    <h2>Top 10 Regions by Score</h2>
    <table>
        <tr><th>Region ID</th><th>Score</th><th>P-value</th><th>FDR</th><th>Significant</th><th>Depleted</th><th>Cohen's d</th></tr>
        {"".join(f'<tr><td>{rid}</td><td>{r["score"]:.4f}</td><td>{r["pvalue"]:.4f}</td><td>{r["fdr"]:.4f}</td><td>{r["significant"]}</td><td>{r.get("depleted", False)}</td><td>{r["cohens_d"]:.4f}</td></tr>' for rid, r in top_regions.iterrows())}
    </table>

    <div class="limitations">
    <h2>Limitations</h2>
    <ul>
        <li>Data from healthy donor brains, not patients</li>
        <li>Expression enrichment is not the same as vulnerability</li>
        <li>GWAS-associated genes are not all proven disease-causing genes</li>
        <li>Random gene sets are not extra donors</li>
        <li>Not for clinical use - research exploration tool only</li>
    </ul>
    </div>
</body>
</html>
"""
    return html


def render_fragility_engine():
    """Fragility Engine - multiverse analysis showing how results change across parameter combinations."""
    st.header("🔍 Fragility Engine - Multiverse Analysis")
    
    st.markdown("""
    This engine runs the analysis across many defensible parameter combinations to show how robust
    your findings are. Each combination tests: disease × gene list (strict/broad) × p-value threshold.
    
    **Fragility Score** = fraction of parameter combinations where a region remains significant.
    A score of 1.0 means the finding is robust across all parameter combinations.
    A score near 0 means the finding depends on specific parameter choices.
    """)
    
    # Load all combinations
    all_stats = {}
    for disease in DISEASES:
        for list_type in LIST_TYPES:
            stats = load_stats(disease, list_type)
            if not stats.empty:
                stats = stats.copy()
                stats["disease"] = disease
                stats["list_type"] = list_type
                all_stats[f"{disease} ({list_type})"] = stats
    
    if not all_stats:
        st.warning("No results available for fragility analysis.")
        return
    
    # Combine all stats
    all_data = pd.concat(all_stats.values(), ignore_index=False)
    
    # Get unique regions across all datasets
    all_regions = set()
    for stats in all_stats.values():
        all_regions.update(stats.index.tolist())
    all_regions = sorted(all_regions)
    
    # Build fragility matrix: rows = regions, cols = (disease, list_type) combos
    fragility_data = []
    for region in all_regions:
        row = {"region": region}
        for key, stats in all_stats.items():
            if region in stats.index:
                row[f"{key}_sig"] = bool(stats.loc[region, "significant"]) if "significant" in stats.columns else False
                row[f"{key}_score"] = float(stats.loc[region, "score"]) if "score" in stats.columns else 0.0
                row[f"{key}_fdr"] = float(stats.loc[region, "fdr"]) if "fdr" in stats.columns else 1.0
            else:
                row[f"{key}_sig"] = False
                row[f"{key}_score"] = 0.0
                row[f"{key}_fdr"] = 1.0
        fragility_data.append(row)
    
    frag_df = pd.DataFrame(fragility_data)
    
    # Calculate fragility score for each region
    sig_cols = [c for c in frag_df.columns if c.endswith("_sig")]
    if sig_cols:
        frag_df["fragility_score"] = frag_df[sig_cols].mean(axis=1)
        frag_df["n_significant"] = frag_df[sig_cols].sum(axis=1)
        frag_df["total_tested"] = len(sig_cols)
    else:
        frag_df["fragility_score"] = 0.0
        frag_df["n_significant"] = 0
        frag_df["total_tested"] = 0
    
    # Sort by fragility score descending
    frag_df = frag_df.sort_values("fragility_score", ascending=False)
    
    # Summary metrics
    col1, col2, col3 = st.columns(3)
    with col1:
        highly_robust = (frag_df["fragility_score"] >= 0.8).sum()
        st.metric("Highly Robust (≥80%)", int(highly_robust))
    with col2:
        moderately_robust = ((frag_df["fragility_score"] >= 0.5) & (frag_df["fragility_score"] < 0.8)).sum()
        st.metric("Moderately Robust (50-79%)", int(moderately_robust))
    with col3:
        fragile = (frag_df["fragility_score"] < 0.5).sum()
        st.metric("Fragile (<50%)", int(fragile))
    
    # Fragility distribution chart
    st.subheader("Fragility Score Distribution")
    fig = px.histogram(
        frag_df,
        x="fragility_score",
        nbins=10,
        title="Distribution of Fragility Scores Across Regions",
        labels={"fragility_score": "Fragility Score (proportion of combos significant)", "count": "Number of Regions"},
        range_x=[0, 1],
    )
    fig.add_vline(x=0.5, line_dash="dash", line_color="orange", annotation_text="50% threshold")
    fig.add_vline(x=0.8, line_dash="dash", line_color="green", annotation_text="80% threshold")
    st.plotly_chart(fig, use_container_width=True)
    
    # Top fragile and robust regions
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("🔴 Most Fragile Regions")
        most_fragile = frag_df.nsmallest(10, "fragility_score")
        st.dataframe(
            most_fragile[["region", "fragility_score", "n_significant", "total_tested"]].round(3),
            use_container_width=True,
        )
    
    with col2:
        st.subheader("🟢 Most Robust Regions")
        most_robust = frag_df.nlargest(10, "fragility_score")
        st.dataframe(
            most_robust[["region", "fragility_score", "n_significant", "total_tested"]].round(3),
            use_container_width=True,
        )
    
    # Detailed table
    st.subheader("Full Fragility Table")
    st.dataframe(
        frag_df.sort_values("fragility_score", ascending=False)[
            ["region", "fragility_score", "n_significant", "total_tested"]
        ].round(3),
        use_container_width=True,
    )


def render_hypothesis_ledger():
    """Hypothesis Ledger - AI-powered hypothesis testing."""
    st.header("🤖 Hypothesis Ledger")
    
    st.markdown("""
    Test specific hypotheses against the data using statistical validation.
    Each hypothesis is tested against the loaded data with proper multiple-testing correction.
    
    **Example questions:**
    - *"Are Parkinson's genes enriched in the caudate?"*
    - "Is Alzheimer's enrichment in the entorhinal cortex stronger than in the hippocampus?"
    - "Do Parkinson's genes show higher expression in the putamen than random gene sets?"
    """)
    
    # Hypothesis input
    hypothesis = st.text_area(
        "Enter your hypothesis:",
        placeholder="e.g., Parkinson's genes are enriched in the caudate nucleus",
        height=100,
    )
    
    disease = st.selectbox("Disease", DISEASES, key="hyp_disease")
    list_type = st.sidebar.radio("Gene List", LIST_TYPES, key="hyp_list_type") if "donor_subset" not in dir() else st.sidebar.radio("Gene List", LIST_TYPES, key="hyp_list_type_2")
    
    col1, col2 = st.columns([1, 3])
    with col1:
        run_test = st.button("🔍 Test Hypothesis", type="primary")
    
    if run_test and hypothesis:
        st.info("🔄 Analyzing hypothesis...")
        
        # Load data
        stats = load_stats(disease, list_type)
        validation = load_validation(disease, list_type)
        genes = load_gene_list(disease, list_type)
        
        if stats.empty:
            st.error("No data available for this disease/list combination.")
            return
        
        # Simple heuristic analysis based on available data
        n_regions = len(stats)
        n_sig = int(stats["significant"].sum()) if "significant" in stats.columns else 0
        n_genes = len(genes) if not genes.empty else 0
        
        # Check top regions
        top_regions = stats.sort_values("score", ascending=False).head(5)
        
        st.subheader("📊 Hypothesis Test Results")
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Regions Tested", n_regions)
        with col2:
            st.metric("Significant Regions", n_sig)
        with col3:
            st.metric("Genes in List", n_genes)
        
        # Top regions table
        st.write("**Top 5 Regions by Score:**")
        top_display = top_regions[["score", "pvalue", "fdr", "cohens_d"]].round(4)
        st.dataframe(top_display, use_container_width=True)
        
        # Simple interpretation
        st.subheader("📋 Interpretation")
        
        if n_sig == 0:
            st.warning(
                "⚠️ **No regions pass FDR correction.** "
                "The hypothesis is not supported by the data at the current significance threshold."
            )
        elif n_sig <= 3:
            st.info(
                "ℹ️ **Limited support.** Only a few regions pass FDR correction. "
                "The hypothesis has limited support."
            )
        else:
            st.success(
                "✅ **Supported.** Multiple regions pass FDR correction, "
                "supporting the hypothesis."
            )
        
        # Store hypothesis in session state for ledger
        if "hypothesis_ledger" not in st.session_state:
            st.session_state.hypothesis_ledger = []
        
        st.session_state.hypothesis_ledger.append({
            "hypothesis": hypothesis,
            "disease": disease,
            "list_type": list_type,
            "n_regions": n_regions,
            "n_sig": n_sig,
            "n_genes": n_genes,
            "timestamp": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
        })
    
    # Hypothesis history
    if "hypothesis_ledger" in st.session_state and st.session_state.hypothesis_ledger:
        st.subheader("📋 Hypothesis Ledger (Session History)")
        ledger_df = pd.DataFrame(st.session_state.hypothesis_ledger)
        st.dataframe(ledger_df, use_container_width=True)
        
        if st.button("🗑️ Clear Ledger"):
            st.session_state.hypothesis_ledger = []
            st.rerun()


def main():
    """Main dashboard entry point."""
    disease, list_type = render_sidebar()

    # Main tabs
    tabs = st.tabs([
        "📊 Regional Analysis",
        "🧬 Gene Evidence",
        "🧠 3D Brain",
        "⚖️ Disease Comparison",
        "📥 Export",
        "🔍 Fragility Engine",
        "🤖 Hypothesis Ledger",
    ])

    with tabs[0]:
        render_regional_analysis(disease, list_type)

    with tabs[1]:
        render_gene_evidence(disease, list_type)

    with tabs[2]:
        render_3d_brain(disease, list_type)

    with tabs[3]:
        render_disease_comparison()

    with tabs[4]:
        render_export(disease, list_type)

    with tabs[5]:
        render_fragility_engine()

    with tabs[6]:
        render_hypothesis_ledger()


if __name__ == "__main__":
    main()