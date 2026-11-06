with open('app.py', 'r', encoding='utf-8') as f:
    content = f.read()
print(f'Total lines: {len(content.splitlines())}')
for name in ['def render_fragility_engine():', 'def render_hypothesis_ledger():', 'def render_sidebar():', 'def render_regional_analysis(disease: str, list_type: str):', 'def render_gene_evidence(disease: str, list_type: str):', 'def render_3d_brain(disease: str, list_type: str):', 'def render_disease_comparison():', 'def render_export(disease: str, list_type: str):', 'def _generate_html_report(disease: str, list_type: str, stats: pd.DataFrame,', 'def main():']:
    idx = content.find(name)
    print(f'{name}: found at {idx}')