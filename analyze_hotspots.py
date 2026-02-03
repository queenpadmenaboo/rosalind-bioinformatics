import os
import re
import pandas as pd
from pathlib import Path

# =================================================================
# MODULE: ADVANCED ANTIBODY HOTSPOT AUDIT
# PURPOSE: Identify PTMs and Structural Risks in Multispecifics
# =================================================================

# --- Path Configuration ---
CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), None)
CATEGORY_FOLDERS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

# --- FULL EXCLUSION LIST ---
EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py', 'analyze_hotspots.py', 'aggregation_predictor.py', 'Aggregation_Risk_Report.xlsx',
    'antibody_diagnostic_tool.py', 'Antibody_Comparison_Report_2025.xlsx', 'build_pnas_shadow.py',
    'compare_hydrophobicity.py', 'developability_hotspots.xlsx', 'mAb_truth_engine',
    'mAb_Truth_Engine_Master.xlsx', 'MISSING_ANTIBODIES_LOG.xlsx', 'pnas_validator.py',
    'PNAS_VS_CALCULATIONS.xlsx'
}

# --- Theory: Post-Translational Modification (PTM) Risks ---
HOTSPOT_MOTIFS = {
    'N_Glycosylation': r'N[^P][ST]',          # Theory: Sugar attachment can block antigen binding
    'Deamidation_NG': r'NG',                  # Theory: Asparagine to Aspartic acid (Charge change)
    'Isomerization_DG': r'DG',                # Theory: Aspartic acid flip (Backbone instability)
    'Hydrophobic_Cluster': r'[VILFW]{5,}',    # Theory: Surface-exposed 'greasy' patches
    'Multispecific_Linker': r'([G]{3,4}S){2,}' # Theory: Linker flexibility/vulnerability
}

def scan_hotspots(sequence, name="Unknown"):
    """Performs scan of the sequence for biophysical risks."""
    seq = str(sequence).upper().strip()
    report = {'Antibody': name}
    total_score = 0
    
    for motif_name, pattern in HOTSPOT_MOTIFS.items():
        matches = len(re.findall(pattern, seq))
        report[motif_name] = matches
        
        # Weighting logic
        if motif_name == 'N_Glycosylation':
            total_score += (matches * 5)
        else:
            total_score += (matches * 2)
            
    report['Total_Hotspot_Score'] = total_score
    return report

def run_hotspot_audit():
    if not ROOT_DIR:
        print("Error: Root path not found.")
        return

    print(f"--- STARTING HOTSPOT AUDIT: {ROOT_DIR} ---")
    all_results = []

    for folder in CATEGORY_FOLDERS:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists():
            continue

        for py_file in folder_path.glob('*.py'):
            if py_file.name in EXCLUDE_FILES:
                continue
            
            try:
                # SAFE TEXT READ
                with open(py_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Match text inside triple quotes
                fasta_match = re.search(r'["\']{3}(.*?)["\']{3}', content, re.DOTALL)
                if not fasta_match:
                    continue
                
                fasta_str = fasta_match.group(1).strip()
                
                # Parse the FASTA headers and sequences
                blocks = re.findall(r'>([^\n]+)\n([A-Z\n\s]+)', fasta_str)
                
                for header, seq_raw in blocks:
                    clean_seq = re.sub(r'\s+', '', seq_raw)
                    analysis = scan_hotspots(clean_seq, py_file.stem)
                    analysis['Chain'] = header.strip()
                    analysis['Format_Folder'] = folder
                    all_results.append(analysis)
            
            except Exception as e:
                print(f"Skipping {py_file.name}: {e}")

    if all_results:
        df = pd.DataFrame(all_results)
        output_file = ROOT_DIR / "developability_hotspots.xlsx"
        df.to_excel(output_file, index=False)
        print(f"SUCCESS: Hotspot Audit complete. Results: {output_file.name}")
    else:
        print("No antibody sequences found to analyze.")

if __name__ == "__main__":
    run_hotspot_audit()