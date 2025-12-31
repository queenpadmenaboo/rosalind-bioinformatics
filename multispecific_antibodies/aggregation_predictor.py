import os
import re
import pandas as pd
import numpy as np
from pathlib import Path

# =================================================================
# MODULE: AGGREGATION & BETA-PROPENSITY PREDICTOR
# THEORY: Beta-sheet propensity + Hydrophobicity = Aggregation Risk.
# This identifies "Sticky" regions that drive protein clumping.
# =================================================================

# --- Biophysical Propensity Scales ---
HYDRO_SCALE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 
    'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 
    'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}
BETA_SCALE = {
    'A': 0.83, 'C': 1.19, 'D': 0.54, 'E': 0.37, 'F': 1.38, 'G': 0.75, 'H': 0.87, 
    'I': 1.60, 'K': 0.74, 'L': 1.30, 'M': 1.05, 'N': 0.89, 'P': 0.55, 'Q': 1.10, 
    'R': 0.93, 'S': 0.75, 'T': 1.19, 'V': 1.70, 'W': 1.37, 'Y': 1.47
}

# --- Dynamic Path Resolution ---
CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), Path("."))

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

FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]
WINDOW = 7

def analyze_sequence(sequence):
    """Calculates APR count and Beta-propensity scores based on theory."""
    sequence = str(sequence).upper().strip()
    h_scores, b_scores = [], []
    apr_count = 0
    
    for i in range(len(sequence) - WINDOW + 1):
        sub = sequence[i : i + WINDOW]
        h = sum(HYDRO_SCALE.get(aa, 0) for aa in sub) / WINDOW
        b = sum(BETA_SCALE.get(aa, 0) for aa in sub) / WINDOW
        h_scores.append(h)
        b_scores.append(b)
        
        # Aggregation Prone Region (APR) Check
        if h > 2.0 and b > 1.2:
            apr_count += 1
            
    avg_beta = np.mean(b_scores) if b_scores else 0
    avg_hydro = np.mean(h_scores) if h_scores else 0
    
    # Theory-based Overall Score calculation
    overall = (apr_count * 10) + (avg_beta * 20) + (avg_hydro * 5)
    return apr_count, round(avg_beta, 3), round(overall, 2)

def run_predictor():
    print(f"--- ANALYZING AGGREGATION RISK: {ROOT_DIR} ---")
    results = []

    for folder in FOLDERS_TO_PROCESS:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists(): continue
        
        for file_path in folder_path.glob("*.py"):
            if file_path.name in EXCLUDE_FILES: continue
            
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                fasta_match = re.search(r'["\']{3}(.*?)["\']{3}', content, re.DOTALL)
                if not fasta_match: continue
                
                blocks = re.findall(r'>([^\n]+)\n([A-Z\n\s]+)', fasta_match.group(1).strip())
                
                for header, seq_raw in blocks:
                    clean_seq = re.sub(r'\s+', '', seq_raw)
                    apr, beta, overall = analyze_sequence(clean_seq)
                    results.append({
                        "Antibody": file_path.stem,
                        "Chain": header.strip(),
                        "Folder": folder,
                        "APR_Count": apr,
                        "Beta_Propensity": beta,
                        "Overall_Score": overall
                    })
            except Exception as e:
                print(f"Skipping {file_path.name}: {e}")

    if results:
        df = pd.DataFrame(results).sort_values(by="Overall_Score", ascending=False)
        output_file = ROOT_DIR / "Aggregation_Risk_Report.xlsx"
        df.to_excel(output_file, index=False)
        print(f"SUCCESS: Aggregation Report created: {output_file.name}")

if __name__ == "__main__":
    run_predictor()