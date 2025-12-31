import os
import re
import pandas as pd
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =================================================================
# MODULE: THERMAL STABILITY & GURUPRASAD INSTABILITY PREDICTOR
# =================================================================

"""
STABILITY ANALYSIS & BIOPHYSICAL THEORY:

1. DEFINITION OF ALIPHATIC: 
   Refers to non-polar, 'oily' amino acid side chains (Alanine, Valine, Isoleucine, 
   and Leucine). In the 'Oil Drop' model, these cluster to form a hydrophobic 
   core, shielding themselves from water to stabilize the antibody.

2. ALIPHATIC INDEX (AI) & Tm:
   AI measures the relative volume occupied by aliphatic side chains. A higher 
   AI indicates a denser 'anchor' in the protein core, directly increasing 
   the Predicted Tm (melting temperature).

3. INSTABILITY INDEX (Guruprasad et al., 1990):
   This methodology predicts stability based on 'Dipeptide' logic.
   - SCORE < 40: Classified as 'Stable'.
   - SCORE > 40: Classified as 'Unstable'.

4. RISK SCORE GRADING:
   - A: High Tm (>75C) & Guruprasad Stable.
   - F: Low Tm AND Very High Instability (>50).
"""

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

def assign_risk_score(tm, instability):
    if tm > 75 and instability < 40: return "A (Excellent)"
    elif tm < 65 and instability > 50: return "F (High Risk)"
    else: return "C (Moderate)"

def calculate_stability_metrics(sequence):
    clean_seq = "".join([aa for aa in str(sequence).upper() if aa in "ACDEFGHIKLMNPQRSTVWY"])
    if len(clean_seq) < 10: return None

    analyser = ProteinAnalysis(clean_seq)
    count = analyser.count_amino_acids()
    total = len(clean_seq)
    
    # Aliphatic Index: (100 * (A + 2.9V + 3.9(I+L))) / total
    aliphatic_index = (100 * (count['A'] + 2.9*count['V'] + 3.9*(count['I'] + count['L']))) / total
    instability_index = analyser.instability_index()
    cys_percent = (count['C'] / total) * 100
    
    # Predicted Tm Estimation logic
    predicted_tm = 55.0 + (0.5 * aliphatic_index) + (2.0 * cys_percent) - (instability_index * 0.1)
    
    return {
        "Aliphatic_Index": round(aliphatic_index, 2),
        "Instability_Index": round(instability_index, 2),
        "Predicted_Tm_C": round(predicted_tm, 1),
        "Risk_Score": assign_risk_score(predicted_tm, instability_index)
    }

def run_stability_predictor():
    print(f"--- ANALYZING THERMAL STABILITY: {ROOT_DIR} ---")
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
                    metrics = calculate_stability_metrics(clean_seq)
                    if metrics:
                        row = {"Antibody": file_path.stem, "Chain": header.strip(), "Folder": folder}
                        row.update(metrics)
                        results.append(row)
            except Exception as e:
                print(f"Skipping {file_path.name}: {e}")

    if results:
        df = pd.DataFrame(results).sort_values(by="Predicted_Tm_C", ascending=False)
        output_file = ROOT_DIR / "Thermal_Stability_Report.xlsx"
        df.to_excel(output_file, index=False)
        print(f"SUCCESS: Stability Report created: {output_file.name}")

if __name__ == "__main__":
    run_stability_predictor()