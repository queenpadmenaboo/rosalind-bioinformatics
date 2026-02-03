"""
BIOINFORMATICS PIPELINE: ANTIBODY STABILITY & VISCOSITY MASTER CALCULATOR
==========================================================================================
UPDATED DEC 2025: Added Homodimer Multiplier Logic for Multispecific Benchmarks.
==========================================================================================
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import string
import os

# ==========================================================================================
# --- Path Configuration ---
# ==========================================================================================

ROOT_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies\shadow_benchmarks")

if not ROOT_DIR.exists():
    raise FileNotFoundError(f"No valid ROOT_DIR found at {ROOT_DIR}")

# --- Project Specific Configuration ---
CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.csv', 'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py', 'antibody_diagnostic_tool.py', 'build_pnas_shadow.py',
    'recover_truth_engine.py', 'merge_truth.py', 'pnas_validator.py'
}

STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

def parse_py_file(filepath: Path) -> tuple:
    """Parses Python files to extract metadata, headers and sequences."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    sequences = {}
    current_header = None
    current_seq = []
    
    # Extract metadata from the triple-quote block
    is_homodimer = "homodimer" in content.lower()
    
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_seq)
            current_header = line[1:]
            current_seq = []
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from', 'FORMAT:')):
            if line.lower() not in ('na', 'n/a'):
                current_seq.append(line)
    
    if current_header:
        sequences[current_header] = ''.join(current_seq)
    
    return sequences, is_homodimer

def analyze_sequence_features(sequence: str) -> dict:
    """Calculates biophysical features and Viscosity Risk."""
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not clean_seq:
        return None
        
    analysis = ProteinAnalysis(clean_seq)
    pi = analysis.isoelectric_point()
    charge_55 = analysis.charge_at_pH(5.5)
    
    # Base risk calculation (based on absolute charge per complex)
    abs_charge = abs(charge_55)
    if abs_charge < 10.0:
        v_risk = "High Risk (Neutral/Sticky)"
    elif abs_charge < 20.0:
        v_risk = "Moderate Risk"
    else:
        v_risk = "Low Risk (Good Repulsion)"

    features = {
        'calculated_pI': pi,
        'net_charge_pH5.5': charge_55,
        'viscosity_risk': v_risk,
        'gravy': analysis.gravy(),
        'instability_index': analysis.instability_index(),
        'aromaticity': analysis.aromaticity()
    }
    
    aa_fractions = analysis.amino_acids_percent
    for aa in STANDARD_AAS:
         features[f'AA_{aa}_fraction'] = aa_fractions.get(aa, 0.0) 
         
    return features

def main():
    print("=" * 75)
    print(f"ANTIBODY FEATURE & VISCOSITY CALCULATION - TARGET: {ROOT_DIR.name}")
    print("=" * 75)
    
    results = []
    
    for folder in CATEGORY_FOLDERS:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists():
            continue
            
        print(f"Processing folder: {folder}")
        for py_file in folder_path.glob('*.py'):
            if py_file.name in EXCLUDE_FILES:
                continue
                
            try:
                sequences, is_homodimer = parse_py_file(py_file)
                if not sequences:
                    continue

                h_seqs = [s for h, s in sequences.items() if any(k in h.lower() for k in ['heavy', 'vh', '_h', 'chain a'])]
                l_seqs = [s for h, s in sequences.items() if any(k in h.lower() for k in ['light', 'vl', '_l', 'chain b'])]

                # Fallback for generic headers
                if not h_seqs and not l_seqs:
                    all_seqs = list(sequences.values())
                    if len(all_seqs) >= 2:
                        h_seqs, l_seqs = [all_seqs[0]], [all_seqs[1]]
                    elif len(all_seqs) == 1:
                        h_seqs, l_seqs = [all_seqs[0]], []

                # --- APPLY STOICHIOMETRY LOGIC ---
                if folder == 'Whole_mAb':
                    # Standard IgG (H2L2)
                    full_sequence = "".join(h_seqs * 2 + l_seqs * 2) 
                    multiplier = 1
                elif is_homodimer:
                    # Acimtamig Style: Double the entire sequence set
                    full_sequence = "".join(h_seqs) + "".join(l_seqs)
                    multiplier = 2
                    print(f"  [METRIC] Applying x2 multiplier for Homodimer: {py_file.stem}")
                else:
                    # Standard Fragments (scFv, VHH, etc.)
                    full_sequence = "".join(h_seqs) + "".join(l_seqs)
                    multiplier = 1

                features = analyze_sequence_features(full_sequence)
                
                if features:
                    # Apply multipliers to extensive properties
                    features['net_charge_pH5.5'] *= multiplier
                    total_len = len(full_sequence) * multiplier
                    
                    # Re-evaluate viscosity risk based on MULTIPLIED charge
                    abs_charge = abs(features['net_charge_pH5.5'])
                    if abs_charge < 10.0:
                        features['viscosity_risk'] = "High Risk (Neutral/Sticky)"
                    elif abs_charge < 20.0:
                        features['viscosity_risk'] = "Moderate Risk"
                    else:
                        features['viscosity_risk'] = "Low Risk (Good Repulsion)"

                    data = {
                        'antibody': py_file.stem,
                        'folder': folder,
                        'is_homodimer': is_homodimer,
                        'total_length': total_len
                    }
                    data.update(features)
                    results.append(data)
            except Exception as e:
                print(f"Skipping {py_file.name} due to error: {e}")
                continue
    
    if results:
        df = pd.DataFrame(results)
        output_file_path = ROOT_DIR / 'sequence_features.xlsx'
        
        with pd.ExcelWriter(output_file_path, engine='xlsxwriter') as writer:
            df.to_excel(writer, sheet_name='Antibody Features', index=False)
            workbook  = writer.book
            worksheet = writer.sheets['Antibody Features']
            
            high_risk_format = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
            worksheet.freeze_panes(1, 0)
            worksheet.autofilter(0, 0, len(df), len(df.columns) - 1)
            
            for i, col in enumerate(df.columns):
                max_len = max(df[col].astype(str).map(len).max(), len(col)) + 4 
                worksheet.set_column(i, i, max_len)

            risk_col_idx = df.columns.get_loc('viscosity_risk')
            worksheet.conditional_format(1, risk_col_idx, len(df), risk_col_idx, {
                'type': 'text', 'criteria': 'containing', 'value': 'High', 'format': high_risk_format
            })

        print(f"SUCCESS: Report saved as {output_file_path.name}")
        print(f"Processed {len(df)} antibodies.")
    
    print("=" * 75)

if __name__ == "__main__":
    main()