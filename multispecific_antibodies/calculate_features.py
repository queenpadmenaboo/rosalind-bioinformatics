"""
BIOINFORMATICS PIPELINE: ANTIBODY STABILITY & VISCOSITY MASTER CALCULATOR
==========================================================================================
VISCOSITY & ELECTROSTATIC THEORY:
1. THE REPULSION PRINCIPLE: 
   Antibodies are large, charged molecules. For high-concentration formulations 
   (e.g., >100 mg/mL), we rely on "Electrostatic Repulsion" to keep them apart. 
   If molecules have the same strong charge, they bounce off each other like magnets.

2. NET CHARGE @ pH 5.5 (Formulation pH):
   Most antibodies are stored in slightly acidic buffers (pH 5.5). 
   - A Net Charge > +20 or < -20 is ideal for "Low Viscosity."
   - A Net Charge near ZERO (-10 to +10) indicates a "High Viscosity Risk" because 
     the molecules lack the force to repel each other, leading to "molecular tangling."

3. ISOELECTRIC POINT (pI):
   The pH where the net charge is zero. If the formulation pH is near the pI, 
   the antibody is "neutral" and highly likely to aggregate or become syrupy.

4. CALCULATION METHOD:
   The script uses the ProtParam module (BioPython) to calculate the pI and then 
   solves the Henderson-Hasselbalch equation for the specific pH of 5.5 to 
   determine the net charge of the entire protein complex.

5. SMART PARSING LOGIC:
   Handles non-standard headers (e.g., Emicizumab, Aducanumab) by looking for 
   expanded keywords (_h, _l, vh, vl) and implementing a fallback that assumes 
   the sequence order (Heavy then Light) if no keywords are found.
==========================================================================================

REQUIRES:
    pip install biopython pandas xlsxwriter
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import string
import os

# ==========================================================================================
# --- Path Configuration (TOGGLE BETWEEN PROJECT AND SHADOW) ---
# ==========================================================================================

# --- OPTION A: PROJECT MODE (Your Designs) ---
# ROOT_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")

# --- OPTION B: SHADOW MODE (The 137 PNAS Benchmarks) ---
ROOT_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies\shadow_benchmarks")

# ==========================================================================================

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
    'recover_truth_engine.py', 'merge_truth.py'
}

STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

def parse_py_file(filepath: Path) -> dict:
    """Parses Python files to extract FASTA headers and sequences."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    sequences = {}
    current_header = None
    current_seq = []
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_seq)
            current_header = line[1:]
            current_seq = []
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            if line.lower() not in ('na', 'n/a'):
                current_seq.append(line)
    
    if current_header:
        sequences[current_header] = ''.join(current_seq)
    return sequences

def analyze_sequence_features(sequence: str) -> dict:
    """Calculates biophysical features and Viscosity Risk based on Electrostatic Theory."""
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not clean_seq:
        return None
        
    analysis = ProteinAnalysis(clean_seq)
    pi = analysis.isoelectric_point()
    charge_55 = analysis.charge_at_pH(5.5)
    
    # VISCOSITY RISK LOGIC
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
                sequences = parse_py_file(py_file)
                if not sequences:
                    continue

                h_seqs = [s for h, s in sequences.items() if any(k in h.lower() for k in ['heavy', 'vh', '_h', 'chain a'])]
                l_seqs = [s for h, s in sequences.items() if any(k in h.lower() for k in ['light', 'vl', '_l', 'chain b'])]

                if not h_seqs and not l_seqs:
                    all_seqs = list(sequences.values())
                    if len(all_seqs) >= 2:
                        h_seqs, l_seqs = [all_seqs[0]], [all_seqs[1]]
                    elif len(all_seqs) == 1:
                        h_seqs, l_seqs = [all_seqs[0]], []

                if folder == 'Whole_mAb':
                    if not h_seqs or not l_seqs:
                        continue
                    full_sequence = "".join(h_seqs * 2 + l_seqs * 2) 
                    num_h, num_l = 2, 2
                else:
                    full_sequence = ''.join(h_seqs) + ''.join(l_seqs)
                    num_h, num_l = len(h_seqs), len(l_seqs)
                
                features = analyze_sequence_features(full_sequence)
                if features:
                    data = {
                        'antibody': py_file.stem,
                        'folder': folder,
                        'heavy_chains': num_h,
                        'light_chains': num_l,
                        'total_length': len(full_sequence)
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