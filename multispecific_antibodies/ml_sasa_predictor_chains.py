"""
Calculate sequence features for antibody files using BioPython using Pandas for output management.

This script traverses specific subdirectories within a determined root path, parses Python files 
within those directories to extract FASTA-formatted sequence strings, calculates various 
biophysical features (pI, gravy, instability, aromaticity, AA fractions) using Bio.SeqUtils.ProtParam,
and compiles the results into a structured Pandas DataFrame, with one row per antibody product.

REQUIRES:
    pip install biopython pandas xlsxwriter openpyxl
"""

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pathlib import Path
import warnings
import re
import os
from Bio.SeqUtils import ProtParamData # For KD scale access
from openpyxl import load_workbook

# Display formatting
pd.options.display.float_format = '{:.2f}'.format

# --- Dynamic Path Resolution ---

CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]

env = os.environ.get("ROOT_DIR")
if env:
    ROOT_DIR = Path(env).expanduser().resolve()
else:
    ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), None)

if ROOT_DIR is None:
    raise FileNotFoundError("No valid ROOT_DIR found. Check paths or environment variable.")

OUTPUT_FILE_PATH = ROOT_DIR / "all_antibody_sasa_chains.xlsx"

# UPDATED EXCLUSION LIST
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

# Define the 20 standard amino acids for consistency
STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

# Tell Python to ignore common warning messages so the output is cleaner
warnings.filterwarnings('ignore')

# --- Utility Functions ---
def parse_fasta_string(fasta_string, source_file_name):
    """
    Parses a long string of text in FASTA format and cuts it into a list of 
    individual sequence dictionaries (Name, Sequence, Source File).
    """
    sequences_list = []
    # Pattern to find header and sequence blocks
    pattern = re.compile(r'>([^\n]+)\n([A-Z\n\s]+)')
    
    for match in pattern.finditer(fasta_string):
        # Remove newlines and spaces from the sequence itself
        clean_seq = re.sub(r'\s+', '', match.group(2))
        sequences_list.append({
            'Sequence_Name': match.group(1).strip(),
            'Sequence': clean_seq,
            'Source_File': source_file_name
        })
    return sequences_list

def calculate_sasa_features(name, sequence, source_file, folder_name):
    """
    Calculates core sequence-based hydrophobicity/SASA input features for a single, full antibody sequence.
    """
    # Ensure BioPython only gets standard amino acids
    clean_seq = "".join([aa for aa in sequence.upper() if aa in STANDARD_AAS])
    if not clean_seq: return None

    X = ProteinAnalysis(clean_seq)

    # Define which single-letter codes are "Polar" (hydrophilic)
    polar_residues = ['Q', 'N', 'S', 'T', 'H', 'K', 'R', 'E', 'D', 'Y', 'C']
    polar_percentage = sum(clean_seq.count(res) for res in polar_residues) / len(clean_seq)
    
    composition = X.amino_acids_percent 
    
    # --- Hydrophobicity Feature Calculations ---
    kd_scale = ProtParamData.kd
    scores = [kd_scale.get(aa, 0.0) for aa in clean_seq]
    profile = pd.Series(scores).rolling(window=9, min_periods=1, center=True).mean()
    
    avg_kd_hydrophobicity = profile.mean()
    max_kd_hydrophobicity = profile.max()

    gravy_score = X.gravy()
        
    features = {
        'Antibody_Name': name,
        'Format_Folder': folder_name,
        'Source_File': source_file,
        'Sequence_Length_Total_AA': len(clean_seq),
        'pI': X.isoelectric_point(),
        'MW': X.molecular_weight(),
        'Polar_Pct': polar_percentage,
        'GRAVY_Score': gravy_score,
        'Avg_KD_Hydrophobicity': avg_kd_hydrophobicity,
        'Max_KD_Hydrophobicity': max_kd_hydrophobicity
    }
    
    # Add individual percentages for all standard AAs
    for aa in STANDARD_AAS:
         features[f'AA_{aa}_Pct'] = composition.get(aa, 0) 
    return features

# --- Main Processing Loop: Reads Files as Text (Avoids Syntax Errors) ---
def process_all_py_files(root_directory: Path):
    """Aggregates sequences on 1H/1L scale and prepares final data."""
    py_files = list(root_directory.rglob("*.py"))
    antibody_files = [f for f in py_files if f.parent.resolve() != root_directory.resolve() and f.name not in EXCLUDE_FILES]
    
    if not antibody_files:
        print(f"ERROR: No antibody definition files found.")
        return pd.DataFrame() 

    print(f"Found {len(antibody_files)} potential antibody files. Starting sequence aggregation...")
    
    all_antibody_products = [] 
    
    for file_path in antibody_files:
        try:
            # SAFE READ: Open file as text to avoid "Invalid Syntax" on Line 2
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Find the string inside the triple quotes
            fasta_match = re.search(r'["\']{3}(.*?)["\']{3}', content, re.DOTALL)
            if not fasta_match:
                continue
            
            fasta_string_found = fasta_match.group(1).strip()
            if not fasta_string_found.startswith('>'):
                continue
            
            # Parse the FASTA string into individual chains
            chains = parse_fasta_string(fasta_string_found, file_path.name)
            
            # --- ANTIBODY ASSEMBLY LOGIC ---
            heavy_seqs = [c['Sequence'] for c in chains if any(x in c['Sequence_Name'].lower() for x in ['heavy', '_h'])]
            light_seqs = [c['Sequence'] for c in chains if any(x in c['Sequence_Name'].lower() for x in ['light', '_l'])]

            if heavy_seqs and light_seqs:
                full_sequence = heavy_seqs[0] + light_seqs[0]
                units_count = "1H/1L (Fv Match)"
            else:
                full_sequence = "".join([c['Sequence'] for c in chains])
                units_count = f"{len(heavy_seqs)}H/{len(light_seqs)}L"

            if full_sequence:
                all_antibody_products.append({
                    'Name': file_path.stem,
                    'Sequence': full_sequence,
                    'Source_File': file_path.name,
                    'Folder_Name': file_path.parent.name,
                    'Chains_Count': units_count
                })

        except Exception as e:
            print(f"Error reading {file_path.name}: {e}")
    
    final_data = []
    for product in all_antibody_products:
        feat = calculate_sasa_features(
            product['Name'], 
            product['Sequence'], 
            product['Source_File'],
            product['Folder_Name']
        )
        if feat:
            feat['Chains_Count'] = product['Chains_Count']
            final_data.append(feat)

    # --- Export to Excel ---
    df = pd.DataFrame(final_data)
    if not df.empty:
        cols = ['Antibody_Name', 'Format_Folder', 'Chains_Count', 'Source_File', 'Sequence_Length_Total_AA', 'MW', 'pI', 'GRAVY_Score']
        df = df[cols + [c for c in df.columns if c not in cols]]
        df = df.round(2)

        df.to_excel(OUTPUT_FILE_PATH, index=False, engine='openpyxl')
        
        # Format Excel Widths
        wb = load_workbook(OUTPUT_FILE_PATH)
        ws = wb.active
        for col_cells in ws.columns:
            column = col_cells[0].column_letter
            max_len = max(len(str(cell.value or "")) for cell in col_cells)
            ws.column_dimensions[column].width = max_len + 2
        wb.save(OUTPUT_FILE_PATH)
        
        print(f"SUCCESS: {len(df)} antibodies processed. Saved to {OUTPUT_FILE_PATH.name}")
    return df

if __name__ == "__main__":
    process_all_py_files(ROOT_DIR)