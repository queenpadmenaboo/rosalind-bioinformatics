"""
Calculate sequence features for antibody files using BioPython using Pandas for output management.

This script traverses specific subdirectories within a determined root path, parses Python files 
within those directories to extract FASTA-formatted sequence strings, calculates various 
biophysical features (pI, gravy, instability, aromaticity, AA fractions) using Bio.SeqUtils.ProtParam,
and compiles the results into a structured Pandas DataFrame. Finally, it saves the compiled 
data to an Excel file with specific column formatting.

The script dynamically determines the project root directory by checking predefined candidate 
paths, ensuring portability across multiple computers.

REQUIRES:
    pip install biopython pandas xlsxwriter
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import string
import os

# --- Path Configuration (Dynamic Root Finding for Multi-Computer Use) ---

# List potential root directories for different computers
CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]

# Check for environment variable override first
env = os.environ.get("ROOT_DIR")

if env:
    ROOT_DIR = Path(env).expanduser().resolve()
else:
    # Iterate through candidates and use the first one that exists
    ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), None)

# Ensure a valid root directory was found
if ROOT_DIR is None:
    raise FileNotFoundError(
        "No valid ROOT_DIR found. Set ROOT_DIR environment variable or ensure one of the candidate paths exists."
    )

# --- Project Specific Configuration ---

CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'categorize_simple.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.csv', 'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv'
}

# Define the 20 standard amino acids for consistency
STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

# ============================
# FUNCTIONS
# ============================

def parse_py_file(filepath: Path) -> dict:
    """
    Parses a Python file assumed to contain sequences in FASTA format within string literals.
    Extracts headers and sequences.
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    sequences = {}
    current_header = None
    current_seq = []
    for line in content.split('\n'):
        line = line.strip()
        # Heuristic to find FASTA headers starting with '>'
        if line.startswith('>'):
            if current_header:
                seq = ''.join(current_seq)
                if seq.upper() not in ('NA', 'N/A', ''):
                    sequences[current_header] = seq
            current_header = line[1:]
            current_seq = []
        # Heuristic to identify sequence lines (not comments, imports, assignments)
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            if line.lower() not in ('na', 'n/a'):
                current_seq.append(line)
    
    # Process the last sequence found
    if current_header:
        seq = ''.join(current_seq)
        if seq.upper() not in ('NA', 'N/A', ''):
            sequences[current_header] = seq
    return sequences


def analyze_sequence_features(sequence: str) -> dict:
    """
    Calculates various biophysical features of a clean protein sequence using BioPython.
    """
    # Clean sequence of non-standard AAs that might cause BioPython errors
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not clean_seq:
        return None
        
    analysis = ProteinAnalysis(clean_seq)
    features = {
        'calculated_pI': analysis.isoelectric_point(),
        'gravy': analysis.gravy(),
        'instability_index': analysis.instability_index(),
        'aromaticity': analysis.aromaticity()
    }
    
    # Add fractions for all standard amino acids
    aa_fractions = analysis.amino_acids_percent
    for aa in STANDARD_AAS:
         features[f'AA_{aa}_fraction'] = aa_fractions.get(aa, 0.0) 
         
    return features


def get_all_antibody_files(folder: Path) -> list:
    """
    Gathers all relevant Python files from specific subfolders within the project root, 
    excluding utility scripts.
    """
    py_files = []
    for subfolder in CATEGORY_FOLDERS:
        subfolder_path = folder / subfolder
        if subfolder_path.exists():
            for f in subfolder_path.glob('*.py'):
                if f.is_file() and f.name not in EXCLUDE_FILES:
                    py_files.append(f)
    return py_files


def main():
    print("=" * 70)
    print("ANTIBODY FEATURE CALCULATION")
    print("=" * 70)
    
    # Use the resolved single ROOT_DIR variable here
    py_files = get_all_antibody_files(ROOT_DIR) 
    
    print(f"Found {len(py_files)} antibody files\n")
    results = []
    
    for py_file in py_files:
        antibody_name = py_file.stem
        folder = py_file.parent.name
        try:
            sequences = parse_py_file(py_file)
        except Exception as e:
            print(f"ERROR parsing {antibody_name}: {e}")
            continue
            
        heavy_seqs = []
        light_seqs = []
        for header, seq in sequences.items():
            header_lower = header.lower()
            if 'heavy' in header_lower:
                heavy_seqs.append(seq)
            elif 'light' in header_lower:
                light_seqs.append(seq)

        if folder == 'Whole_mAb':
            # Whole mAbs typically have 2 heavy and 2 light chains in the final product
            if len(heavy_seqs) == 0 or len(light_seqs) == 0:
                print(f"ERROR: {antibody_name} missing chains - heavy: {len(heavy_seqs)}, light: {len(light_seqs)}")
                continue
            # Assume one representative H and one L sequence were found, multiply up for MW calculation
            full_sequence = "".join(heavy_seqs * 2 + light_seqs * 2) 
            num_heavy = 2
            num_light = 2
        else:
            # For fragments, just concatenate whatever sequences are present once
            full_sequence = ''.join(heavy_seqs) + ''.join(light_seqs)
            num_heavy = len(heavy_seqs)
            num_light = len(light_seqs)
            
        total_length = len(full_sequence)
        features = analyze_sequence_features(full_sequence)
        
        if features:
            antibody_data = {
                'antibody': antibody_name,
                'folder': folder,
                'heavy_chains': num_heavy,
                'light_chains': num_light,
                'total_length_aa': total_length,
            }
            antibody_data.update(features)
            results.append(antibody_data)
    
    print("Sample results (first 5):")
    print("-" * 70)
    if results:
        # Display summary columns neatly for the console output
        summary_cols = ['antibody', 'calculated_pI', 'total_length_aa', 'gravy', 'instability_index']
        print(pd.DataFrame(results[:5])[summary_cols].to_string(index=False))
    if len(results) > 5:
        print(f"  ... and {len(results) - 5} more\n")
    
    # Use the resolved single ROOT_DIR variable to define the output path
    output_file_path = ROOT_DIR / 'sequence_features.xlsx'
    
    df = pd.DataFrame(results)

    # Add the sum column: Sum across columns (axis=1) for all AA fraction columns
    aa_fraction_cols = [f'AA_{aa}_fraction' for aa in STANDARD_AAS]
    df['Total_AA_Fraction_Sum'] = df[aa_fraction_cols].sum(axis=1)
    
    # --- CRITICAL FIXES FOR FORMATTING START (Pandas/XlsxWriter interactions) ---

    # 1. Ensure the sum column is treated as a float in Pandas before writing
    df['Total_AA_Fraction_Sum'] = df['Total_AA_Fraction_Sum'].astype(float)

    # Initialize the Excel writer engine
    writer = pd.ExcelWriter(output_file_path, engine='xlsxwriter')
    
    # 2. Write to Excel *without* the global float_format argument in to_excel
    df.to_excel(writer, sheet_name='Antibody Features', index=False)
    
    # Get the xlsxwriter workbook and worksheet objects to apply specific formatting
    workbook  = writer.book
    worksheet = writer.sheets['Antibody Features']
    
    # 3. Define a specific format for the sum column: 3 decimal places
    three_decimal_format = workbook.add_format({'num_format': '0.000'})
    
    # 4. Find the index of the 'Total_AA_Fraction_Sum' column dynamically
    sum_col_index = df.columns.get_loc('Total_AA_Fraction_Sum')
    
    # 5. Apply the 3 decimal format to the entire column (col_start, col_end, width, format)
    worksheet.set_column(sum_col_index, sum_col_index, 20, three_decimal_format)

    # Original code to auto-fit other column widths (applied globally)
    for i, col in enumerate(df.columns):
        # Set all columns to a standard width of 20
        worksheet.set_column(i, i, 20) 
        
    # --- CRITICAL FIXES FOR FORMATTING END ---

    # Finalize and close the file writer
    writer.close()
    
    print(f"Results saved to: {output_file_path.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()
