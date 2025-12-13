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
# Added this line to format how all Pandas DataFrames display floats (2 decimal places)
pd.options.display.float_format = '{:.2f}'.format
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pathlib import Path
import warnings
import re
import importlib.util 
from Bio import BiopythonDeprecationWarning
import sys
import os
import numpy as np
from Bio.SeqUtils import ProtParamData # For KD scale access
from openpyxl import load_workbook 


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
    raise FileNotFoundError(
        "No valid ROOT_DIR found. Set ROOT_DIR environment variable or ensure one of the candidate paths exists."
    )

# Changed output file extension to .xlsx
OUTPUT_FILENAME = "all_antibody_sasa_chains.xlsx" 
OUTPUT_FILE_PATH = ROOT_DIR / OUTPUT_FILENAME

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.csv', 'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py'
}


# Define the 20 standard amino acids for consistency
STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

# Tell Python to ignore common warning messages so the output is cleaner
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# --- Utility Functions ---
def parse_fasta_string(fasta_string, source_file_name):
    """
    Parses a long string of text in FASTA format and cuts it into a list of 
    individual sequence dictionaries (Name, Sequence, Source File).
    """
    sequences_list = []
    fasta_string = fasta_string.strip('\"\'') 
    pattern = re.compile(r'>([^\n]+)\n([A-Z\n]+)')
    
    for match in pattern.finditer(fasta_string):
        name = match.group(1).strip()
        sequence = match.group(2).replace('\n', '').strip()
        
        sequences_list.append({
            'Sequence_Name': name,
            'Sequence': sequence,
            'Source_File': source_file_name
        })
    return sequences_list

def calculate_sasa_features(name, sequence, source_file, folder_name):
    """
    Calculates core sequence-based hydrophobicity/SASA input features for a single, full antibody sequence.
    """
    X = ProteinAnalysis(sequence)

    # Define which single-letter codes are "Polar" (hydrophilic)
    polar_residues = ['Q', 'N', 'S', 'T', 'H', 'K', 'R', 'E', 'D', 'Y', 'C']
    polar_count = sum(sequence.count(res) for res in polar_residues)
    
    polar_percentage = (polar_count / len(sequence)) 
    composition = X.amino_acids_percent 
    
    # --- Hydrophobicity Feature Calculations ---
    kd_scale = ProtParamData.kd
    scores = [kd_scale.get(aa, 0.0) for aa in sequence]
    profile = pd.Series(scores).rolling(window=9, min_periods=1, center=True).mean()
    
    avg_kd_hydrophobicity = profile.mean()
    max_kd_hydrophobicity = profile.max()

    gravy_score = X.gravy()
        
    features = {
        'Antibody_Name': name,
        'Format_Folder': folder_name,
        'Source_File': source_file,
        'Sequence_Length_Total_AA': len(sequence),
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


# --- Main Processing Loop: Reads Python Variables from Files ---
def process_all_py_files(root_directory: Path):
    """
    Scans your specified folder structure, finds the data files, 
    aggregates sequences by antibody product, and calculates features per product.
    """
    
    py_files = list(root_directory.rglob("*.py"))
    
    # Filter out files in the main folder itself, focusing only on subfolders
    antibody_files = [f for f in py_files if f.parent.resolve() != root_directory.resolve() and f.name not in EXCLUDE_FILES]
    
    if not antibody_files:
        print(f"ERROR: No antibody definition files found in subdirectories of {root_directory}.")
        return pd.DataFrame() 

    print(f"Found {len(antibody_files)} potential antibody files in subfolders. Starting sequence aggregation...")
    
    all_antibody_products = [] # This will store the single combined sequences
    
    for file_path in antibody_files:
        module_name = file_path.stem
        folder_name = file_path.parent.name # e.g., 'Whole_mAb'
        antibody_name = file_path.stem      # e.g., 'acimtamig'

        try:
            # Safely load the Python file as a module
            spec = importlib.util.spec_from_file_location(module_name, file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            fasta_string_found = None
            
            # Iterate through all defined objects to find the sequence string variable
            for attr_name in dir(module):
                attr = getattr(module, attr_name)
                if isinstance(attr, str) and attr.strip().startswith('>') and '\n>' in attr:
                    fasta_string_found = attr
                    break 

            if not fasta_string_found:
                 print(f"    WARNING: No valid sequence string found in {file_path.name}.")
                 continue
            
            # Parse the FASTA string into individual chains
            chains = parse_fasta_string(fasta_string_found, file_path.name)
            
            # --- AGGREGATION LOGIC (UPDATED FILTER) ---
            # Updated these lines to include checks for "_h" and "_l" for H1/L1 conventions
            heavy_seqs = [c['Sequence'] for c in chains if 'heavy' in c['Sequence_Name'].lower() or '_h' in c['Sequence_Name'].lower()]
            light_seqs = [c['Sequence'] for c in chains if 'light' in c['Sequence_Name'].lower() or '_l' in c['Sequence_Name'].lower()]

            if folder_name == 'Whole_mAb':
                # For a whole mAb, use 2 heavy and 2 light chains
                full_sequence = "".join(heavy_seqs * 2 + light_seqs * 2) 
            else:
                # For all other formats, simply concatenate the chains found once
                full_sequence = ''.join(heavy_seqs) + ''.join(light_seqs)
                
            if full_sequence:
                all_antibody_products.append({
                    'Name': antibody_name,
                    'Sequence': full_sequence,
                    'Source_File': file_path.name,
                    'Folder_Name': folder_name
                })
            else:
                # This message will ideally no longer appear for emicizumab.py
                print(f"    SKIPPED: {file_path.name} was skipped because no valid sequence was generated during aggregation.")
            # --- END AGGREGATION LOGIC ---

        except Exception as e:
            print(f"    ERROR: Failed to load/execute {file_path.name}. Error: {e}")
    
    print(f"Total {len(all_antibody_products)} full antibody products extracted for feature calculation.")

    # Process all aggregated sequences
    results = []
    for product_data in all_antibody_products:
        try:
            features = calculate_sasa_features(
                product_data['Name'], 
                product_data['Sequence'], 
                product_data['Source_File'],
                product_data['Folder_Name']
            )
            results.append(features)
        except Exception as e:
            print(f"    ERROR: Failed to calculate features for {product_data['Name']}. Error: {e}")

    return pd.DataFrame(results)


# --- Execution ---
if __name__ == "__main__":
    # Ensure pandas displays all columns in the console output summary
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)

    if ROOT_DIR:
        print(f"--- Starting Feature Calculation in: {ROOT_DIR} ---")
        df_features = process_all_py_files(ROOT_DIR)

        if not df_features.empty:
            # Round the *entire* DataFrame to 2 decimal places before saving/printing
            df_features = df_features.round(2)

            print(f"\nSUCCESS: Processed {len(df_features)} antibodies.")
            print("\n" + "=" * 60)
            print("ML-SASA Feature Comparison Table Summary (First 5 Rows)")
            print("=" * 60)
            print(df_features.head())
            print("\n")
            
            # --- FINAL EXCEL WRITING LOGIC (With Filters, Widths, and Freeze) ---
            # Save the results to the determined path as a normal Excel file first
            df_features.to_excel(OUTPUT_FILE_PATH, index=False, engine='openpyxl')

            # Use openpyxl to apply formatting after pandas writes the data
            wb = load_workbook(OUTPUT_FILE_PATH)
            ws = wb.active
            
            # FREEZE THE TOP ROW (FREEZE PANES ABOVE CELL A2)
            ws.freeze_panes = 'A2'

            # Apply Auto-filter to the header row (row index 1 in Excel)
            # FIX applied here: Accessing shape elements by index correctly (0 for rows, 1 for columns)
            max_row = df_features.shape[0] + 1
            max_col = df_features.shape[1]
            filter_range = f"A1:{ws.cell(row=max_row, column=max_col).coordinate}"
            ws.auto_filter.ref = filter_range

            # Automatically adjust column widths
            for col_cells in ws.columns:
                max_length = 0
                column = col_cells[0].column_letter # FIX: Access letter from the first cell in the tuple
                for cell in col_cells:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = (max_length + 2)
                ws.column_dimensions[column].width = adjusted_width
            
            wb.save(OUTPUT_FILE_PATH)
            # --- END FINAL EXCEL WRITING LOGIC ---

            print(f"Final results saved to '{OUTPUT_FILE_PATH.name}' with filters, frozen header, and optimized column widths.")
        else:
            print("\nNo features calculated. DataFrame is empty.")
    else:
        print("Script failed to find a valid root directory for processing.")
