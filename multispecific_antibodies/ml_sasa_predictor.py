import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pathlib import Path
import warnings
import re
import importlib.util 
from Bio import BiopythonDeprecationWarning
import sys
import os

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

OUTPUT_FILENAME = "all_antibody_sasa_features.csv"

# Suppress the Biopython deprecation warnings
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# --- Utility Functions (Parser, Calculator) ---
def parse_fasta_string(fasta_string, source_file_name):
    """Parses a multi-entry FASTA string into a list of dictionaries."""
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

def calculate_sasa_features(name, sequence, source_file):
    """Calculates core SASA input features for a single sequence."""
    
    X = ProteinAnalysis(sequence)
    polar_residues = ['Q', 'N', 'S', 'T', 'H', 'K', 'R', 'E', 'D']
    polar_count = sum(sequence.count(res) for res in polar_residues)
    polar_percentage = (polar_count / len(sequence)) * 100
    composition = X.amino_acids_percent 
    
    return {
        'Sequence_Name': name,
        'Source_File': source_file,
        'Sequence_Length': len(sequence),
        'pI': X.isoelectric_point(),
        'MW': X.molecular_weight(),
        'Polar_Pct': polar_percentage,
        'A_Pct': composition.get('A', 0) * 100,
        'L_Pct': composition.get('L', 0) * 100,
    }

# --- Main Processing Loop: Reads Python Variables from Files ---
def process_all_py_files(root_directory: Path):
    """Recursively finds all .py files and extracts sequence string variables."""
    
    # Find all .py files recursively from the root
    py_files = list(root_directory.rglob("*.py"))
    
    # CRITICAL FIX: Filter out any file located in the ROOT_DIR itself.
    # We only want files that reside in a SUBDIRECTORY.
    antibody_files = [f for f in py_files if f.parent.resolve() != root_directory.resolve()]
    
    if not antibody_files:
        print(f"ERROR: No antibody definition files found in subdirectories of {root_directory}.")
        return pd.DataFrame() 

    print(f"Found {len(antibody_files)} potential antibody files in subfolders. Starting extraction...")
    
    all_sequences_to_process = []
    
    for file_path in antibody_files:
        module_name = file_path.stem
        
        try:
            # Safely load the Python file as a module
            spec = importlib.util.spec_from_file_location(module_name, file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            sequence_string_found = False
            
            # Iterate through all defined objects to find the sequence string
            for attr_name in dir(module):
                attr = getattr(module, attr_name)
                
                # Check if it's a string and contains a FASTA header
                if isinstance(attr, str) and attr.strip().startswith('>') and '\n>' in attr:
                    
                    # Parse and collect sequences
                    parsed_sequences = parse_fasta_string(attr, file_path.name)
                    all_sequences_to_process.extend(parsed_sequences)
                    sequence_string_found = True
                    break 

            if not sequence_string_found:
                 # This warning is for files in subfolders that don't contain sequences, which is fine
                 print(f"    WARNING: No valid sequence string found in {file_path.name}.")
            
        except Exception as e:
            # Catch errors if the file is truly corrupt or has syntax errors
            print(f"    ERROR: Failed to load/execute {file_path.name}. Error: {e}")
    
    print(f"Total {len(all_sequences_to_process)} sequences extracted for feature calculation.")

    # Process all extracted sequences
    results = []
    for seq_data in all_sequences_to_process:
        try:
            features = calculate_sasa_features(seq_data['Sequence_Name'], seq_data['Sequence'], seq_data['Source_File'])
            results.append(features)
        except Exception as e:
            print(f"    WARNING: Skipping {seq_data['Sequence_Name']}. Calculation error: {e}")
            
    return pd.DataFrame(results)

# --- Execution ---
if __name__ == '__main__':
    
    print("--- Starting ML-SASA Feature Extraction ---")
    
    # 1. Run the main processor
    results_df = process_all_py_files(ROOT_DIR)
    
    if results_df.empty:
        print("Feature calculation finished, but no valid data was generated.")
        sys.exit(0)
    
    # 2. Final data preparation
    cols = ['Sequence_Name', 'Source_File', 'Sequence_Length', 'pI', 'MW', 'Polar_Pct', 'A_Pct', 'L_Pct']
    results_df = results_df.reindex(columns=cols)
    
    # 3. Print a summary
    print(f"\nSUCCESS: Processed {len(results_df)} sequences.")
    print("--- ML-SASA Feature Comparison Table Summary (First 5 Rows) ---")
    print(results_df.head().round(2).to_string())
    
    # 4. Save to CSV file
    results_df.to_csv(OUTPUT_FILENAME, index=False, float_format='%.2f')
    print(f"\nFinal results saved to '{OUTPUT_FILENAME}'")