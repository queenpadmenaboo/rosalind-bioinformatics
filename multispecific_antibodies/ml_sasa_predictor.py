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

from Bio.SeqUtils import ProtParamData                      # Contains the actual numeric values for hydrophobicity scales
import numpy as np                                          # A powerful library for math operations and averages


ROOT_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies") 

OUTPUT_FILENAME = ROOT_DIR / "all_antibody_sasa_features.csv"

OUTPUT_FILENAME_STR = str(OUTPUT_FILENAME)


EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.csv', 'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py'
}

# Tell Python to ignore common warning messages so the output is cleaner
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# --- Utility Functions ---
def parse_fasta_string(fasta_string, source_file_name):
    """
    Parses a long string of text in FASTA format and cuts it into a list of 
    individual sequence dictionaries (Name, Sequence, Source File).
    
    Args:
        fasta_string (str): Raw FASTA format text containing one or more sequences 
                            from a .py file variable.
        source_file_name (str): The name of the file this sequence came from.

    Returns:
        list of dicts: E.g., [{'Sequence_Name': 'Acimtamig_Heavy_Chain_1', 'Sequence': 'QVQLVQSGAEVKKPG...', 'Source_File': 'acimtamig.py'}]
    """
    sequences_list = []
    fasta_string = fasta_string.strip('\"\'') 
    # The 'pattern' variable looks for a '>' character followed by text (the name) 
    # and then a lot of A-Z letters (the sequence).
    pattern = re.compile(r'>([^\n]+)\n([A-Z\n]+)')
    
    # Loop through all matches found in the big text block
    for match in pattern.finditer(fasta_string):
        name = match.group(1).strip()
        sequence = match.group(2).replace('\n', '').strip()
        
        # Add the structured data to our list
        sequences_list.append({
            'Sequence_Name': name,
            'Sequence': sequence,
            'Source_File': source_file_name
        })
    return sequences_list

def calculate_sasa_features(name, sequence, source_file):
    """
    Calculates core sequence-based hydrophobicity/SASA input features for a single sequence.
    
    NOTE: These calculations use sequence rules, NOT a 3D model of the protein.

    Args:
        name (str): The name of the sequence (e.g., 'Acimtamig_Heavy_Chain_1').
        sequence (str): The raw amino acid sequence (e.g., 'QVQLVQSGAEVKKPG...').
        source_file (str): The name of the file it originated from.

    Returns:
        dict: A dictionary of calculated features ready for a Pandas DataFrame.
              E.g., {'Sequence_Length': 120, 'pI': 8.25, 'GRAVY_Score': -0.15, ...}
    """
    # Create a ProteinAnalysis object from the sequence; this object holds many useful tools
    X = ProteinAnalysis(sequence)

    # Define which single-letter codes are "Polar" (hydrophilic)
    polar_residues = ['Q', 'N', 'S', 'T', 'H', 'K', 'R', 'E', 'D', 'Y', 'C']
    polar_count = sum(sequence.count(res) for res in polar_residues)
    
    polar_percentage = (polar_count / len(sequence)) 
    # Get the percentage breakdown of every single amino acid type (A, L, G, etc.)
    composition = X.amino_acids_percent 
    
    # --- New Hydrophobicity Feature Calculations ---
    
    # Load the Kyte-Doolittle scale (a dictionary mapping letters to a slipperiness number)
    kd_scale = ProtParamData.kd
    # Create a list of 'slipperiness scores' for our specific sequence
    scores = [kd_scale.get(aa, 0.0) for aa in sequence]
    # Use Pandas to calculate the 'rolling average' over a window of 9 residues. 
    # This helps find local hydrophobic "patches".
    profile = pd.Series(scores).rolling(window=9, min_periods=1, center=True).mean()
    
    avg_kd_hydrophobicity = profile.mean()
    max_kd_hydrophobicity = profile.max()

    # Calculate the GRAVY score (Grand Average of Hydropathicity)
    gravy_score = X.gravy()
        
    return {
        'Sequence_Name': name,
        'Source_File': source_file,
        'Sequence_Length': len(sequence),
        'pI': X.isoelectric_point(),
        'MW': X.molecular_weight(),
        'Polar_Pct': polar_percentage,
        # Get individual percentages for all nonpolar residues (used for hydrophobicity learning)
        'A_Pct': composition.get('A', 0),
        'F_Pct': composition.get('F', 0),
        'G_Pct': composition.get('G', 0),
        'I_Pct': composition.get('I', 0),
        'L_Pct': composition.get('L', 0),
        'M_Pct': composition.get('M', 0),
        'P_Pct': composition.get('P', 0),
        'V_Pct': composition.get('V', 0),  
        'W_Pct': composition.get('W', 0),
        'GRAVY_Score': gravy_score,
        'Avg_KD_Hydrophobicity': avg_kd_hydrophobicity,
        'Max_KD_Hydrophobicity': max_kd_hydrophobicity
    }

# --- Main Processing Loop: Reads Python Variables from Files ---
def process_all_py_files(root_directory: Path):
    """
    Scans your specified folder structure, finds the data files, 
    and organizes all sequences for processing.
    
    Args:
        root_directory (Path): The starting folder path to scan for .py files.

    Returns:
        pd.DataFrame: A data table (DataFrame) containing all calculated features 
                      for every single sequence found.
    """
    
    # Find all files ending in .py starting from the root folder
    py_files = list(root_directory.rglob("*.py"))
    
    # Filter out files in the main folder itself, focusing only on subfolders
    antibody_files = [f for f in py_files if f.parent.resolve() != root_directory.resolve()]
    
    if not antibody_files:
        print(f"ERROR: No antibody definition files found in subdirectories of {root_directory}.")
        return pd.DataFrame() # Return an empty table if nothing is found

    print(f"Found {len(antibody_files)} potential antibody files in subfolders. Starting extraction...")
    
    all_sequences_to_process = []
    
    for file_path in antibody_files:
        module_name = file_path.stem
        
        try:
            # We treat the .py file like a Python module to access the variables inside it
            spec = importlib.util.spec_from_file_location(module_name, file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            
            sequence_string_found = False
            
            # Look at every variable defined in the file we just loaded
            for attr_name in dir(module):
                attr = getattr(module, attr_name)
                
                # Check if the variable is a string AND looks like our FASTA format
                if isinstance(attr, str) and attr.strip().startswith('>') and '\n>' in attr:
                    
                    # Send the large sequence string to our parser function
                    parsed_sequences = parse_fasta_string(attr, file_path.name)
                    all_sequences_to_process.extend(parsed_sequences)
                    sequence_string_found = True
                    break 

            if not sequence_string_found:
                 print(f"    WARNING: No valid sequence string found in {file_path.name}.")
            
        except Exception as e:
            # If the .py file itself has a syntax error, this catches it
            print(f"    ERROR: Failed to load/execute {file_path.name}. Error: {e}")
    
    print(f"Total {len(all_sequences_to_process)} sequences extracted for feature calculation.")

    # Process all extracted sequences one by one
    results = []
    for seq_data in all_sequences_to_process:
        try:
            # Run our calculation function defined earlier
            features = calculate_sasa_features(seq_data['Sequence_Name'], seq_data['Sequence'], seq_data['Source_File'])
            results.append(features)
        except Exception as e:
            print(f"    WARNING: Skipping {seq_data['Sequence_Name']}. Calculation error: {e}")
            
    # Turn the final list of results into a pandas DataFrame (our data table)
    return pd.DataFrame(results)

# --- Execution (The main part that runs when you start the script) ---
if __name__ == '__main__':
    
    print("--- Starting ML-SASA Feature Extraction ---")
    
    # 1. Run the main processing function
    results_df = process_all_py_files(ROOT_DIR)
    
    if results_df.empty:
        print("Feature calculation finished, but no valid data was generated.")
        sys.exit(0) # Stop the script here if there is no data
    
    # 2. Final data preparation: Define the exact order of columns for the final file
    cols = ['Sequence_Name', 'Source_File', 'Sequence_Length', 'pI', 'MW', 'Polar_Pct', 
            'A_Pct', 'L_Pct', 'V_Pct', 'I_Pct', 'P_Pct', 'F_Pct', 'W_Pct', 'M_Pct', 'G_Pct',
            'GRAVY_Score', 'Avg_KD_Hydrophobicity', 'Max_KD_Hydrophobicity']
    # Reorganize the DataFrame columns into the order specified above
    results_df = results_df.reindex(columns=cols)
    
    # 3. Print a summary to the console so we know it worked
    print(f"\nSUCCESS: Processed {len(results_df)} sequences.")
    print("--- ML-SASA Feature Comparison Table Summary (First 5 Rows) ---")
    # Show the first 5 rows nicely formatted
    print(results_df.head().round(2).to_string())
    
    # 4. Save the final data table to a CSV file (Comma Separated Values)
    results_df.to_csv(OUTPUT_FILENAME, index=False, float_format='%.2f')
    print(f"\nFinal results saved to '{OUTPUT_FILENAME}'")
