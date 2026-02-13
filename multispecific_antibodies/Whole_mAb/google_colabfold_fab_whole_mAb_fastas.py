# ============================================================================
# COLABFOLD FAB BUILDER - BATCH MODE (temp_fastas only)
# ============================================================================

import os
import pandas as pd
import re
import warnings
import logging
from pathlib import Path
from tqdm import tqdm
import subprocess
import json
from google.colab import files

# Suppress library noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# --- COLAB PATHS ---
BASE_DIR = Path("/content/antibody_data")
OUTPUT_ROOT = Path("/content/PDB_Output_ColabFold_Fab_Structures")
CSV_PATH = Path("/content/TheraSAbDab_SeqStruc_07Dec2025.csv")
TEMP_FASTA_DIR = Path("/content/temp_fastas")
COLABFOLD_BIN = "/usr/local/bin/colabfold_batch"

# IUPAC Standard 20 Amino Acids
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# ============================================================================
# CONSTANTS — CH1 ONLY + CL ONLY (Fab level modeling)
# ============================================================================

# Human Constant Region Sequences (Standard Reference Library)
CONSTANTS = {
    'HEAVY_CH1': {
        'G1': (
            "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
        ),
        'G2': (
            "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSNFGTQTYTCNVDHKPSNTKVDKTV"
        ),
        'G4': (
            "ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSSLGTKTYTCNVDHKPSNTKVDKRV"
        ),
    },

    'LIGHT_CL': {
        'Kappa': (
            "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
        ),
        'Lambda': (
            "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
        ),
    }
}

def clean_sequence(seq):
    """Removes any characters not in the standard 20 AA alphabet."""
    return "".join([aa for aa in seq.upper() if aa in AA_ALPHABET])

def detect_isotype(file_path, lookup_dict):
    """Scans the source file for isotype keywords."""
    stem_lower = file_path.stem.lower()
    # Priority 1: CSV Lookup
    if stem_lower in lookup_dict:
        h_raw, l_raw = lookup_dict[stem_lower]
        h_iso = 'G1'
        if 'G2' in str(h_raw): h_iso = 'G2'
        elif 'G4' in str(h_raw): h_iso = 'G4'
        l_iso = 'Lambda' if 'LAMBDA' in str(l_raw).upper() else 'Kappa'
        return h_iso, l_iso

    # Priority 2: Keyword Scan
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read().upper()
        h_iso = 'G1'
        if 'G2' in content or 'IGG2' in content: h_iso = 'G2'
        elif 'G4' in content or 'IGG4' in content: h_iso = 'G4'
        l_iso = 'Lambda' if 'LAMBDA' in content else 'Kappa'
        return h_iso, l_iso
    except:
        return 'G1', 'Kappa'

# EXTRACT SEQUENCES and return ONLY ONE Fab pair (no duplication)
def extract_chains_dynamic(file_path):
    found_sequences = []

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()

        # Handle Python triple-quoted strings: remove everything before first FASTA entry
        # This handles: antibody_name = """   or   antibody_name = '''
        if '"""' in content or "'''" in content:
            # Find first > after any triple quotes
            first_fasta = content.find('>')
            if first_fasta > 0:
                content = content[first_fasta:]

        # Split into blocks based on the FASTA-style headers
        blocks = content.split('>')
        for block in blocks[1:]:  # Skip any text before the first '>'
            lines = block.strip().split('\n')
            if len(lines) < 2:
                continue

            header = lines[0].upper()
            # Combine all lines after the header line
            raw_data = "".join(lines[1:]).strip()

            # Remove Python artifacts (quotes, commas, semicolons, spaces, newlines)
            # Keep only uppercase letters
            clean_seq = re.sub(r'[^A-Z]', '', raw_data.upper())

            # Validate sequence quality
            if clean_seq and len(clean_seq) > 20:
                # Validate it's actually amino acids (not random Python code)
                invalid = set(clean_seq) - set(AA_ALPHABET)
                if invalid:
                    continue

                # Check for heavy chain indicators
                if any(keyword in header for keyword in ['HEAVY', 'HC', 'VH', 'H_CHAIN', '_H_', '_HC_']):
                    chain_type = 'H'
                # Check for light chain indicators
                elif any(keyword in header for keyword in ['LIGHT', 'LC', 'VL', 'L_CHAIN', 'KAPPA', 'LAMBDA', '_L_', '_LC_']):
                    chain_type = 'L'
                # Fallback: Default to light chain if header is ambiguous
                else:
                    chain_type = 'L'

                found_sequences.append((chain_type, clean_seq))

    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

    heavies = [s for t, s in found_sequences if t == 'H']
    lights = [s for t, s in found_sequences if t == 'L']

    # Validate found chains
    if not heavies or not lights:
        print(f"Warning: {file_path.stem} missing chains")
        return None

    heavy_fv = heavies[0]
    light_fv = lights[0]

    return heavy_fv, light_fv

# CREATE EXACTLY 2 CHAINS (Fab)
def create_fab_fasta(heavy_fv: str, light_fv: str, h_iso: str, l_iso: str, output_fasta: Path):

    heavy_fab = heavy_fv + CONSTANTS['HEAVY_CH1'][h_iso]
    light_fab = light_fv + CONSTANTS['LIGHT_CL'][l_iso]

    with open(output_fasta, 'w') as f:
        f.write(">Fab\n")
        f.write(f"{heavy_fab}:{light_fab}\n")

def run_pipeline():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    os.makedirs(TEMP_FASTA_DIR, exist_ok=True)

    # Load Master Isotype Data from CSV
    isotype_lookup = {}
    if CSV_PATH.exists():
        try:
            df = pd.read_csv(CSV_PATH)
            isotype_lookup = dict(zip(df['Therapeutic'].str.lower(), zip(df['CH1 Isotype'], df['VD LC'])))
            print(f"--- LOADED {len(isotype_lookup)} ISOTYPES FROM CSV ---")
        except Exception as e:
            print(f"--- ERROR LOADING CSV: {e}. Falling back to keyword search only. ---")
    else:
        print("--- NO CSV FOUND. Using keyword search for isotypes. ---")

    # Get antibody files   
    antibody_files = sorted([
        f for f in BASE_DIR.rglob("*.py") 
        if not f.name.startswith("._")
    ])
    
    if not antibody_files:
        print("[ERROR] No antibody files found")
        return

    print(f"\nSTARTING COLABFOLD ASSEMBLY: {len(antibody_files)} antibody file(s)\n")

    for f_path in tqdm(antibody_files, desc="Building Fab Structures with ColabFold"):
        tqdm.write(f"\n{'='*60}")
        tqdm.write(f"FOLDING: {f_path.stem}")
        tqdm.write(f"{'='*60}")
        
        # Extract sequences
        result = extract_chains_dynamic(f_path)
        if not result:
            tqdm.write(f"[ERROR] No valid sequences found in {f_path.stem}")
            continue

        heavy_fv, light_fv = result
        
        # Detect isotype
        h_iso, l_iso = detect_isotype(f_path, isotype_lookup)
        tqdm.write(f"Isotype: Heavy={h_iso}, Light={l_iso}")
        tqdm.write(f"Heavy chain: {len(heavy_fv)} aa")
        tqdm.write(f"Light chain: {len(light_fv)} aa")
        
        # Create Fab FASTA
        fasta_path = TEMP_FASTA_DIR / f"{f_path.stem}.fasta"
        create_fab_fasta(heavy_fv, light_fv, h_iso, l_iso, fasta_path)
        tqdm.write(f"Created FASTA: {fasta_path}")
        
        # Create antibody-specific output directory
        antibody_output = OUTPUT_ROOT / f_path.stem
        os.makedirs(antibody_output, exist_ok=True)  

# ============================================================================
# STEP 1: UPLOAD FILES
# ============================================================================
print("="*60)
print("STEP 1: UPLOAD YOUR FILES")
print("="*60)
print("\n1. Upload antibody.py files")
print("2. Upload your CSV file (TheraSAbDab_SeqStruc_07Dec2025.csv)")
print("\nClick 'Choose Files' below:\n")

uploaded = files.upload()

# Create directories and move files
os.makedirs(BASE_DIR, exist_ok=True)

for filename, content in uploaded.items():
    if filename.endswith('.py'):
        # Save antibody file
        dest = BASE_DIR / filename
        with open(dest, 'wb') as f:
            f.write(content)
        print(f"✓ Saved antibody file: {dest}")
    elif filename.endswith('.csv'):
        # Save CSV
        with open(CSV_PATH, 'wb') as f:
            f.write(content)
        print(f"✓ Saved CSV: {CSV_PATH}")

print("\n" + "="*60)
print("FILES UPLOADED SUCCESSFULLY")
print("="*60)

# ============================================================================
# STEP 2: RUN PIPELINE
# ============================================================================
print("\n\nStarting pipeline in 3 seconds...\n")
import time
time.sleep(3)

if __name__ == "__main__":
    run_pipeline()