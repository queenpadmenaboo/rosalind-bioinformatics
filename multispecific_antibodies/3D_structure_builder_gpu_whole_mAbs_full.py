import os
# 1. GPU SECURITY BYPASS (Required for 4080 SUPER / Torch 2.6+)
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
import pandas as pd
import torch
torch.set_float32_matmul_precision('highest')

import re
import warnings
import logging
from pathlib import Path
from tqdm import tqdm

# Suppress library noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# Load IgFold
from igfold import IgFoldRunner

# Chain ID map for assembly: A, B, C, D...
CHAIN_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

# IUPAC Standard 20 Amino Acids (IgFold requirement)
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# --- CONFIG & EXCLUSIONS ---
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")
OUTPUT_ROOT = BASE_DIR / "PDB_Output_Files_GPU_Full"
CSV_PATH = Path(r"C:\Users\bunsr\TheraSAbDab_SeqStruc_07Dec2025.csv")

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py', 'analyze_hotspots.py', 'aggregation_predictor.py', 'Aggregation_Risk_Report.xlsx',
    'antibody_diagnostic_tool.py', 'Antibody_Comparison_Report_2025.xlsx', 'build_pnas_shadow.py',
    'compare_hydrophobicity.py', 'developability_hotspots.xlsx', 'mAb_truth_engine.py',
    'mAb_Truth_Engine_Master.xlsx', 'MISSING_ANTIBODIES_LOG.xlsx', 'pnas_validator.py',
    'PNAS_VS_CALCULATIONS.xlsx', '3D_structure_builder_gpu.py', 'format_pairing_predictor.py',
    'immunogenicity_predictor.py', 'thermal_stability_predictor.py', 'mAb_chains_logic.py',
    'Antibody_Chain_Count_Summary.xlsx', '3D_structure_builder_gpu_whole_mAbs.py',
    '3D_structure_builder_gpu_whole_mAbs_full.py', 'non_human_antibodies_to_remove.txt',
    'whole_mAbs_folder_check.py', 'whole_mAbs_isotypes_check.py', 'whole_mAb_antibody_list.txt'
}

# Human Constant Region Sequences (Standard Reference Library)
CONSTANTS = {
    'HEAVY': {
        # Full CH1-Hinge-CH2-CH3
        'G1': ("ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV"
               "EPKSCDKTHTCPPCP"
               "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAK"
               "GQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"),
        'G2': ("ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSNFGTQTYTCNVDHKPSNTKVDKTV"
               "ERKCCVECPPCP"
               "APPVAGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTFRVVSVLTVVHQDWLNGKEYKCKVSNKGLPAPIEKTISKTK"
               "GQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDISVEWESNGQPENNYKTTPPMLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"),
        'G4': ("ASTKGPSVFPLAPCSRSTSESTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSSLGTKTYTCNVDHKPSNTKVDKRV"
               "ESKYGPPCPSCP"
               "APEFLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSQEDPEVQFNWYVDGVEVHNAKTKPREEQFNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKGLPSSIEKTISKAK"
               "GQPREPQVYTLPPSQEEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSRLTVDKSRWQEGNVFSCSVMHEALHNHYTQKSLSLSLGK")
    },
    'LIGHT': {
        'Kappa': "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
        'Lambda': "GQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPVKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS"
    }
}

def clean_sequence(seq):
    """Removes any characters not in the standard 20 AA alphabet."""
    return "".join([aa for aa in seq.upper() if aa in AA_ALPHABET])

# Function to detect which 'tail' to glue onto the variable 'head'
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

def extract_chains_dynamic(file_path):
    """Extracts H and L sequences and returns them as pairs."""
    found_sequences = []
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
       # Split into blocks based on the FASTA-style headers
        blocks = content.split('>')
        for block in blocks[1:]: # Skip text before the first '>'
            lines = block.strip().split('\n')
            if len(lines) < 2:
                continue
            
            header = lines[0].upper()
            # Combine all lines after the header line
            raw_data = "".join(lines[1:]).strip()
            # Remove any non-alpha characters (quotes, spaces, commas, etc.)
            clean_seq = re.sub(r'[^A-Z]', '', raw_data.upper())
            
            # Validate sequence quality: Must have content, not be "NA", and exceed 20 residues (minimum Fv fragment size)
            if clean_seq and clean_seq != "NA" and len(clean_seq) > 20:
                """
                ROBUST CHAIN TYPE DETECTION
                ---------------------------
                Problem: Antibody files use inconsistent naming conventions:
                - "Heavy_Chain" vs "HC" vs "VH" vs "H_1"
                - "Light_Chain" vs "LC" vs "VL" vs "Kappa" vs "Lambda"
                
                Solution: Check multiple keyword patterns after uppercase conversion (Line 67)
                
                Heavy Chain Keywords: HEAVY, HC, VH, H_CHAIN, _H_, _HC_
                Light Chain Keywords: LIGHT, LC, VL, L_CHAIN, KAPPA, LAMBDA, _L_, _LC_
                
                Default: If no match, assume Light chain (conservative fallback)
                """
    
                # Check for heavy chain indicators (multiple patterns)
                if any(keyword in header for keyword in ['HEAVY', 'HC', 'VH', 'H_CHAIN', '_H_', '_HC_']):
                    chain_type = 'H'
                
                # Check for light chain indicators (multiple patterns including isotypes)
                elif any(keyword in header for keyword in ['LIGHT', 'LC', 'VL', 'L_CHAIN', 'KAPPA', 'LAMBDA', '_L_', '_LC_']):
                    chain_type = 'L'
                
                # Fallback: Default to light chain if header is ambiguous
                else:
                    chain_type = 'L'
                
                # Store tuple of (chain_type, sequence) for downstream H-L pairing
                found_sequences.append((chain_type, clean_seq))
                          
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    
    heavies = [s for t, s in found_sequences if t == 'H']
    lights = [s for t, s in found_sequences if t == 'L']
    pairs = list(zip(heavies, lights))


    # If Whole_mAb has only 1 pair, duplicate it for 2H:2L biological assembly
    if len(pairs) == 1:
        pairs = pairs * 2  # Creates [(VH, VL), (VH, VL)]
    return pairs if pairs else None
    

def run_pipeline():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    runner = IgFoldRunner()

    isotype_lookup = {}
    # Load Master Isotype Data from CSV
    try:
        df = pd.read_csv(CSV_PATH)
        isotype_lookup = dict(zip(df['Therapeutic'].str.lower(), zip(df['CH1 Isotype'], df['VD LC'])))
        print(f"--- LOADED {len(isotype_lookup)} ISOTYPES FROM CSV ---")
    except Exception as e:
        print(f"--- ERROR LOADING CSV: {e}. Falling back to keyword search only. ---")
       
    antibody_files = [
        f for f in BASE_DIR.rglob("*.py") 
        if f.name not in EXCLUDE_FILES 
        and "PDB_Output_Files_GPU" not in f.parts
        and "PDB_Output_Files_GPU_Full" not in f.parts
        and not f.name.startswith("._")
    ]
    
    print(f"\nSUCCESS: STARTING ASSEMBLY: {len(antibody_files)} antibody files in Whole_mAb.")

    for f_path in tqdm(antibody_files, desc="Building 3D Structures of Whole mAbs"):
        tqdm.write(f"--- FOLDING: {f_path.stem} ---")
        target_dir = OUTPUT_ROOT
        pdb_path = target_dir / f"{f_path.stem}.pdb"

        if pdb_path.exists():       
            pdb_path.unlink()

        # Get list of Whole mAb sequences
        pairs = extract_chains_dynamic(f_path)
        if not pairs:
            continue
                    
        # Detect isotype once per file to select tails
        h_iso, l_iso = detect_isotype(f_path, isotype_lookup)   
        
        # Loop through each sequence, fold individually, and assemble
        for i, (heavy_fv, light_fv) in enumerate(pairs):
            temp_pdb = f"temp_{f_path.stem}_pair_{i}.pdb"
            # Glue the constant tails to the variable heads
            full_heavy = heavy_fv + CONSTANTS['HEAVY'][h_iso]
            full_light = light_fv + CONSTANTS['LIGHT'][l_iso]   

            try:
                # 1. Fold the full-length chains
                runner.fold(
                    pdb_file=temp_pdb,
                    sequences={"H": full_heavy, "L": full_light}, 
                    do_refine=False,    # Refinement on full mAb is unstable/slow
                    do_renum=False
                )

                # 2. Append to master PDB with 150A Offset and correct Chain ID
                if Path(temp_pdb).exists():
                    with open(temp_pdb, 'r') as f_in, open(pdb_path, 'a') as f_out:
                        for line in f_in:
                            if line.startswith("ATOM"):
                                # Get original chain ID ('H' or 'L')
                                original_chain = line[21].strip()
        
                                # Assign new chain ID for Whole_mAb (A,B,C,D pattern)
                                if original_chain == 'H':
                                    new_chain_id = CHAIN_LABELS[i * 2]      # Pair 0 -> A, Pair 1 -> C
                                else:  # original_chain == 'L'
                                    new_chain_id = CHAIN_LABELS[i * 2 + 1]  # Pair 0 -> B, Pair 1 -> D
                            
                                # Apply 150A offset for spatial separation
                                x_orig = float(line[30:38])
                                x_new = x_orig + (i * 150.0)
                            
                                # Write: [cols 1-21] + new_chain_ID + [cols 23-30] + shifted_X + [cols 39-end]
                                f_out.write(line[:21] + new_chain_id + line[22:30] + f"{x_new:8.3f}" + line[38:])

                        f_out.write("TER\n")
                        
            except Exception as e:
                tqdm.write(f"\n[ERROR] folding {f_path.stem}: {e}")
            finally:
                if Path(temp_pdb).exists(): Path(temp_pdb).unlink()
                # Clean up IgFold's temporary fasta files
                tf = temp_pdb.replace(".pdb", ".fasta")
                if Path(tf).exists(): Path(tf).unlink()

        # Free GPU memory after each multi-chain antibody
        torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== 3D STRUCTURE WHOLE MAB ASSEMBLY COMPLETE ===")