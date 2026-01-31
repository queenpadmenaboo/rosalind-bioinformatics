import os
import pandas as pd
import re
import warnings
import logging
from pathlib import Path
from tqdm import tqdm
import subprocess
import json

# Suppress library noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# --- CONFIG & EXCLUSIONS ---
BASE_DIR = Path("/mnt/c/Users/bunsr/rosalind-bioinformatics/multispecific_antibodies/Whole_mAb")
OUTPUT_ROOT = BASE_DIR / "PDB_Output_ColabFold_Structures"
CSV_PATH = Path("/mnt/c/Users/bunsr/TheraSAbDab_SeqStruc_07Dec2025.csv")
TEMP_FASTA_DIR = BASE_DIR / "temp_fastas"

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
    'whole_mAbs_folder_check.py', 'whole_mAbs_isotypes_check.py', 'whole_mAb_antibody_list.txt',
    '3D_structure_builder_colabfold_whole_mAbs.py', 'colabfold_gpu_whole_mAbs.py'
}

# IUPAC Standard 20 Amino Acids
AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# ============================================================================
# TEST MODE CONFIGURATION
# ============================================================================
# Set TEST_MODE = True to run only first few antibodies for validation
# Set TEST_MODE = False to run full batch (all antibodies)
TEST_MODE = True
MAX_TEST_ANTIBODIES = 2  # Only used when TEST_MODE = True
# ============================================================================

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
    """
    Extracts H and L sequences from Python files containing FASTA format.
    Handles format: variable_name = '''\\n>Header\\nSequence\\n'''
    """
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
            if clean_seq and clean_seq != "NA" and len(clean_seq) > 20:
                # Validate it's actually amino acids (not random Python code)
                invalid_chars = set(clean_seq) - set(AA_ALPHABET)
                if invalid_chars:
                    print(f"Warning: Invalid amino acids {invalid_chars} in {file_path.stem}, skipping this sequence")
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
    
    # Separate and pair chains
    heavies = [s for t, s in found_sequences if t == 'H']
    lights = [s for t, s in found_sequences if t == 'L']
    
    # Validate we found chains
    if not heavies or not lights:
        print(f"Warning: {file_path.stem} missing chains (H={len(heavies)}, L={len(lights)})")
        return None
    
    pairs = list(zip(heavies, lights))
    
    # If Whole_mAb has only 1 pair, duplicate it for 2H:2L biological assembly
    if len(pairs) == 1:
        pairs = pairs * 2  # Creates [(VH, VL), (VH, VL)]
    
    return pairs if pairs else None

def create_multimer_fasta(pairs, h_iso, l_iso, output_fasta):
    """
    Creates a single-entry FASTA file for ColabFold multimer prediction.
    Sequences are separated by a colon (:) for complex assembly.
    """
    all_chains = []
    for i, (heavy_fv, light_fv) in enumerate(pairs):
        # Construct full-length chains
        full_heavy = heavy_fv + CONSTANTS['HEAVY'][h_iso]
        full_light = light_fv + CONSTANTS['LIGHT'][l_iso]
        
        all_chains.append(full_heavy)
        all_chains.append(full_light)
    
    # CRITICAL: Join sequences with ':' so ColabFold treats them as one multimer
    complex_sequence = ":".join(all_chains)
    
    with open(output_fasta, 'w') as f:
        # Use one header for the entire complex
        f.write(f">Whole_mAb_Assembly\n{complex_sequence}\n")

def run_colabfold(fasta_path, output_dir, antibody_name):
    """
    Runs ColabFold batch prediction using colabfold_batch command.
    Optimized for 4080 Super: Fast batch processing (3-7 min per antibody).
    """
    cmd = [
        "colabfold_batch",
        str(fasta_path),
        str(output_dir),
        "--num-models", "1",  # Use best model only (not all 5)
        "--num-recycle", "3", 
        "--model-type", "alphafold2_multimer_v3",  # Use multimer for assembly
        "--max-msa", "32:64",
        "--amber",
        "--use-gpu-relax",
        "--rank", "plddt",  # Rank by pLDDT confidence
    ]
    
    try:
        print(f"Running ColabFold for {antibody_name}...")
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"[SUCCESS] ColabFold completed for {antibody_name}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[FAILED] ColabFold failed for {antibody_name}")
        print(f"Error: {e.stderr}")
        return False

def rename_output_pdb(output_dir, antibody_name, final_output_path):
    """
    ColabFold outputs: antibody_relaxed_rank_001_*.pdb (with AMBER relaxation)
    Rename to: antibody_name.pdb
    """
    pdb_files = list(output_dir.glob("*_relaxed_rank_001_*.pdb"))
    
    if pdb_files:
        source_pdb = pdb_files[0]
        source_pdb.rename(final_output_path)
        print(f"[SAVED] {final_output_path}")
        
        # Clean up intermediate files
        for pattern in ["*.a3m", "*.json", "*_rank_00[2-5]_*.pdb"]:
            for f in output_dir.glob(pattern):
                f.unlink()
    else:
        print(f"[ERROR] No output PDB found for {antibody_name}")

def run_pipeline():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
    os.makedirs(TEMP_FASTA_DIR, exist_ok=True)

    # Load Master Isotype Data from CSV
    isotype_lookup = {}
    try:
        df = pd.read_csv(CSV_PATH)
        isotype_lookup = dict(zip(df['Therapeutic'].str.lower(), zip(df['CH1 Isotype'], df['VD LC'])))
        print(f"--- LOADED {len(isotype_lookup)} ISOTYPES FROM CSV ---")
    except Exception as e:
        print(f"--- ERROR LOADING CSV: {e}. Falling back to keyword search only. ---")
       
    antibody_files = [
        f for f in BASE_DIR.rglob("*.py") 
        if f.name not in EXCLUDE_FILES 
        and "PDB_Output_Files" not in str(f)
        and not f.name.startswith("._")
    ]
    
    # Apply test mode limit if enabled
    if TEST_MODE:
        antibody_files = antibody_files[:MAX_TEST_ANTIBODIES]
        print(f"\n[WARNING] TEST MODE ENABLED: Processing only {len(antibody_files)} antibodies")
        print(f"    To run full batch, set TEST_MODE = False in script\n")
    
    print(f"\nSUCCESS: STARTING COLABFOLD ASSEMBLY: {len(antibody_files)} antibody files in Whole_mAb.\n")

    for f_path in tqdm(antibody_files, desc="Building 3D Structures with ColabFold"):
        tqdm.write(f"\n{'='*60}")
        tqdm.write(f"FOLDING: {f_path.stem}")
        tqdm.write(f"{'='*60}")
        
        # Get sequences
        pairs = extract_chains_dynamic(f_path)
        if not pairs:
            tqdm.write(f"[ERROR] No valid sequences found in {f_path.stem}")
            continue
        
        # Detect isotype
        h_iso, l_iso = detect_isotype(f_path, isotype_lookup)
        tqdm.write(f"Isotype: Heavy={h_iso}, Light={l_iso}")
        tqdm.write(f"Chains: {len(pairs)} H-L pairs")
        
        # Create multimer FASTA
        fasta_path = TEMP_FASTA_DIR / f"{f_path.stem}.fasta"
        create_multimer_fasta(pairs, h_iso, l_iso, fasta_path)
        
        # Create antibody-specific output directory
        antibody_output = OUTPUT_ROOT / f_path.stem
        os.makedirs(antibody_output, exist_ok=True)
        
        # Run ColabFold
        success = run_colabfold(fasta_path, antibody_output, f_path.stem)
        
        if success:
            # Move and rename final PDB
            final_pdb = OUTPUT_ROOT / f"{f_path.stem}.pdb"
            rename_output_pdb(antibody_output, f_path.stem, final_pdb)
        
        # Cleanup temp FASTA
        if fasta_path.exists():
            fasta_path.unlink()
        
        tqdm.write("")  # Blank line between antibodies

    print("\n" + "="*60)
    print("=== 3D STRUCTURE WHOLE MAB ASSEMBLY COMPLETE ===")
    print("="*60)

if __name__ == "__main__":
    run_pipeline()