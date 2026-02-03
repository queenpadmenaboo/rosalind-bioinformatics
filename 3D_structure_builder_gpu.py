import os
# 1. GPU SECURITY BYPASS (Required for 4080 SUPER / Torch 2.6+)
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"
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

# --- CONFIG & EXCLUSIONS ---
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")
OUTPUT_ROOT = BASE_DIR / "PDB_Output_Files_GPU"

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
    'Antibody_Chain_Count_Summary.xlsx'
}

print("--- INITIALIZING 4080 SUPER ENGINE (STABLE 3.10 MODE) ---")
runner = IgFoldRunner()

def extract_chains_dynamic(file_path):
    """
    1. Reads the file.
    2. Identifies headers (>).
    3. Filters out "NA" and empty sequences.
    4. Dynamically assigns H1, L1, H2, L2 etc.
    """
    found_seq_dict = {}
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
        
        h_idx, l_idx = 1, 1
        for block in blocks:
            lines = block.strip().split('\n')
            if len(lines) < 2: continue
            
            header = lines[0].lower()
            raw_seq = "".join(lines[1:]).strip().upper()
            clean_seq = re.sub(r'[^A-Z]', '', raw_seq).upper()
            
            # The "NA" protector
            if not clean_seq or clean_seq == "NA":
                continue
            
            # Assign keys based on header keywords
            if any(h in header for h in ["heavy", "vh", "vhh"]):
                key = f"H{h_idx}"
                h_idx += 1
            else:
                key = f"L{l_idx}"
                l_idx += 1
            
            found_seq_dict[key] = clean_seq
    except:
        pass
    return found_seq_dict if found_seq_dict else None    

def run_pipeline():
    os.makedirs(OUTPUT_ROOT, exist_ok=True)
        
    antibody_files = [
        f for f in BASE_DIR.rglob("*.py") 
        if f.name not in EXCLUDE_FILES 
        and "shadow_benchmarks" not in f.parts 
        and "PDB_Output_Files_GPU" not in f.parts
        and not f.name.startswith("._")
    ]
    
    print(f"\nSUCCESS: Found {len(antibody_files)} valid antibody files.")

    for f_path in tqdm(antibody_files, desc="Folding"):
        relative_path = f_path.parent.relative_to(BASE_DIR)
        target_dir = OUTPUT_ROOT / relative_path
        target_dir.mkdir(parents=True, exist_ok=True)

        ab_name = f_path.stem 
        pdb_path = target_dir / f"{ab_name}.pdb"

        # --- THE SKIP CHECK (Safe to comment out if folder is deleted) ---
        # if pdb_path.exists():
        #      continue 

        seq_dict = extract_chains_dynamic(f_path)
        
        if not seq_dict:
            continue

        try:
            # The runner folds all chains in the dict into one unified PDB
            runner.fold(
                pdb_file=str(pdb_path),
                sequences=seq_dict,
                do_refine=False, 
                do_renum=False
            )
            
            # CRITICAL CHECK: Ensure the file actually has data
            if pdb_path.exists() and pdb_path.stat().st_size < 100:
                print(f"\n[WARNING] {ab_name} produced an empty PDB. GPU Attention error?")
                pdb_path.unlink() # Delete the bad empty file
        
        except Exception as e:
            
            if "Invalid or missing coordinate" in str(e):
                 print(f"\n[GPU ERROR] {ab_name} failed coordinate generation (Check VRAM/Weights).")
            else:
                 print(f"\n[ERROR] Skipping {ab_name}: {e}")

        torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== PROCESSING COMPLETE ===")