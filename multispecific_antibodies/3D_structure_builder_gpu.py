import os
# 1. GPU SECURITY BYPASS (Required for 4080 SUPER / Torch 2.6+)
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

import re
import torch
import warnings
import logging
from pathlib import Path
from tqdm import tqdm

# Suppress library noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

import transformers
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer

# PATCH FOR PYTORCH 2.6+
try:
    torch.serialization.add_safe_globals([transformers.tokenization_utils.Trie, BertConfig, BertTokenizer])
except:
    pass

from igfold import IgFoldRunner
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- CONFIG & EXCLUSIONS ---
# Targeted specifically to your bioinformatics directory
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")
OUTPUT_DIR = BASE_DIR / "PDB_Output_Files_GPU"

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
    'PNAS_VS_CALCULATIONS.xlsx', '3D_structure_builder_gpu.py'
}

print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

def extract_chains_raw_text(file_path):
    """Parses .py files for FASTA sequences, ignoring code syntax."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Matches headers like >Heavy and captures sequence until next quote or header
    blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
    
    chains = {"H": None, "L": None}
    for block in blocks:
        lines = block.strip().split('\n')
        if not lines: continue
        
        header = lines[0].lower()
        sequence = re.sub(r'[^A-Z]', '', "".join(lines[1:]).upper())
        
        if not sequence: continue

        if any(h in header for h in ["heavy", "vh", "vhh"]):
            chains["H"] = sequence
        elif any(l in header for l in ["light", "vl"]):
            chains["L"] = sequence
            
    return chains if (chains["H"] and chains["L"]) else None

def run_pipeline():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # SEARCH LOGIC: 
    # 1. rglob("*.py") finds files in ALL subfolders
    # 2. filename check blocks your utility scripts
    # 3. 'parts' check blocks the shadow_benchmarks folder
    antibody_files = [
        f for f in BASE_DIR.rglob("*.py") 
        if f.name not in EXCLUDE_FILES 
        and "shadow_benchmarks" not in f.parts 
        and not f.name.startswith("._")
    ]
    
    print(f"\nSUCCESS: Found {len(antibody_files)} valid antibody files.")
    print(f"SHIELD ACTIVE: {len(EXCLUDE_FILES)} utility files ignored.\n")

    for f_path in tqdm(antibody_files, desc="Processing"):
        ab_name = f_path.stem 
        seq_dict = extract_chains_raw_text(f_path)
        
        # Skip if the file is just code (no sequences found)
        if not seq_dict:
            continue

        pdb_path = OUTPUT_DIR / f"{ab_name}.pdb"
        sasa_path = pdb_path.with_suffix(".sasa.txt")

        # Skip already completed work
        if sasa_path.exists():
            continue

        try:
            # 1. FOLD
            runner.fold(
                pdb_file=str(pdb_path),
                sequences=seq_dict,
                do_refine=False,
                do_renum=False
            )

            # 2. CALCULATE FEATURES (SASA) [2025-12-30]
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(ab_name, str(pdb_path))
            sr = ShrakeRupley()
            sr.compute(struct, level="S")
            
            with open(sasa_path, "w") as f:
                f.write(f"{struct.sasa:.2f}")

        except Exception as e:
            print(f"\n[ERROR] Skipping {ab_name}: {e}")

        # Clear VRAM for the 4080 SUPER
        torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== PROCESSING COMPLETE ===")