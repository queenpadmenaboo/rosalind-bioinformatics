import os
# 1. FORCE SECURITY BYPASS (Must be before torch is imported)
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

# 2. PATCH PYTORCH 2.6+ 
try:
    torch.serialization.add_safe_globals([transformers.tokenization_utils.Trie, BertConfig, BertTokenizer])
except:
    pass

from igfold import IgFoldRunner
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- CONFIG ---
INPUT_BASE = Path("shadow_benchmarks")
OUTPUT_BASE = Path("PDB_Output_Files_GPU")
FOLDERS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

def extract_chains_raw_text(file_path):
    """Reads .py files as RAW TEXT. Pulls sequences by Arm (_1, _2)."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regex: Find all FASTA blocks starting with >
    blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
    
    chains_by_number = {}
    for block in blocks:
        lines = block.strip().split('\n')
        if not lines: continue
        
        header = lines[0].strip()
        sequence = re.sub(r'[^A-Z]', '', "".join(lines[1:]).upper())
        
        if not sequence: continue

        num_match = re.search(r'_(\d+)', header)
        pair_num = num_match.group(1) if num_match else "1"
        
        if pair_num not in chains_by_number:
            chains_by_number[pair_num] = {}
        
        h_labels = ["heavy", "vhh", "vh"]
        l_labels = ["light", "vl"]
        
        header_lower = header.lower()
        if any(label in header_lower for label in h_labels):
            chains_by_number[pair_num]["H"] = sequence
        elif any(label in header_lower for label in l_labels):
            chains_by_number[pair_num]["L"] = sequence
            
    return [c for c in chains_by_number.values() if "H" in c and "L" in c]

def run_pipeline():
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    
    for folder_name in FOLDERS:
        folder_path = INPUT_BASE / folder_name
        output_folder = OUTPUT_BASE / folder_name
        
        if not folder_path.exists():
            continue
            
        os.makedirs(output_folder, exist_ok=True)

        # --- FIX: Only look at files in the TOP LEVEL of the folder ---
        # This stops the script from finding 937 files if only 471 exist.
        antibody_files = [f for f in folder_path.iterdir() if f.is_file() and f.suffix == ".py"]
        
        print(f"\nFOLDER: {folder_name} | TRUE COUNT: {len(antibody_files)} files")

        for f_path in tqdm(antibody_files, desc=f"Folding {folder_name}"):
            ab_name = f_path.stem 
            paired_list = extract_chains_raw_text(f_path)
            
            for i, seq_dict in enumerate(paired_list):
                pair_id = f"{ab_name}_Pair_{i+1}"
                pdb_path = output_folder / f"{pair_id}.pdb"
                sasa_path = pdb_path.with_suffix(".sasa.txt")

                if sasa_path.exists():
                    continue

                try:
                    runner.fold(
                        pdb_file=str(pdb_path),
                        sequences=seq_dict,
                        do_refine=False,
                        do_renum=False
                    )

                    # [2025-12-30] CALCULATE FEATURES (SASA)
                    parser = PDBParser(QUIET=True)
                    struct = parser.get_structure(pair_id, str(pdb_path))
                    sr = ShrakeRupley()
                    sr.compute(struct, level="S")
                    
                    with open(sasa_path, "w") as f:
                        f.write(f"{struct.sasa:.2f}")

                except Exception as e:
                    print(f"\n[ERROR] Failed {pair_id}: {e}")

            torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== PROCESSING COMPLETE ===")