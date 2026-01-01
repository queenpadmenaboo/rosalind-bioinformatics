import os
# 1. FORCE SECURITY BYPASS (Must be before torch is imported)
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

import re
import torch
import warnings
import logging
from pathlib import Path
from tqdm import tqdm

# Suppress noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

import transformers
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer

# 2. PATCH PYTORCH 2.6+ SECURITY
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
    """
    Reads .py files as RAW TEXT to bypass SyntaxErrors.
    Pulls sequences and groups them by Arm (_1, _2).
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regex: Find all FASTA blocks starting with >
    blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
    
    chains_by_number = {}
    for block in blocks:
        lines = block.strip().split('\n')
        if not lines: continue
        
        header = lines[0].strip()
        # Clean sequence: Keep ONLY A-Z letters
        sequence = re.sub(r'[^A-Z]', '', "".join(lines[1:]).upper())
        
        if not sequence: continue

        # Detect Arm 1 vs Arm 2 via _1 or _2 suffix in header
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

        # GATHER FILES AND PRINT COUNT FOR VERIFICATION
        antibody_files = list(folder_path.glob("*.py"))
        print(f"\nFOLDER: {folder_name} | FOUND: {len(antibody_files)} files")

        # Progress bar for the files in this specific folder
        for f_path in tqdm(antibody_files, desc=f"Folding {folder_name}"):
            # PULL NAME FROM FILENAME
            ab_name = f_path.stem 
            
            # Scrape pairs (Arm 1, Arm 2) from the text
            paired_list = extract_chains_raw_text(f_path)
            
            for i, seq_dict in enumerate(paired_list):
                # Format: filename_Pair_1
                pair_id = f"{ab_name}_Pair_{i+1}"
                pdb_path = output_folder / f"{pair_id}.pdb"
                sasa_path = pdb_path.with_suffix(".sasa.txt")

                if sasa_path.exists():
                    continue

                try:
                    # GPU FOLDING TASK
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

            # Keep the 4080 SUPER VRAM clear
            torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== PROCESSING COMPLETE. ALL PDBs AND SASA FILES GENERATED ===")