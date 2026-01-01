import os
# 1. FORCE GPU SECURITY BYPASS
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

try:
    torch.serialization.add_safe_globals([transformers.tokenization_utils.Trie, BertConfig, BertTokenizer])
except: pass

from igfold import IgFoldRunner
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- PATH CONFIG ---
# We point to the MAIN folder and will tell the script NOT to look inside subfolders
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")
OUTPUT_DIR = BASE_DIR / "PDB_Output_Files_GPU"

print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

def extract_chains_raw_text(file_path):
    """Reads .py files as text to find H and L chains."""
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
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
    
    # --- THE CRITICAL FIX ---
    # .iterdir() only looks at files in the MAIN folder. 
    # It will NOT enter "shadow_benchmarks" or any other subfolder.
    antibody_files = [
        f for f in BASE_DIR.iterdir() 
        if f.is_file() and f.suffix == ".py" and not f.name.startswith("._")
    ]
    
    print(f"\nSUCCESS: Found {len(antibody_files)} main antibody files in the root folder.")
    print(f"IGNORING: All subfolders (including shadow_benchmarks).\n")

    for f_path in tqdm(antibody_files, desc="Folding Progress"):
        ab_name = f_path.stem 
        seq_dict = extract_chains_raw_text(f_path)
        
        if not seq_dict:
            continue

        pdb_path = OUTPUT_DIR / f"{ab_name}.pdb"
        sasa_path = pdb_path.with_suffix(".sasa.txt")

        if sasa_path.exists():
            continue

        try:
            # GPU FOLDING
            runner.fold(
                pdb_file=str(pdb_path),
                sequences=seq_dict,
                do_refine=False,
                do_renum=False
            )

            # [2025-12-30] calculate_features (SASA)
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure(ab_name, str(pdb_path))
            sr = ShrakeRupley()
            sr.compute(struct, level="S")
            
            with open(sasa_path, "w") as f:
                f.write(f"{struct.sasa:.2f}")

        except Exception as e:
            print(f"\n[ERROR] Failed {ab_name}: {e}")

        # Keep 4080 SUPER VRAM clean
        torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n=== PROCESSING COMPLETE ===")