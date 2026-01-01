import os
# 1. FORCE SECURITY WALL DOWN GLOBALLY
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

import re
import glob
import torch
import warnings
import logging
from pathlib import Path
from io import StringIO
import importlib.util

# Suppress noise
logging.getLogger().setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# 2. ALLOWLIST SPECIFIC CLASSES FOR PYTORCH 2.6
import transformers
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer
try:
    torch.serialization.add_safe_globals([transformers.tokenization_utils.Trie, BertConfig, BertTokenizer])
except:
    pass # Handle if already added

from igfold import IgFoldRunner
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- CONFIG ---
INPUT_BASE = Path("shadow_benchmarks")
OUTPUT_BASE = Path("PDB_Output_Files_GPU")
FOLDERS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

def parse_fasta_to_pairs(fasta_string):
    """Uses the exact logic from your other file."""
    chains_by_number = {}
    handle = StringIO(fasta_string)
    for record in SeqIO.parse(handle, "fasta"):
        seq_id = record.id.lower()
        seq = str(record.seq)
        
        number_match = re.search(r'_(\d+)$', record.id)
        pair_num = number_match.group(1) if number_match else "1"
        
        if "heavy" in seq_id or "vhh" in seq_id:
            chains_by_number.setdefault(pair_num, {})["H"] = seq
        elif "light" in seq_id:
            chains_by_number.setdefault(pair_num, {})["L"] = seq
            
    return [c for c in chains_by_number.values() if "H" in c or "L" in c]

def run_pipeline():
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    
    for folder_name in FOLDERS:
        folder_path = INPUT_BASE / folder_name
        output_folder_path = OUTPUT_BASE / folder_name
        if not folder_path.exists(): continue
        os.makedirs(output_folder_path, exist_ok=True)

        print(f"\nProcessing Folder: {folder_name}")
        antibody_files = list(folder_path.glob("*.py"))

        for file_path in antibody_files:
            antibody_name = file_path.stem
            
            try:
                # DYNAMIC IMPORT: Exactly how your other file does it
                spec = importlib.util.spec_from_file_location(antibody_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                if hasattr(module, antibody_name):
                    sequence_string = getattr(module, antibody_name).strip()
                else:
                    continue

                paired_list = parse_fasta_to_pairs(sequence_string)

                for i, sequences_dict in enumerate(paired_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    output_path = output_folder_path / f"{pair_id}.pdb"
                    sasa_path = output_path.with_suffix(".sasa.txt")

                    if sasa_path.exists():
                        continue

                    print(f"Folding {pair_id} on GPU...")
                    runner.fold(
                        pdb_file=str(output_path),
                        sequences=sequences_dict,
                        do_refine=False,
                        do_renum=False
                    )

                    # [2025-12-30] CALCULATE FEATURES (SASA)
                    p = PDBParser(QUIET=True)
                    struct = p.get_structure(pair_id, str(output_path))
                    sr = ShrakeRupley()
                    sr.compute(struct, level="S")
                    
                    with open(sasa_path, "w") as f:
                        f.write(f"{struct.sasa:.2f}")

                torch.cuda.empty_cache()

            except Exception as e:
                print(f"Failed {antibody_name}: {e}")

if __name__ == "__main__":
    run_pipeline()
    print("\n=== FINISHED ALL FOLDERS ===")