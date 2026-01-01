import os
import re
import torch
import torch.serialization
import collections
from pathlib import Path
from tqdm import tqdm
from io import StringIO
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# =================================================================
# 1. THE "SMART" SECURITY OVERRIDE (For PyTorch 2.6+)
# This fixes the "Multiple values for keyword argument" error.
# =================================================================
original_load = torch.serialization.load

def patched_load(*args, **kwargs):
    # This prevents the conflict by safely overwriting the flag
    kwargs.pop('weights_only', None) 
    return original_load(*args, weights_only=False, **kwargs)

# Apply the patch to the core engine
torch.load = patched_load
torch.serialization.load = patched_load

# =================================================================
# 2. INITIALIZE ENGINE
# =================================================================
from igfold import IgFoldRunner
print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

# --- SETTINGS ---
INPUT_BASE = Path("shadow_benchmarks")
OUTPUT_BASE = Path("PDB_Output_Files_GPU")
FOLDERS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

def extract_fasta_via_regex(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        match = re.search(r'["\']{3}([\s\S]*?>[\s\S]*?)["\']{3}', content)
        return match.group(1).strip() if match else None
    except Exception:
        return None

def parse_chains(fasta_str):
    records = list(SeqIO.parse(StringIO(fasta_str), "fasta"))
    chains = {} 
    for r in records:
        num_match = re.search(r'_(\d+)$', r.id)
        num = num_match.group(1) if num_match else "1"
        if num not in chains: chains[num] = {}
        if "heavy" in rid or "vhh" in rid: chains[num]["H"] = str(r.seq)
        elif "light" in rid: chains[num]["L"] = str(r.seq)
    return [c for c in chains.values() if "H" in c or "L" in c]

def run_pipeline():
    os.makedirs(OUTPUT_BASE, exist_ok=True)
    for folder_name in FOLDERS:
        folder_path = INPUT_BASE / folder_name
        output_folder = OUTPUT_BASE / folder_name
        if not folder_path.exists(): continue
        os.makedirs(output_folder, exist_ok=True)
        
        files = list(folder_path.glob("*.py"))
        for f_path in tqdm(files, desc=f"Processing {folder_name}"):
            ab_name = f_path.stem
            fasta_str = extract_fasta_via_regex(f_path)
            if not fasta_str: continue
            
            tasks = parse_chains(fasta_str)
            for i, seq_dict in enumerate(tasks):
                pdb_fn = output_folder / f"{ab_name}_Pair_{i+1}.pdb"
                sasa_fn = pdb_fn.with_suffix(".sasa.txt")
                
                if sasa_fn.exists(): continue 
                
                try:
                    runner.fold(str(pdb_fn), seq_dict, do_refine=False, do_renum=False)
                    
                    # [2025-12-30 Instruction: calculate_features]
                    parser = PDBParser(QUIET=True)
                    struct = parser.get_structure(ab_name, str(pdb_fn))
                    sr = ShrakeRupley()
                    sr.compute(struct, level="S")
                    
                    with open(sasa_fn, "w") as sf:
                        sf.write(f"{struct.sasa:.2f}")
                except Exception as e:
                    print(f"\n[ERROR] {ab_name}: {e}")
            torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()