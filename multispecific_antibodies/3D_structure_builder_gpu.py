import os
import re
import torch
import torch.serialization
import collections
from pathlib import Path
from tqdm import tqdm
from io import StringIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# 1. THE "NUCLEAR" SECURITY FIX (For PyTorch 2.6+)
# Overrides weights_only to False to stop the UnpicklingError
original_load = torch.serialization.load
def patched_load(*args, **kwargs):
    kwargs.pop('weights_only', None) 
    return original_load(*args, weights_only=False, **kwargs)

torch.load = patched_load
torch.serialization.load = patched_load

# 2. START GPU ENGINE
from igfold import IgFoldRunner
print("--- INITIALIZING 4080 SUPER ENGINE ---")
runner = IgFoldRunner()

# --- CONFIG ---
INPUT_BASE = Path("shadow_benchmarks")
OUTPUT_BASE = Path("PDB_Output_Files_GPU")
FOLDERS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

def extract_tasks_directly(file_path):
    """
    ULTRA-AGGRESSIVE EXTRACTION:
    Finds every header starting with > and grabs every letter following it.
    Ignores all Python quotes, commas, variable names, and case.
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Grab blocks: everything from > to the next > (or end of file)
    blocks = re.findall(r'>(.*?)(?=>|\Z)', content, re.DOTALL)
    
    chains = {}
    for block in blocks:
        lines = block.strip().split('\n')
        header = lines[0].strip().lower()
        # Keep ONLY letters A-Z (removes quotes, commas, spaces, numbers)
        sequence = re.sub(r'[^A-Z]', '', "".join(lines[1:]).upper())
        
        if not sequence: continue

        # Detect Arm 1 vs Arm 2 (_1, _2)
        # Suffix $ removed to handle cases where quotes/commas follow the header
        num_match = re.search(r'_(\d+)', header)
        n = num_match.group(1) if num_match else "1"
        
        if n not in chains: 
            chains[n] = {}
        
        # Case-insensitive chain identification
        h_labels = ["heavy", "vhh", "vh"]
        l_labels = ["light", "vl"]
        
        if any(label in header for label in h_labels):
            chains[n]["H"] = sequence
        elif any(label in header for label in l_labels):
            chains[n]["L"] = sequence
            
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
            
            # Extract sequences from the .py file content
            tasks = extract_tasks_directly(f_path)
            if not tasks: 
                continue # If list is empty, skip to next file
            
            for i, seq_dict in enumerate(tasks):
                pdb_fn = output_folder / f"{ab_name}_Pair_{i+1}.pdb"
                sasa_fn = pdb_fn.with_suffix(".sasa.txt")
                
                # Resumable: skip if SASA result already exists
                if sasa_fn.exists(): continue 
                
                try:
                    # STEP A: GPU 3D FOLDING (Using 4080 SUPER)
                    runner.fold(str(pdb_fn), seq_dict, do_refine=False, do_renum=False)
                    
                    # STEP B: CALCULATE FEATURES (SASA) - Per Dec 30 Instruction
                    parser = PDBParser(QUIET=True)
                    struct = parser.get_structure(ab_name, str(pdb_fn))
                    sr = ShrakeRupley()
                    sr.compute(struct, level="S")
                    
                    # Save numeric SASA result
                    with open(sasa_fn, "w") as sf:
                        sf.write(f"{struct.sasa:.2f}")
                        
                except Exception as e:
                    print(f"\n[ERROR] {ab_name} (Pair {i+1}): {e}")

            # Keep 4080 Super VRAM clean
            torch.cuda.empty_cache()

if __name__ == "__main__":
    run_pipeline()
    print("\n--- ALL FOLDING AND SASA CALCULATIONS COMPLETE ---")