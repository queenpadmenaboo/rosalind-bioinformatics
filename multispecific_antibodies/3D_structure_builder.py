import os
import glob
from io import StringIO
from pathlib import Path
import importlib.util
import re
import warnings
import logging
logging.getLogger().setLevel(logging.ERROR)

# Suppress warnings
warnings.filterwarnings("ignore")

import torch
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer
from igfold import IgFoldRunner
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import transformers

# Allowlist the specific transformer classes that PyTorch 2.6 now blocks by default
torch.serialization.add_safe_globals([transformers.tokenization_utils.Trie])


# Allowlist BertConfig for safe checkpoint loading (IgFold / AntiBERTy)
torch.serialization.add_safe_globals([BertConfig, BertTokenizer])

# --- Configuration Variables ---
USE_REFINEMENT = False
OUTPUT_BASE_DIR = Path("PDB_Output_Files")
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

print("Initializing IgFold...", flush=True)
runner = IgFoldRunner()
print("IgFold ready.\n", flush=True)

def parse_fasta_string_to_pairs(fasta_string):
    chains_by_number = {}       # if python sees _1: it creates a "Bin 1"; puts the Heavy chain in Bin 1["H"] and the Light chain in Bin 1["L"].
    handle = StringIO(fasta_string)
    for record in SeqIO.parse(handle, "fasta"):
        seq_id = record.id
        seq = str(record.seq)

        # Extract pair number from end (e.g., "_1" or "_2")
        number_match = re.search(r'_(\d+)$', seq_id)
        pair_num = number_match.group(1) if number_match else "1"

        is_heavy = "heavy" in seq_id.lower()
        is_light = "light" in seq_id.lower()

        if is_heavy:
            chains_by_number.setdefault(pair_num, {})["H"] = seq
        elif is_light:
            chains_by_number.setdefault(pair_num, {})["L"] = seq

    # Return only valid pairs
    return [chains for chains in chains_by_number.values() if "H" in chains and "L" in chains]  # Bin 1: Does it have an H and L? Bin 2: Does it have an H and an L? # of "Valid Pairs" is determined

def process_directory(base_dir, subfolders):
    base_path = Path(base_dir)

    # Dictionary of notes for each folder type
    FOLDER_NOTES = {
        "Bispecific_mAb": "NOTE: Expecting 2 different H/L pairs (Pair 1 and Pair 2) for asymmetric arms.",
        "Bispecific_scFv": "NOTE: Processing single-chain variable fragments; pairing depends on FASTA IDs.",
        "Other_Formats": "NOTE: Processing non-standard antibody formats.",
        "Whole_mAb": "NOTE: Standard IgG. Expecting 1 unique H/L pair. (Identity math: 1 model = 2 identical arms)."
    }

    for folder_name in subfolders:
        output_folder_path = OUTPUT_BASE_DIR / folder_name
        os.makedirs(output_folder_path, exist_ok=True)

    total_processed = 0
    total_skipped = 0

    for folder_name in subfolders:
        folder_path = base_path / folder_name
        if not folder_path.exists():
            print(f"Input folder not found: {folder_path}, skipping.", flush=True)
            continue

        # --- Print Folder Note ---
        note = FOLDER_NOTES.get(folder_name, "Processing files...")
        print(f"\n{'='*60}")
        print(f"FOLDER: {folder_name}")
        print(f"{note}")
        print(f"{'='*60}\n", flush=True)

        antibody_files = glob.glob(str(folder_path / "*.py"))

        for file_path in antibody_files:
            file_path = Path(file_path)
            antibody_name = file_path.stem

            try:
                spec = importlib.util.spec_from_file_location(antibody_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                if hasattr(module, antibody_name):
                    sequence_string = getattr(module, antibody_name)
                    sequence_string = sequence_string.strip() 
                else:
                    print(f"Failed to find variable '{antibody_name}' in {file_path}. Skipping.", flush=True)
                    continue

                paired_sequences_list = parse_fasta_string_to_pairs(sequence_string)
                output_folder_path = OUTPUT_BASE_DIR / folder_name

                for i, sequences_dict in enumerate(paired_sequences_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    output_path = output_folder_path / f"{pair_id}.pdb"
                    
                    if output_path.exists():
                        print(f"Skipping {pair_id} (already exists)", flush=True)
                        total_skipped += 1
                        continue
                    
                    print(f"Predicting structure for {pair_id}...", flush=True)

                    runner.fold(
                        pdb_file=str(output_path),
                        sequences=sequences_dict,
                        do_refine=USE_REFINEMENT,
                        do_renum=False,     # Keep this False for Windows compatibility
                    )
                    print(f"Saved {output_path}", flush=True)
                    

                    try:
                        p = PDBParser(QUIET=True)
                        struct = p.get_structure(pair_id, str(output_path))
                        
                        sr = ShrakeRupley()
                        sr.compute(struct, level="S") # "S" for total structure SASA
    
                        total_sasa = struct.sasa
                        print(f"SASA for {pair_id}: {total_sasa:.2f} Å²", flush=True)
    
                        # Optional: Save SASA to a text file next to the PDB
                        with open(output_path.with_suffix(".sasa.txt"), "w") as f:
                            f.write(f"Total SASA: {total_sasa:.2f} Å²\n")
        
                    except Exception as sasa_error:
                        print(f"SASA calculation failed for {pair_id}: {sasa_error}", flush=True)
    
                    total_processed += 1

            except Exception as e:
                print(f"Failed to process {file_path}: {e}", flush=True)

    print(f"\n=== DONE: {total_processed} structures created, {total_skipped} skipped ===", flush=True)


if __name__ == "__main__":
    process_directory(base_dir=".", subfolders=FOLDERS_TO_PROCESS)