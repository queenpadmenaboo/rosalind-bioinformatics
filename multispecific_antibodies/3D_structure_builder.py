import os
import glob
from io import StringIO
from pathlib import Path
import importlib.util
import re

import torch
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer
from igfold import IgFoldRunner
from Bio import SeqIO

# Allowlist BertConfig for safe checkpoint loading (IgFold / AntiBERTy)
torch.serialization.add_safe_globals([BertConfig, BertTokenizer])

# --- Configuration Variables ---
USE_REFINEMENT = False
OUTPUT_BASE_DIR = Path("PDB_Output_Files")
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

runner = IgFoldRunner()

def parse_fasta_string_to_pairs(fasta_string):
    chains_by_number = {}
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

    # Return only complete pairs
    return [chains for chains in chains_by_number.values() if "H" in chains and "L" in chains]

def process_directory(base_dir, subfolders):
    base_path = Path(base_dir)

    for folder_name in subfolders:
        output_folder_path = OUTPUT_BASE_DIR / folder_name
        os.makedirs(output_folder_path, exist_ok=True)

    for folder_name in subfolders:
        folder_path = base_path / folder_name
        if not folder_path.exists():
            print(f"Input folder not found: {folder_path}, skipping.")
            continue

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
                    print(f"Failed to find variable '{antibody_name}' in {file_path}. Skipping.")
                    continue

                paired_sequences_list = parse_fasta_string_to_pairs(sequence_string)
                output_folder_path = OUTPUT_BASE_DIR / folder_name

                for i, sequences_dict in enumerate(paired_sequences_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    print(f"Predicting structure for {pair_id}...")

                    output_path = output_folder_path / f"{pair_id}.pdb"
                    predicted_structure = runner.fold(
                        pdb_file=str(output_path),
                        sequences=sequences_dict,
                        do_refine=USE_REFINEMENT,
                    )
                    print(f"Saved {output_path}")
                    

            except Exception as e:
                print(f"Failed to process {file_path}: {e}")


if __name__ == "__main__":
    process_directory(base_dir=".", subfolders=FOLDERS_TO_PROCESS)
