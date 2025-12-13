import os
import glob
from io import StringIO
from pathlib import Path
import importlib.util

import torch
from transformers.models.bert.configuration_bert import BertConfig
from transformers.models.bert.tokenization_bert import BertTokenizer
from igfold import IgFoldRunner
from Bio import SeqIO
import os

# Allowlist BertConfig for safe checkpoint loading (IgFold / AntiBERTy)
torch.serialization.add_safe_globals([BertConfig, BertTokenizer])

# --- Configuration Variables ---
USE_REFINEMENT = True
OUTPUT_BASE_DIR = Path("PDB_Output_Files")
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]
# Add any other folders you need to process here

runner = IgFoldRunner()

def parse_fasta_string_to_pairs(fasta_string):
    chains_by_number = {}
    handle = StringIO(fasta_string)
    # MINIMAL CHANGE 1: Use fasta-blast format to ignore comments/leading newlines
    for record in SeqIO.parse(handle, "fasta-blast"): 
        seq_id = record.id
        seq = str(record.seq)

        # Robust checks for both full names and _h/_l suffixes
        is_heavy = "_Heavy_Chain" in seq_id or "_h" in seq_id.lower() or "_h1" in seq_id.lower() or "heavy" in seq_id.lower()
        is_light = "_Light_Chain" in seq_id or "_l" in seq_id.lower() or "_l1" in seq_id.lower() or "light" in seq_id.lower()

        if is_heavy or is_light:
            # Simple approach: assume single H/L pair per file if names are unique enough
            chain_type = "H" if is_heavy else "L"
            # Use a generic key if complex numbering isn't required
            chains_by_number.setdefault("1", {})[chain_type] = seq

    # Return only complete pairs
    return [chains for chains in chains_by_number.values() if "H" in chains and "L" in chains]


def process_directory(base_dir, subfolders):
    base_path = Path(base_dir)

    # Ensure the output folders exist
    for folder_name in subfolders:
        output_folder_path = OUTPUT_BASE_DIR / folder_name
        os.makedirs(output_folder_path, exist_ok=True)

    # Process input folders if they exist
    for folder_name in subfolders:
        folder_path = base_path / folder_name
        if not folder_path.exists():
            print(f"Input folder not found: {folder_path}, skipping.")
            continue

        antibody_files = glob.glob(str(folder_path / "*.py"))

        for file_path in antibody_files:
            file_path = Path(file_path)
            antibody_name = file_path.stem

            # Skip files that are in your general EXCLUDE_FILES list if needed
            # (You would need to import EXCLUDE_FILES here if you want that functionality)

            try:
                # Safely load the Python file as a module using importlib
                spec = importlib.util.spec_from_file_location(antibody_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)

                # Get the sequence string variable directly from the module object
                if hasattr(module, antibody_name):
                    sequence_string = getattr(module, antibody_name)
                else:
                    print(f"Failed to find variable '{antibody_name}' in {file_path}. Skipping.")
                    continue


                paired_sequences_list = parse_fasta_string_to_pairs(sequence_string)
                output_folder_path = OUTPUT_BASE_DIR / folder_name

                for i, sequences_dict in enumerate(paired_sequences_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    print(f"Predicting structure for {pair_id}...")

                    # MINIMAL CHANGE 2: Call .fold (not .predict) AND include pdb_file=None argument
                    predicted_structure = runner.fold(
                        pdb_file=None,
                        sequences=sequences_dict,
                        do_refine=USE_REFINEMENT,
                    )
                    
                    # MINIMAL CHANGE 3: Add safety check to prevent NoneType crash
                    if predicted_structure is not None:
                        output_path = output_folder_path / f"{pair_id}.pdb"
                        predicted_structure.save(str(output_path))
                        print(f"Saved {output_path}")
                    else:
                        print(f"Prediction failed for {pair_id}, no structure returned by runner.fold().")

            except Exception as e:
                print(f"Failed to process {file_path}: {e}")


if __name__ == "__main__":
    # Ensure your current working directory is correct when calling the script
    process_directory(base_dir=".", subfolders=FOLDERS_TO_PROCESS)