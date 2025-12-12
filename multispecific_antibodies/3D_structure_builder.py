import os
import glob
from io import StringIO

import torch
from transformers.models.bert.configuration_bert import BertConfig
from igfold import IgFoldRunner
from Bio import SeqIO

# Allowlist BertConfig for safe checkpoint loading (IgFold / AntiBERTy)
torch.serialization.add_safe_globals([BertConfig])

# --- Configuration Variables ---
USE_REFINEMENT = True
OUTPUT_BASE_DIR = "PDB_Output_Files"
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

runner = IgFoldRunner()


def parse_fasta_string_to_pairs(fasta_string):
    chains_by_number = {}
    handle = StringIO(fasta_string)
    for record in SeqIO.parse(handle, "fasta"):
        seq_id = record.id
        seq = str(record.seq)

        if "_Heavy_Chain" in seq_id:
            chain_num = seq_id.split("_Heavy_Chain_")[-1]
            chains_by_number.setdefault(chain_num, {})["H"] = seq
        elif "_Light_Chain" in seq_id:
            chain_num = seq_id.split("_Light_Chain_")[-1]
            chains_by_number.setdefault(chain_num, {})["L"] = seq

    return [chains for chains in chains_by_number.values() if "H" in chains and "L" in chains]


def process_directory(base_dir, subfolders):
    # Ensure the output folders exist
    for folder_name in subfolders:
        output_folder_path = os.path.join(OUTPUT_BASE_DIR, folder_name)
        os.makedirs(output_folder_path, exist_ok=True)

    # Process input folders if they exist
    for folder_name in subfolders:
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.exists(folder_path):
            print(f"Input folder not found: {folder_path}, skipping.")
            continue

        antibody_files = glob.glob(os.path.join(folder_path, "*.py"))

        for file_path in antibody_files:
            antibody_name, _ = os.path.splitext(os.path.basename(file_path))
            try:
                with open(file_path, "r") as f:
                    file_content = f.read()
                # Execute the file to get the sequence variable
                exec(file_content, globals())
                sequence_string = globals()[antibody_name]

                paired_sequences_list = parse_fasta_string_to_pairs(sequence_string)
                output_folder_path = os.path.join(OUTPUT_BASE_DIR, folder_name)

                for i, sequences_dict in enumerate(paired_sequences_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    predicted_structure = runner.predict(
                        sequences=sequences_dict,
                        do_refine=USE_REFINEMENT,
                    )
                    output_path = os.path.join(output_folder_path, f"{pair_id}.pdb")
                    predicted_structure.save(output_path)
                    print(f"Saved {output_path}")

            except Exception as e:
                print(f"Failed to process {file_path}: {e}")


if __name__ == "__main__":
    process_directory(base_dir=".", subfolders=FOLDERS_TO_PROCESS)
