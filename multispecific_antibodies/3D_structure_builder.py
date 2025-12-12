import os
import glob
from io import StringIO
from igfold import IgFoldRunner
from Bio import SeqIO

# --- Configuration Variables ---
# Use use_openmm=True for refinement, which prioritizes OpenMM over PyRosetta
USE_REFINEMENT = True
# Name of the main directory for PDB output files
OUTPUT_BASE_DIR = "PDB_Output_Files"
# Subfolders within the 'multispecific_antibodies' directory for processing
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Whole_mAb"]
# -------------------------------

# Initialize the IgFold runner without the 'device' argument
runner = IgFoldRunner()

def parse_fasta_string_to_pairs(fasta_string):
    """
    Parses a multi-FASTA string into a list of dictionaries,
    each representing an H/L pair for IgFold prediction.
    Chains are grouped based on the numeric identifier in the header (e.g., Chain_1).
    """
    chains_by_number = {}
    handle = StringIO(fasta_string)
    for record in SeqIO.parse(handle, "fasta"):
        seq_id = record.id
        seq = str(record.seq)

        if '_Heavy_Chain' in seq_id:
            chain_num = seq_id.split('_Heavy_Chain_')[-1]
            if chain_num not in chains_by_number:
                chains_by_number[chain_num] = {}
            chains_by_number[chain_num]['H'] = seq
        elif '_Light_Chain' in seq_id:
            chain_num = seq_id.split('_Light_Chain_')[-1]
            if chain_num not in chains_by_number:
                chains_by_number[chain_num] = {}
            chains_by_number[chain_num]['L'] = seq
        else:
            continue

    paired_structures = []
    for chain_num, chains in chains_by_number.items():
        if 'H' in chains and 'L' in chains:
            paired_structures.append({'H': chains['H'], 'L': chains['L']})
        else:
            continue
            
    return paired_structures

def process_directory(base_dir, subfolders):
    """Iterates through specified subfolders to process all antibody .py files."""
    for folder_name in subfolders:
        folder_path = os.path.join(base_dir, folder_name)
        output_folder_path = os.path.join(OUTPUT_BASE_DIR, folder_name)
        
        if not os.path.exists(folder_path):
            continue
            
        os.makedirs(output_folder_path, exist_ok=True)
        antibody_files = glob.glob(os.path.join(folder_path, '*.py'))

        for file_path in antibody_files:
            antibody_name, _ = os.path.splitext(os.path.basename(file_path))

            try:
                with open(file_path, 'r') as f:
                    file_content = f.read()
                
                exec(file_content, globals())
                sequence_string = globals()[antibody_name]

                paired_sequences_list = parse_fasta_string_to_pairs(sequence_string)

                for i, sequences_dict in enumerate(paired_sequences_list):
                    pair_id = f"{antibody_name}_Pair_{i+1}"
                    
                    # Call predict without 'device' argument
                    predicted_structure = runner.predict(
                        sequences=sequences_dict,
                        do_refine=USE_REFINEMENT,
                    )
                    
                    output_path = os.path.join(output_folder_path, f"{pair_id}.pdb")
                    predicted_structure.save(output_path)

            except Exception as e:
                continue

if __name__ == "__main__":
    process_directory(base_dir=".", subfolders=FOLDERS_TO_PROCESS)