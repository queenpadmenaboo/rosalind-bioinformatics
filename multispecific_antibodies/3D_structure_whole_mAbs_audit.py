import os
from pathlib import Path
from Bio.PDB import PDBParser

# Path to your IgFold/Theoretical PDBs
PDB_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb\PDB_Output_Files_GPU_Full")
OUTPUT_FILE = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb\Whole_mAb_Structure_Audit.txt")

def audit_pdb_contents():
    parser = PDBParser(QUIET=True)

    # Check if directory exists first
    if not PDB_DIR.exists():
        print(f"ERROR: Path does not exist: {PDB_DIR}")
        return
    
    pdb_files = list(PDB_DIR.rglob("*.pdb"))
    
    if not pdb_files:
        print(f"Zero PDBs found in {PDB_DIR}")
        return
    
    print(f"Found {len(pdb_files)} PDB files. Auditing structure...")
    
    with open(OUTPUT_FILE, "w") as f:
        header = f"{'Filename':<35} | {'Chains':<15} | {'Residues':<10}"
        f.write("FULL MAB STRUCTURE AUDIT REPORT\n" + "="*70 + "\n")
        f.write(header + "\n" + "-"*70 + "\n")
        
        print(header)
        print("-" * 60)

        for pdb_path in pdb_files: 
            try:
                structure = parser.get_structure("check", str(pdb_path))
                chains = [chain.id for chain in structure.get_chains()]
                residue_count = len(list(structure.get_residues()))
                
                result_line = f"{pdb_path.name:<35} | {str(chains):<15} | {residue_count:<10}"
                
                # Print to terminal AND write to file
                print(result_line)
                f.write(result_line + "\n")
                
            except Exception as e:
                err_msg = f"Error reading {pdb_path.name}: {e}"
                print(err_msg)
                f.write(err_msg + "\n")

    print(f"\n=== DONE! Full audit saved to: {OUTPUT_FILE.name} ===")

if __name__ == "__main__":
    audit_pdb_contents()