from pathlib import Path

BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")

# glob instead of rglob to strictly capture the 300 monospecifics in this folder
antibody_files = [
    f for f in BASE_DIR.glob("*.py")
    if "PDB_Output_Files_GPU" not in f.parts
    and "PDB_Output_Files_GPU_Full" not in f.parts
    and not f.name.startswith("._")  
]

antibody_names = sorted([f.stem for f in antibody_files])

print(f"Found {len(antibody_names)} antibodies in Whole_mAb folder:")
print("-" * 60)
for name in antibody_names:
    print(name)

output_path = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\whole_mab_antibody_list.txt"
with open(output_path, "w") as f:
    f.write("\n".join(antibody_names))

print(f"\n Clean list of {len(antibody_names)} entries saved to: {output_path}")