from pathlib import Path

BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")

# Get all .py files in Whole_mAb folder
antibody_files = [
    f for f in BASE_DIR.rglob("*.py")
    if "PDB_Output_Files_GPU" not in f.parts  # Skip output folder
    and not f.name.startswith("._")           # Skip hidden files
]

# Extract names
antibody_names = sorted([f.stem for f in antibody_files])

print(f"Found {len(antibody_names)} antibodies in Whole_mAb folder:")
print("-" * 60)
for name in antibody_names:
    print(name)

# Save to file
with open(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\whole_mab_antibody_list.txt", "w") as f:
    f.write("\n".join(antibody_names))
print(f"\n List saved to: whole_mab_antibody_list.txt")