import pandas as pd
import re
from pathlib import Path

# 1. SCAN THE ACTUAL FOLDER INSTEAD OF THE TXT FILE
# This ensures it only looks at the 305 files actually sitting in the directory
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")
folder_names = [f.stem.lower() for f in BASE_DIR.glob("*.py")]

print(f"Scanning folder: Found {len(folder_names)} actual files.")

# 2. Load CSV and Match
df = pd.read_csv(r"C:\Users\bunsr\TheraSAbDAb_SeqStruc_07Dec2025.csv")
matched = df[df['Therapeutic'].str.lower().isin(folder_names)]

print(f"Matched {len(matched)}/{len(folder_names)} antibodies in folder to CSV records")

# Match (case-insensitive)
matched = df[df['Therapeutic'].str.lower().isin(folder_names)]

print(f"Matched {len(matched)}/{len(folder_names)} antibodies to CSV")
print("\n" + "="*60)
print("CH1 ISOTYPE DISTRIBUTION:")
print("="*60)
print(matched['CH1 Isotype'].value_counts())
print("\n" + "="*60)
print("LIGHT CHAIN DISTRIBUTION:")
print("="*60)
print(matched['VD LC'].value_counts())
print("\n" + "="*60)
print("FORMAT DISTRIBUTION:")
print("="*60)
print(matched['Format'].value_counts())

# Show any antibodies NOT matched
not_matched = set(folder_names) - set(matched['Therapeutic'].str.lower())
if not_matched:
    print("\n" + "="*60)
    print(f"{len(not_matched)} antibodies in folder NOT found in CSV:")
    print("="*60)
    for name in sorted(not_matched):
        print(f"  - {name}")

# ============================================================
# IDENTIFY NON-HUMAN ANTIBODIES
# ============================================================
print("\n" + "="*60)
print("IDENTIFYING NON-HUMAN ANTIBODIES:")
print("="*60)

non_human = matched[
    (matched['Format'].str.contains('Mouse', na=False)) |
    (matched['Format'].str.contains('Canine', na=False)) |
    (matched['CH1 Isotype'].str.contains('G2a', na=False)) |
    (matched['CH1 Isotype'].str.contains('G2b', na=False))
]

print(f"\nFound {len(non_human)} non-human antibodies:")
if len(non_human) > 0:
    print(non_human[['Therapeutic', 'Format', 'CH1 Isotype']].to_string(index=False))
    
    # Save list of non-human antibodies to remove
    non_human_names = non_human['Therapeutic'].str.lower().tolist()
    with open(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\non_human_antibodies_to_remove.txt", "w") as f:
        f.write("\n".join(non_human_names))
    
    print(f"\n✅ Saved to: non_human_antibodies_to_remove.txt")
else:
    print("✅ No non-human antibodies found!")

# ============================================================
# CHECK BISPECIFIC/FUSION FILES FOR ACTUAL SEQUENCE COUNTS
# ============================================================
print("\n" + "="*60)
print("CHECKING BISPECIFIC/FUSION FILES:")
print("="*60)

BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")

# Get antibodies marked as "Bispecific Whole mAb" or "Fusion" in CSV
bispec_or_fusion = df[
    (df['Format'].str.contains('Bispecific', na=False)) |
    (df['Format'].str.contains('Fusion', na=False))
]['Therapeutic'].str.lower().tolist()

print(f"\nChecking {len(bispec_or_fusion)} files labeled as Bispecific/Fusion:")

for name in bispec_or_fusion:
    file_path = BASE_DIR / f"{name}.py"
    
    if not file_path.exists():
        continue
    
    # Count sequences
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    blocks = content.split('>')
    sequences = []
    
    for block in blocks[1:]:
        lines = block.strip().split('\n')
        if len(lines) < 2:
            continue
        
        header = lines[0].upper()
        raw_data = "".join(lines[1:]).strip()
        clean_seq = re.sub(r'[^A-Z]', '', raw_data.upper())
        
        if clean_seq and clean_seq != "NA" and len(clean_seq) > 20:
            chain_type = 'H' if 'HEAVY' in header else 'L'
            sequences.append(chain_type)
    
    heavy_count = sequences.count('H')
    light_count = sequences.count('L')
    
    csv_format = df[df['Therapeutic'].str.lower() == name]['Format'].values[0]
    
    print(f"\n{name}:")
    print(f"  CSV Format: {csv_format}")
    print(f"  Sequences in file: {heavy_count}H + {light_count}L")
    
    if heavy_count == 1 and light_count == 1:
        print(f"  KEEP - Standard whole mAb structure (will be duplicated to 2H:2L)")
    elif heavy_count == 2 and light_count == 2:
        print(f"  MOVE to Bispecific_mAb/ - Has 2 different H+L pairs")
    else:
        print(f"  UNUSUAL - Review manually")

# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "="*60)
print("SUMMARY:")
print("="*60)
print(f"Total antibodies: {len(matched)}")
print(f"Non-human (should remove): {len(non_human)}")
print(f"Bispecific/Fusion to review: {len(bispec_or_fusion)}")
print(f"Clean human dataset: {len(matched) - len(non_human)}")