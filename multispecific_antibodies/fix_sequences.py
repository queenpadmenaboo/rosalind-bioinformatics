"""
Fixes antibody .py files by replacing sequences with CSV sequences.

PURPOSE:
    - Finds antibody files with sequence mismatches
    - Replaces .py sequences with correct sequences from CSV
    - Creates backup of original files before fixing

USAGE:
    python fix_sequences.py
"""

import csv
import shutil
from pathlib import Path
from typing import Dict

# ============================
# CONFIGURATION
# ============================

CSV_PATH = Path(r"C:\Users\meeko\TheraSAbDab_SeqStruc_08Dec2025.csv")
ANTIBODY_FOLDER = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")

CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

# ============================
# FUNCTIONS
# ============================

def load_csv_data(csv_file: Path) -> Dict[str, Dict]:
    # Reads TheraSAbDab CSV and returns dict of antibody data
    antibodies = {}
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            therapeutic_raw = row['Therapeutic'].strip()
            key = therapeutic_raw.lower().replace(' ', '_').replace('-', '_')
            antibodies[key] = {
                'name': therapeutic_raw,
                'heavy1': row['HeavySequence'].strip(),
                'light1': row['LightSequence'].strip(),
                'heavy2': row['HeavySequence(ifbispec)'].strip(),
                'light2': row['LightSequence(ifbispec)'].strip(),
            }
    return antibodies


def find_antibody_file(antibody_name: str) -> Path:
    # Finds the .py file for an antibody in any of the 4 folders
    for folder in CATEGORY_FOLDERS:
        filepath = ANTIBODY_FOLDER / folder / f"{antibody_name}.py"
        if filepath.exists():
            return filepath
    # Check root folder
    filepath = ANTIBODY_FOLDER / f"{antibody_name}.py"
    if filepath.exists():
        return filepath
    return None


def fix_antibody_file(filepath: Path, csv_data: Dict) -> bool:
    # Replaces sequences in .py file with CSV sequences
    # Returns True if fixed, False if error
    
    antibody_name = filepath.stem
    
    # Create backup
    backup_path = filepath.with_suffix('.py.bak')
    shutil.copy(filepath, backup_path)
    
    # Build new file content
    # Format: antibody_name = """
    # >Header_Heavy_Chain_1
    # SEQUENCE
    # >Header_Light_Chain_1
    # SEQUENCE
    # """
    
    heavy1 = csv_data['heavy1'] if csv_data['heavy1'] else 'na'
    light1 = csv_data['light1'] if csv_data['light1'] else 'na'
    heavy2 = csv_data['heavy2'] if csv_data['heavy2'] else ''
    light2 = csv_data['light2'] if csv_data['light2'] else ''
    
    # Build FASTA content
    lines = []
    lines.append(f'{antibody_name} = """')
    
    # Heavy chain 1
    lines.append(f'>{antibody_name.capitalize()}_Heavy_Chain_1')
    lines.append(heavy1)
    
    # Light chain 1
    lines.append(f'>{antibody_name.capitalize()}_Light_Chain_1')
    lines.append(light1)
    
    # Heavy chain 2 (bispecifics)
    if heavy2:
        lines.append(f'>{antibody_name.capitalize()}_Heavy_Chain_2')
        lines.append(heavy2)
    
    # Light chain 2 (bispecifics)
    if light2:
        lines.append(f'>{antibody_name.capitalize()}_Light_Chain_2')
        lines.append(light2)
    
    lines.append('"""')
    
    # Write new content
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    return True


def main():
    print("=" * 70)
    print("FIX ANTIBODY SEQUENCES")
    print("=" * 70)
    
    # List of antibodies to fix (from validation report)
    antibodies_to_fix = [
        'belimumab', 'obinutuzumab', 'trastuzumab', 'acasunlimab',
        'amivantamab', 'amostomig', 'arumakimig', 'besufetamig',
        'cevostamab', 'ciduvectamig', 'davutamig', 'denecimig',
        'dilpacimab', 'elranatamab', 'etentamig', 'etuptamig',
        'fazpilodemab', 'gremubamab', 'ivicentamab', 'linclatamig',
        'linvoseltamab', 'marlotamig', 'navicixizumab', 'blinatumomab'
    ]
    
    print(f"Antibodies to fix: {len(antibodies_to_fix)}\n")
    
    # Load CSV data
    csv_antibodies = load_csv_data(CSV_PATH)
    print(f"Loaded {len(csv_antibodies)} antibodies from CSV\n")
    
    # Confirm
    print("This will:")
    print("  1. Create .bak backup of each file")
    print("  2. Replace sequences with CSV sequences")
    print()
    response = input("Proceed? (yes/no): ").strip().lower()
    if response != 'yes':
        print("Cancelled.")
        return
    
    print("\n" + "-" * 70)
    
    fixed = 0
    errors = 0
    
    for antibody_name in antibodies_to_fix:
        filepath = find_antibody_file(antibody_name)
        
        if not filepath:
            print(f"  ERROR: File not found for {antibody_name}")
            errors += 1
            continue
        
        if antibody_name not in csv_antibodies:
            print(f"  ERROR: {antibody_name} not in CSV")
            errors += 1
            continue
        
        csv_data = csv_antibodies[antibody_name]
        
        try:
            fix_antibody_file(filepath, csv_data)
            print(f"  Fixed: {antibody_name}")
            fixed += 1
        except Exception as e:
            print(f"  ERROR fixing {antibody_name}: {e}")
            errors += 1
    
    print("\n" + "=" * 70)
    print(f"Fixed: {fixed}")
    print(f"Errors: {errors}")
    print("=" * 70)
    print("\nRun validate_antibody_sequences.py to verify fixes.")


if __name__ == "__main__":
    main()