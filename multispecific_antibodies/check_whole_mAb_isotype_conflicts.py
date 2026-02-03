import pandas as pd
from pathlib import Path

# --- DIRECT PATHS ---
CSV_PATH = Path(r"C:\Users\bunsr\TheraSAbDab_SeqStruc_07Dec2025.csv")
MAB_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")

def run_isotype_inventory():
    if not CSV_PATH.exists():
        print(f"ERROR: CSV not found at {CSV_PATH}")
        return
    
    # 1. Load Master Data
    df = pd.read_csv(CSV_PATH)
    df['match_name'] = df['Therapeutic'].str.strip().str.lower()

    # 2. Get EVERY .py file in the folder (No exclusions)
    local_files = {f.stem.lower(): f.name for f in MAB_DIR.glob("*.py")}
    local_names = set(local_files.keys())
    
    # 3. Match against CSV
    inventory = df[df['match_name'].isin(local_names)].copy()
    matched_names = set(inventory['match_name'])
    
    # 4. Identify Flags (Files in folder NOT in CSV)
    ghost_files = local_names - matched_names

    print(f"\n--- AUDIT RESULTS: {len(local_names)} FILES FOUND ---")
    
    if ghost_files:
        print(f"\n[FLAG] {len(ghost_files)} FILES NOT MATCHING CSV DATA:")
        for ghost in sorted(ghost_files):
            print(f"  - {local_files[ghost]}")
    
    print("-" * 50)
    print("CURRENT FOLDER DATA (FROM CSV):")
    
    if not inventory.empty:
        # 5. Output raw counts and the full data table
        print("\nHEAVY CHAIN ISOTYPES:")
        print(inventory['CH1 Isotype'].value_counts())
        
        print("\nLIGHT CHAIN TYPES:")
        print(inventory['VD LC'].value_counts())

        print("\nFULL DATA MAP:")
        print(inventory[['Therapeutic', 'CH1 Isotype', 'VD LC']].sort_values('Therapeutic').to_string(index=False))
        
        # Save the master sync file
        inventory[['Therapeutic', 'CH1 Isotype', 'VD LC']].to_csv(MAB_DIR / "mAb_Data_Sync_Master.csv", index=False)
    else:
        print("No matching records found in CSV for files in this directory.")

if __name__ == "__main__":
    run_isotype_inventory()
 
import os
from pathlib import Path

MAB_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb")

def find_hidden_files():
    # Get everything in the folder
    all_files = list(MAB_DIR.iterdir())
    py_files = list(MAB_DIR.glob("*.py"))
    
    print(f"Total items in folder: {len(all_files)}")
    print(f"Total .py files found: {len(py_files)}")
    
    if len(all_files) != len(py_files):
        print("\nNON-PYTHON FILES DETECTED:")
        for f in all_files:
            if f.suffix != '.py':
                print(f"  - {f.name} ({f.suffix})")

if __name__ == "__main__":
    find_hidden_files()