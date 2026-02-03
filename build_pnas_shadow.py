import pandas as pd
from pathlib import Path

# =================================================================
# PROJECT: Shadow Builder (2025 Edition)
# PURPOSE: Convert Excel benchmarks into individual .py sequence files
#          Distributed across categorized therapeutic folders.
# =================================================================

# --- SETTINGS ---
BASE_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")
MASTER_FILE = BASE_DIR / "mAb_Truth_Engine_Master.xlsx"
SHADOW_BASE = BASE_DIR / "shadow_benchmarks"

def main():
    if not MASTER_FILE.exists():
        print(f"ERROR: Cannot find {MASTER_FILE}. Run mAb_truth_engine.py first!")
        return

    # 1. Load Data
    df = pd.read_excel(MASTER_FILE)
    
    # 2. MAPPING
    name_col = 'Name'
    vh_col = 'VH'
    vl_col = 'VL'
    folder_col = 'Target_Folder'

    # 3. Extraction Loop
    print(f"--- DISTRIBUTING {len(df)} BENCHMARKS TO CATEGORIZED SHADOW RIGS ---")
    
    count = 0
    for i, row in df.iterrows():
        # Sanitize Name
        raw_name = str(row[name_col]).strip()
        if raw_name == 'nan' or not raw_name:
            continue
            
        clean_name = "".join([c if c.isalnum() else "_" for c in raw_name])
        
        # Determine the correct sub-folder
        target_folder_name = str(row.get(folder_col, "Whole_mAb"))
        category_dir = SHADOW_BASE / target_folder_name
        
        if not category_dir.exists():
            category_dir.mkdir(parents=True, exist_ok=True)
        
        vh = str(row[vh_col]).strip()
        vl = str(row[vl_col]).strip()
        
        # 4. WRITE THE .PY FILE
        file_path = category_dir / f"{clean_name}.py"
        
        try:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(f'"""\n')
                f.write(f'FORMAT: {row["Format"]}\n') # Respects calculate_features requirement
                f.write(f'> {clean_name}_VH\n{vh}\n')
                f.write(f'> {clean_name}_VL\n{vl}\n')
                f.write(f'"""\n')
            count += 1 # Added to ensure success message is accurate
        except Exception as e:
            print(f"Skipping {clean_name} due to write error: {e}")

    print(f"\nSUCCESS: Distributed {count} shadow files across:")
    for folder in SHADOW_BASE.iterdir():
        if folder.is_dir():
            file_count = len(list(folder.glob('*.py')))
            print(f" - {folder.name}: {file_count} files")

    print(f"\nACTION: You can now run calculate_features.py on the '{SHADOW_BASE.name}' directory.")

if __name__ == "__main__":
    main()