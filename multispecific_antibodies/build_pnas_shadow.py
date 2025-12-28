import pandas as pd
from pathlib import Path

# --- SETTINGS ---
BASE_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")
MASTER_FILE = BASE_DIR / "mAb_Truth_Engine_Master.xlsx"
SHADOW_DIR = BASE_DIR / "shadow_benchmarks" / "Whole_mAb"

def main():
    if not MASTER_FILE.exists():
        print(f"ERROR: Cannot find {MASTER_FILE}")
        return

    # 1. Load Data
    df = pd.read_excel(MASTER_FILE)
    
    # 2. Setup Shadow Dir
    if not SHADOW_DIR.exists():
        SHADOW_DIR.mkdir(parents=True)
    
    # 3. DIRECT MAPPING (Based on your terminal output)
    name_col = 'Name'
    vh_col = 'VH'
    vl_col = 'VL'

    # 4. Extraction Loop
    print(f"--- EXTRACTING {len(df)} BENCHMARKS TO SHADOW RIG ---")
    
    count = 0
    for i, row in df.iterrows():
        # Sanitize Name
        raw_name = str(row[name_col]).strip()
        if raw_name == 'nan' or not raw_name:
            continue
            
        # Clean filename: only allow alphanumeric and underscores
        clean_name = "".join([c if c.isalnum() else "_" for c in raw_name])
        
        vh = str(row[vh_col]).strip()
        vl = str(row[vl_col]).strip()
        
        # Write the .py file in the format your scripts love
        file_path = SHADOW_DIR / f"{clean_name}.py"
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(f'""" PNAS BENCHMARK: {clean_name} """\n')
            f.write(f'> {clean_name}_H\n{vh}\n> {clean_name}_L\n{vl}\n')
        count += 1

    print(f"\nSUCCESS: Created {count} files in {SHADOW_DIR}")

if __name__ == "__main__":
    main()