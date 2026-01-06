import os
import re
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# --- CONFIG ---
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")

ANTIBODY_FOLDERS = [
    "Bispecific_mAb",
    "Bispecific_scFv",
    "Other_Formats",
    "Whole_mAb"
]

def extract_chains_and_count(file_path):
    """
    ---------------------------------------------------------------------------
    [REFERENCE NOTES: STOICHIOMETRY & ARCHITECTURE]
    ---------------------------------------------------------------------------
    WHOLE mAb:
    - Molar Ratio: 2 Heavy Chains : 2 Light Chains (2:2 or 1:1 Hc:Lc).
    - Mass Contribution: Two ~50 kDa heavy chains and two ~25 kDa light chains.
    - Assembly: Held by disulfide bonds, ensuring each heavy chain pairs with one light chain.
    - Function: This structure creates two identical antigen-binding sites (Fab regions). 

    IgG-LIKE BsAbs:
    - Stoichiometry: Four Chains. Often need two unique HCs and two unique LCs (HC1, HC2, LC1, LC2).

    COMMON BISPECIFIC scFv FORMATS:
    - Tandem scFv (taFv/scFv2): Two different scFvs connected by a linker (1:1 stoichiometry).
    - Diabody: Short linker favors inter-chain pairing; dimeric complex with two binding sites.
    - BiTEs: Links an anti-tumor antigen scFv and an anti-CD3 (T-cell) scFv.
    - scFv-IgG Formats: scFvs appended to or inserted into an intact IgG (e.g., 1+1, 2+2).

    KEY STRUCTURAL FEATURES:
    - Linker Length & Orientation: Critical for folding, stability, and preventing aggregation.
    - Variable Domains: VH and VL domains forming the minimal antigen-binding unit.
    - Engineered Disulfide Bonds: Cys residues added for stability; can sometimes cause aggregation.
    ---------------------------------------------------------------------------
    """
    found_seqs = []
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Count sequence blocks starting with >
        blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
        
        for block in blocks:
            lines = block.strip().split('\n')
            if len(lines) < 2: continue
            
            raw_seq = "".join(lines[1:]).strip().upper()
            clean_seq = re.sub(r'[^A-Z]', '', raw_seq)
            
            if clean_seq and clean_seq != "NA":
                found_seqs.append(clean_seq)
    except:
        pass
    return len(found_seqs)

def run_chain_count_summary():
    master_rows = []

    for folder_name in ANTIBODY_FOLDERS:
        folder_path = BASE_DIR / folder_name
        if not folder_path.exists(): continue
        
        py_files = list(folder_path.glob("*.py"))
        
        for py_path in tqdm(py_files, desc=f"Counting {folder_name}"):
            ab_name = py_path.stem
            
            # COUNT THE CHAINS IN THE .PY FILE
            chain_count = extract_chains_and_count(py_path)

            # MAP THE MULTIPLIER (TRUTH)
            # ONLY apply x2 if it's in the Whole_mAb folder AND has 2 chains.
            # Everything else (scFvs, Bispecifics, etc.) stays at 1.
            if folder_name == "Whole_mAb" and chain_count == 2:
                multiplier = 2
            else:
                multiplier = 1
            
            master_rows.append({
                "Antibody_Name": ab_name,
                "Folder": folder_name,
                "Chains_Counted": chain_count,
                "Biological_Multiplier": multiplier
            })

    # Output to Excel
    df = pd.DataFrame(master_rows)
    df.to_excel(BASE_DIR / "Antibody_Chain_Count_Summary.xlsx", index=False)
    print(f"\nDONE. Summary created at: {BASE_DIR} / Antibody_Chain_Count_Summary.xlsx")

if __name__ == "__main__":
    run_chain_count_summary()