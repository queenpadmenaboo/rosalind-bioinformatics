import os
import pandas as pd
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =================================================================
# PROJECT: mAb Truth Engine
# SOURCE: Thera-SAbDab Master (Oxford Protein Informatics Group)
# =================================================================

def calculate_biophysics(vh_seq, vl_seq):
    """Calculates theoretical pI and GRAVY for the variable domains."""
    try:
        if pd.isna(vh_seq) or pd.isna(vl_seq):
            return None, None
            
        # Combine sequences for the Fv complex
        full_seq = str(vh_seq).strip().upper() + str(vl_seq).strip().upper()
        # Clean: BioPython fails on non-standard AA or whitespace
        clean = ''.join(aa for aa in full_seq if aa in "ACDEFGHIKLMNPQRSTVWY")
        
        if len(clean) < 50:
            return None, None
            
        analysis = ProteinAnalysis(clean)
        return round(analysis.isoelectric_point(), 2), round(analysis.gravy(), 3)
    except Exception:
        return None, None

def run_truth_engine():
    print("--- STARTING mAb TRUTH ENGINE (SAbDab 2025 SYNC) ---")

    # --- UPDATED PATHS TO MATCH YOUR DEC 2025 DOWNLOAD ---
    base_path = r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"
    # Pointing exactly to where your file is located
    SABDAB_FILE = r"C:\Users\meeko\TheraSAbDab_SeqStruc_08Dec2025.csv" 

    if not os.path.exists(SABDAB_FILE):
        print(f">>> [CRITICAL] Could not find {SABDAB_FILE}. Please check file path.")
        return

    try:
        # Thera-SAbDab often uses tabs or commas; engine='python' handles formatting quirks
        df = pd.read_csv(SABDAB_FILE, engine='python')
        print(f">>> [SUCCESS] Loaded {len(df)} therapeutics from SAbDab.")
    except Exception as e:
        print(f">>> [ERROR] CSV Parse failed: {e}")
        return

    # 1. 2025 STATUS LOGIC (Dynamic Mapping)
    # Mapping the Feb '25 headers found in Oxford's latest exports
    trial_col = next((c for c in df.columns if 'Highest Clinical Trial' in c), 'Highest_Clinical_Trial')
    status_col = next((c for c in df.columns if 'Estimated Status' in c), 'Estimated_Status')
    format_col = next((c for c in df.columns if 'Format' in c), 'Format')
    
    def assign_2025_status(row):
        trial = str(row.get(trial_col, "")).upper()
        est = str(row.get(status_col, "")).upper()
        
        if "APPROV" in trial:
            return "SUCCESS_APPROVED"
        if "DISCON" in est or "WITHDRAWN" in est:
            return "FAILURE_DISCONTINUED"
        if "PHASE" in trial:
            return f"ACTIVE_{trial}"
        return "INVESTIGATIONAL"

    df['Final_Status_2025'] = df.apply(assign_2025_status, axis=1)

    # --- NEW: MULTI-FOLDER CATEGORIZATION ---
    def assign_target_folder(row):
        fmt = str(row.get(format_col, "")).lower()
    
        # 1. Catch Homodimers/Tandems first - route to Other_Formats
        if any(k in fmt for k in ['homodimer', 'tandem', 'dual-variable']):
            return 'Other_Formats' 
        
        # 2. Standard Bispecifics
        if 'bispecific' in fmt or 'bsab' in fmt:
            return 'Bispecific_mAb'
        
        # 3. Standard Fragments
        if 'scfv' in fmt:
            return 'Bispecific_scFv'
    
        # 4. Catch-all for other fragments
        if any(k in fmt for k in ['vhh', 'nanobody', 'fab', 'fragment']):
            return 'Other_Formats'

        # 5. Default
        return 'Whole_mAb'

    df['Target_Folder'] = df.apply(assign_target_folder, axis=1)

    # 2. SEQUENCE CLEANING
    vh_col = next((c for c in df.columns if 'Heavy Chain Sequence' in c or 'VH' in c.upper() or 'HEAVY' in c.upper()), 'Heavy Sequence')
    vl_col = next((c for c in df.columns if 'Light Chain Sequence' in c or 'VL' in c.upper() or 'LIGHT' in c.upper()), 'Light Sequence')

    print(f">>> Calculating biophysics using BioPython...")
    results = df.apply(lambda x: calculate_biophysics(x[vh_col], x[vl_col]), axis=1)
    df[['Calculated_pI', 'Calculated_GRAVY']] = pd.DataFrame(results.tolist(), index=df.index)

    # 3. STANDARDIZATION FOR PIPELINE
    # Renaming so 'Name', 'VH', and 'VL' are consistent for build_pnas_shadow.py
    main_name_col = next((c for c in df.columns if 'Therapeutic' in c or 'name' in c.lower()), 'Therapeutic Name')
    df = df.rename(columns={main_name_col: 'Name', vh_col: 'VH', vl_col: 'VL'})

    # 4. EXPORT MASTER
    output_path = os.path.join(base_path, "mAb_Truth_Engine_Master.xlsx")
    df.to_excel(output_path, index=False)
    
    print("-" * 50)
    print(f"TRUTH ENGINE SYNCED TO 2025")
    print(f"Total Benchmarks: {len(df)}")
    print(f"File Saved: {output_path}")
    print(f"Folders Mapped: {df['Target_Folder'].unique()}")
    print("-" * 50)

if __name__ == "__main__":
    run_truth_engine()