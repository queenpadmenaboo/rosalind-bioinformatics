import os
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =================================================================
# PROJECT: mAb Truth Engine
# DIRECTORY: C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies
# =================================================================

def calculate_biophysics(vh_seq, vl_seq):
    """Calculates theoretical pI and GRAVY for the IgG model (2H + 2L)."""
    try:
        if pd.isna(vh_seq) or pd.isna(vl_seq):
            return None, None
            
        # Construct the full protein complex (2 Heavy + 2 Light)
        full_seq = (str(vh_seq).strip() * 2) + (str(vl_seq).strip() * 2)
        
        # Clean sequence: Remove whitespace and non-standard AA
        clean = ''.join(aa for aa in full_seq.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
        
        if len(clean) < 50:
            return None, None
            
        analysis = ProteinAnalysis(clean)
        return round(analysis.isoelectric_point(), 2), round(analysis.gravy(), 3)
    except Exception:
        return None, None

def run_truth_engine():
    print("--- STARTING mAb TRUTH ENGINE ---")

    base_path = r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"
    
    FILE_S1 = os.path.join(base_path, "pnas.1616408114.sd01.xlsx")
    FILE_S2 = os.path.join(base_path, "pnas.1616408114.sd02.xlsx")
    FILE_S3 = os.path.join(base_path, "pnas.1616408114.sd03.xlsx")

    try:
        print(">>> Attempting handshake with PNAS supplements...")
        # skiprows=1 is often needed for PNAS files as they have a title row
        df_list = pd.read_excel(FILE_S1, engine='openpyxl')
        df_seqs = pd.read_excel(FILE_S2, engine='openpyxl')
        df_labs = pd.read_excel(FILE_S3, engine='openpyxl')
        print(">>> [SUCCESS] All 3 Excel sheets parsed.")
    except Exception as e:
        print(f">>> [CRITICAL ERROR] Handshake failed: {e}")
        return

    # 3. DATA ALIGNMENT
    for df in [df_list, df_seqs, df_labs]:
        name_col = next((c for c in df.columns if 'name' in c.lower()), None)
        if name_col:
            df['join_key'] = df[name_col].astype(str).str.lower().str.strip()

    # Merge the 3 datasets
    df_merged = pd.merge(df_seqs, df_labs, on='join_key', suffixes=('', '_lab'))
    df_final = pd.merge(df_merged, df_list, on='join_key', suffixes=('', '_list'))

    # 4. EXECUTE CALCULATIONS
    vh_col = next((c for c in df_final.columns if 'VH' in c.upper()), None)
    vl_col = next((c for c in df_final.columns if 'VL' in c.upper()), None)

    if vh_col and vl_col:
        print(f">>> Calculating biophysics for {len(df_final)} antibodies...")
        results = df_final.apply(lambda x: calculate_biophysics(x[vh_col], x[vl_col]), axis=1)
        df_final[['Calculated_pI', 'Calculated_GRAVY']] = pd.DataFrame(results.tolist(), index=df_final.index)
        
        # --- STEP 4.5: RENAME FOR DOWNSTREAM PIPELINE ---
        # This ensures 'Name', 'VH', and 'VL' exist for the shadow builder
        main_name_col = next((c for c in df_final.columns if 'name' in c.lower()), 'Antibody Name')
        df_final = df_final.rename(columns={main_name_col: 'Name', vh_col: 'VH', vl_col: 'VL'})
    else:
        print("[ERROR] Could not find VH/VL sequence columns.")

    # 5. EXPORT
    output_path = os.path.join(base_path, "mAb_Truth_Engine_Master.xlsx")
    df_final.to_excel(output_path, index=False)
    
    print("-" * 50)
    print(f"TRUTH ENGINE COMPLETE")
    print(f"File Saved: {output_path}")
    print("-" * 50)

if __name__ == "__main__":
    run_truth_engine()