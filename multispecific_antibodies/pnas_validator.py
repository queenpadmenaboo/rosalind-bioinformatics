import pandas as pd
from pathlib import Path

# --- PATHS ---
BASE_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")
PNAS_MASTER = BASE_DIR / "mAb_Truth_Engine_Master.xlsx"
YOUR_PREDS = BASE_DIR / "shadow_benchmarks" / "sequence_features.xlsx"

def run_pnas_validation():
    print("--- PNAS BENCHMARK VALIDATION ENGINE ---")
    
    if not PNAS_MASTER.exists() or not YOUR_PREDS.exists():
        print("ERROR: Missing Master or Prediction file. Run your other scripts first!")
        return

    # 1. Load Data
    df_lab = pd.read_excel(PNAS_MASTER)
    df_calc = pd.read_excel(YOUR_PREDS)
    
    # 2. Merge (The Handshake)
    # We use 'antibody' from your calc and 'Name' from PNAS
    df_merged = pd.merge(
        df_calc, 
        df_lab, 
        left_on='antibody', 
        right_on='Name', 
        how='inner',
        suffixes=('_MINE', '_PNAS')
    )
    
    # 3. Cleanup & Logic Check
    # Locate the actual Viscosity column from PNAS (usually contains '150')
    visc_col = next((c for c in df_merged.columns if 'Viscosity' in c and '150' in str(c)), None)
    
    if visc_col:
        # Sort by the thickest viscosity to see the worst offenders
        df_merged = df_merged.sort_values(by=visc_col, ascending=False)
        
        # Add a quick Theory Match column
        # Theory: Absolute Net Charge > 20 should result in Low Viscosity
        def validate_theory(row):
            if abs(row['net_charge_pH5.5']) > 20 and row[visc_col] < 20:
                return "Theory Matches"
            elif abs(row['net_charge_pH5.5']) < 10 and row[visc_col] > 20:
                return "Theory Matches (Sticky)"
            else:
                return "Outlier (Check Hydrophobicity)"

        df_merged['PNAS_Validation_Verdict'] = df_merged.apply(validate_theory, axis=1)

    # 4. Save Final Report
    output_path = BASE_DIR / "PNAS_VS_CALCULATIONS.xlsx"
    df_merged.to_excel(output_path, index=False)
    
    print(f"\nSUCCESS: Comparison complete.")
    print(f"Mapped: {len(df_merged)} antibodies.")
    print(f"Final Report: {output_path.name}")

if __name__ == "__main__":
    run_pnas_validation()