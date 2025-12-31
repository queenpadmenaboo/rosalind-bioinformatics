import pandas as pd
from pathlib import Path

# --- PATHS ---
BASE_DIR = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")
PNAS_MASTER = BASE_DIR / "mAb_Truth_Engine_Master.xlsx"
YOUR_PREDS = BASE_DIR / "shadow_benchmarks" / "sequence_features.xlsx"

def run_pnas_validation():
    print("--- PNAS BENCHMARK VALIDATION ENGINE (2025 STATUS READY) ---")
    
    if not PNAS_MASTER.exists() or not YOUR_PREDS.exists():
        print("ERROR: Missing files. Run your other scripts first!")
        return

    # 1. Load Data
    df_lab = pd.read_excel(PNAS_MASTER)
    df_calc = pd.read_excel(YOUR_PREDS)
    
    # 2. Merge (The Handshake)
    # We use 'inner' to only see antibodies where we have both sequence and clinical data
    df_merged = pd.merge(
        df_calc, 
        df_lab, 
        left_on='antibody', 
        right_on='Name', 
        how='inner',
        suffixes=('_MINE', '_PNAS')
    )
    
    # 3. Dynamic Column Identification
    # SAbDab/PNAS often have slightly different naming for trial status
    status_col = next((c for c in df_merged.columns if 'Final_Status_2025' in c), 'Final_Status_2025')
    # Try to find the Experimental Viscosity column from PNAS dataset
    visc_col = next((c for c in df_merged.columns if 'Viscosity' in str(c) and '150' in str(c)), None)

    # 4. Cleanup & Comparison Logic
    if visc_col:
        # Sort so the highest viscosity (worst performers) are at the top
        df_merged = df_merged.sort_values(by=visc_col, ascending=False)
    
    def check_theory_alignment(row):
        calc_risk = str(row.get('viscosity_risk', ''))
        clinical = str(row.get(status_col, '')).upper()
        
        # THEORY: Approved drugs should generally be "Low Risk" (Charge > 20)
        if "APPROV" in clinical and "Low Risk" in calc_risk:
            return "MATCH: Validated Theory"
        if "DISCON" in clinical and "High Risk" in calc_risk:
            return "MATCH: Failure Predicted"
        return "MISMATCH / INVESTIGATIONAL"

    df_merged['Theory_Validation'] = df_merged.apply(check_theory_alignment, axis=1)

    # 5. Save Final Report
    output_path = BASE_DIR / "PNAS_VS_CALCULATIONS.xlsx"
    
    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        df_merged.to_excel(writer, index=False, sheet_name='Comparison')
        wb, ws = writer.book, writer.sheets['Comparison']
        
        # --- FORMATTING ---
        # 1. Approved/Validated = Green
        green_fmt = wb.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})
        # 2. High Viscosity Risk = Red
        red_fmt = wb.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})

        # Apply Highlight to Theory Validation
        theory_idx = df_merged.columns.get_loc('Theory_Validation')
        ws.conditional_format(1, theory_idx, len(df_merged), theory_idx, {
            'type': 'text', 'criteria': 'containing', 'value': 'MATCH', 'format': green_fmt
        })

        # Highlight calculated High Risk
        risk_idx = df_merged.columns.get_loc('viscosity_risk')
        ws.conditional_format(1, risk_idx, len(df_merged), risk_idx, {
            'type': 'text', 'criteria': 'containing', 'value': 'High Risk', 'format': red_fmt
        })

        # Freeze Headers & Add Filters
        ws.autofilter(0, 0, len(df_merged), len(df_merged.columns) - 1)
        ws.freeze_panes(1, 0)
        
        # Auto-fit column widths
        for i, col in enumerate(df_merged.columns):
            ws.set_column(i, i, max(df_merged[col].astype(str).map(len).max(), len(col)) + 2)

    print(f"\nSUCCESS: Comparison complete with 2025 status highlighting.")
    print(f"Total Matches Analyzed: {len(df_merged)}")
    print(f"Final Report: {output_path.name}")

if __name__ == "__main__":
    run_pnas_validation()