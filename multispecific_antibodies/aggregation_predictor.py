import os
import glob
import importlib.util
from pathlib import Path
from io import StringIO
import numpy as np
import pandas as pd
from Bio import SeqIO

# --- Biophysical Propensity Scales ---
HYDRO_SCALE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4, 'H': -3.2, 
    'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 
    'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}
BETA_SCALE = {
    'A': 0.83, 'C': 1.19, 'D': 0.54, 'E': 0.37, 'F': 1.38, 'G': 0.75, 'H': 0.87, 
    'I': 1.60, 'K': 0.74, 'L': 1.30, 'M': 1.05, 'N': 0.89, 'P': 0.55, 'Q': 1.10, 
    'R': 0.93, 'S': 0.75, 'T': 1.19, 'V': 1.70, 'W': 1.37, 'Y': 1.47
}

FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]
WINDOW = 7

def analyze_sequence(sequence):
    h_scores, b_scores = [], []
    apr_count = 0
    for i in range(len(sequence) - WINDOW + 1):
        sub = sequence[i : i + WINDOW]
        h = sum(HYDRO_SCALE.get(aa, 0) for aa in sub) / WINDOW
        b = sum(BETA_SCALE.get(aa, 0) for aa in sub) / WINDOW
        h_scores.append(h)
        b_scores.append(b)
        if h > 2.0 and b > 1.2:
            apr_count += 1
    avg_beta = np.mean(b_scores) if b_scores else 0
    avg_hydro = np.mean(h_scores) if h_scores else 0
    overall = (apr_count * 10) + (avg_beta * 20) + (avg_hydro * 5)
    return apr_count, round(avg_beta, 3), round(overall, 2)

def run_predictor():
    results = []
    base_path = Path(".")

    for folder in FOLDERS_TO_PROCESS:
        folder_path = base_path / folder
        if not folder_path.exists(): continue
        
        py_files = glob.glob(str(folder_path / "*.py"))
        for file_path in py_files:
            file_path = Path(file_path)
            ab_name = file_path.stem
            try:
                spec = importlib.util.spec_from_file_location(ab_name, file_path)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                if hasattr(module, ab_name):
                    fasta_string = getattr(module, ab_name)
                    handle = StringIO(fasta_string.strip())
                    for record in SeqIO.parse(handle, "fasta"):
                        apr, beta, overall = analyze_sequence(str(record.seq))
                        results.append({
                            "Antibody": ab_name,
                            "Chain": record.id,
                            "Folder": folder,
                            "APR_Count": apr,
                            "Beta_Propensity": beta,
                            "Overall_Score": overall
                        })
            except Exception as e:
                print(f"Error in {ab_name}: {e}")

    df = pd.DataFrame(results)
    if not df.empty:
        # Sort by score: Riskiest candidates at the top
        df = df.sort_values(by="Overall_Score", ascending=False)
        
        output_file = "Aggregation_Risk_Report.xlsx"
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Aggregation Risk', index=False)

        workbook  = writer.book
        worksheet = writer.sheets['Aggregation Risk']

        # 1. ADD FILTER BUTTONS TO THE TOP ROW
        # Defines the range for the filter from first cell (0,0) to the end of data
        worksheet.autofilter(0, 0, len(df), len(df.columns) - 1)

        # 2. FREEZE THE TOP ROW (Headers stay visible when scrolling)
        worksheet.freeze_panes(1, 0)

        # 3. AUTO-FIT COLUMN WIDTHS
        for i, col in enumerate(df.columns):
            # Calculate width based on max string length in column or header name
            max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, max_len)

        writer.close()
        print(f"SUCCESS: {output_file} created with filter buttons and frozen header.")
    else:
        print("No antibody sequences found to process.")

if __name__ == "__main__":
    run_predictor()