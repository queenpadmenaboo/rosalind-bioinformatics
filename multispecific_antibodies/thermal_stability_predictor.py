import os
import glob
import importlib.util
from pathlib import Path
from io import StringIO
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

"""
STABILITY ANALYSIS & BIOPHYSICAL THEORY:

1. DEFINITION OF ALIPHATIC: 
   Refers to non-polar, 'oily' amino acid side chains (Alanine, Valine, Isoleucine, 
   and Leucine). In the 'Oil Drop' model, these cluster to form a hydrophobic 
   core, shielding themselves from water to stabilize the antibody.

2. ALIPHATIC INDEX (AI) & Tm:
   AI measures the relative volume occupied by aliphatic side chains. A higher 
   AI indicates a denser 'anchor' in the protein core, directly increasing 
   the Predicted Tm (melting temperature).

3. INSTABILITY INDEX (Guruprasad et al., 1990):
   This methodology predicts stability based on 'Dipeptide' logic.
   - DIPEPTIDE LOGIC: The algorithm looks at pairs of neighbors.
   - SCORE < 40: Classified as 'Stable'.
   - SCORE > 40: Classified as 'Unstable' (shorter half-life/prone to degradation).

4. RISK SCORE GRADING:
   - A: High Tm (>75C) & Guruprasad Stable. (Best for manufacturing)
   - B: Good Tm but slightly higher Instability Index.
   - C: Low Tm or High Instability.
   - F: Low Tm AND Very High Instability (>50). (High risk of aggregation/failure)
"""

# Configuration
FOLDERS_TO_PROCESS = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

def assign_risk_score(tm, instability):
    """Assigns a letter grade based on thermal and chemical stability."""
    if tm > 75 and instability < 40:
        return "A (Excellent)"
    elif tm > 70 and instability < 45:
        return "B (Good)"
    elif tm < 65 and instability > 50:
        return "F (High Risk)"
    else:
        return "C (Moderate)"

def calculate_stability_metrics(sequence):
    analyser = ProteinAnalysis(sequence)
    count = analyser.count_amino_acids()
    total = len(sequence)
    
    # 1. Aliphatic Index
    a, b = 2.9, 3.9 
    aliphatic_index = (100 * (count['A'] + a*count['V'] + b*(count['I'] + count['L']))) / total
    
    # 2. Oily Residue Percentage
    oily_percent = ((count['A'] + count['V'] + count['I'] + count['L']) / total) * 100
    
    # 3. Cysteine Content
    cys_percent = (count['C'] / total) * 100
    
    # 4. Guruprasad Instability Index
    instability_index = analyser.instability_index()
    
    # 5. Predicted Tm
    predicted_tm = 55.0 + (0.5 * aliphatic_index) + (2.0 * cys_percent) - (instability_index * 0.1)
    
    # 6. Status and Risk
    status = "Stable" if instability_index < 40 else "Unstable"
    risk_grade = assign_risk_score(predicted_tm, instability_index)

    return {
        "Aliphatic_Index": round(aliphatic_index, 2),
        "Oily_Residue_Percent": round(oily_percent, 2),
        "Cys_Percent": round(cys_percent, 2),
        "Instability_Index": round(instability_index, 2),
        "Predicted_Tm_C": round(predicted_tm, 1),
        "Stability_Status": status,
        "Risk_Score": risk_grade
    }

def run_stability_predictor():
    results = []
    base_path = Path(".")

    for folder in FOLDERS_TO_PROCESS:
        folder_path = base_path / folder
        if not folder_path.exists(): continue
        
        print(f"Scanning: {folder}...")
        for file_path in glob.glob(str(folder_path / "*.py")):
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
                        row_data = {
                            "Antibody": ab_name,
                            "Chain": record.id,
                            "Folder": folder,
                            "Seq_Length": len(record.seq)
                        }
                        metrics = calculate_stability_metrics(str(record.seq))
                        row_data.update(metrics)
                        results.append(row_data)
            except Exception as e:
                print(f"Error in {ab_name}: {e}")

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values(by="Predicted_Tm_C", ascending=False)
        output_file = "Thermal_Stability_Report.xlsx"
        writer = pd.ExcelWriter(output_file, engine='xlsxwriter')
        df.to_excel(writer, sheet_name='Stability Analysis', index=False)

        workbook = writer.book
        worksheet = writer.sheets['Stability Analysis']
        worksheet.autofilter(0, 0, len(df), len(df.columns) - 1)
        worksheet.freeze_panes(1, 0)
        
        # Adding some conditional formatting to make Grade F red
        red_format = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
        worksheet.conditional_format(1, 10, len(df), 10, {
            'type': 'text', 'criteria': 'containing', 'value': 'F', 'format': red_format
        })

        for i, col in enumerate(df.columns):
            max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            worksheet.set_column(i, i, max_len)
        writer.close()

        print(f"\n" + "="*50)
        print(f"GURUPRASAD STABILITY SUMMARY")
        print(f"="*50)
        print(f"Total Chains Analyzed: {len(df)}")
        print(f"Average Predicted Tm:  {df['Predicted_Tm_C'].mean():.2f}°C")
        print(f"Grade A Candidates:    {len(df[df['Risk_Score'].str.startswith('A')])}")
        print(f"Grade F (High Risk):   {len(df[df['Risk_Score'].str.startswith('F')])}")
        print(f"Report Saved:          {output_file}")
        print(f"="*50 + "\n")

if __name__ == "__main__":
    run_stability_predictor()