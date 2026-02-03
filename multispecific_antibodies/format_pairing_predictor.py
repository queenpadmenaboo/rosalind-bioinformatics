import pandas as pd
import re
from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ======================================================================================
# MODULE: FORMAT-SPECIFIC PAIRING PREDICTOR (ASYMMETRY & COMPLEXITY)
# ======================================================================================
# THEORY & CALCULATIONS:
#
# 1. ARM SIZE ASYMMETRY (MW Delta):
#    Calculates mass difference between arms. High delta enables SEC purification.
#
# 2. CHARGE/HYDROPHOBICITY DIFFERENCES:
#    Uses pH 5.5 titration and GRAVY to determine IEX/HIC separation efficiency.
#
# 3. PAIRING COMPLEXITY & YIELD PREDICTION:
#    Estimates % Yield based on physical divergence. Low delta = high homodimer junk.
#
# 4. LINKER CONSTRAINTS (Short vs Long):
#    Linkers < 12aa force Diabody formation (inter-chain pairing).
# ======================================================================================

CANDIDATE_ROOTS = [
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p for p in CANDIDATE_ROOTS if p.exists()), Path("."))
FORMULATION_PH = 5.5

def detect_linker_constraints(sequence):
    gs_linkers = re.findall(r'((?:G{3,5}S){1,})', sequence)
    if gs_linkers:
        longest_linker = max(gs_linkers, key=len)
        l_len = len(longest_linker)
        l_count = len(gs_linkers)
        if l_len < 12: 
            return f"SHORT ({l_len}aa) - Diabody Risk", l_len, l_count
        return f"LONG ({l_len}aa) - Flexible", l_len, l_count
    return "None", 0, 0

def get_physics(sequence):
    clean = re.sub(r'[^A-Z]', '', sequence.upper())
    if not clean: return None
    analysis = ProteinAnalysis(clean)
    l_desc, l_len, l_count = detect_linker_constraints(clean)
    return {
        "mw": analysis.molecular_weight(),
        "pI": analysis.isoelectric_point(),
        "charge": analysis.charge_at_pH(FORMULATION_PH),
        "hydro": analysis.gravy(),
        "linker_desc": l_desc,
        "linker_len": l_len,
        "linker_count": l_count
    }

def run_predictor():
    print(f"--- EXECUTING FORMAT PAIRING PREDICTION ---")
    results = []
    folders = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

    for folder in folders:
        f_path = ROOT_DIR / folder
        if not f_path.exists(): continue
        for py_file in f_path.glob("*.py"):
            try:
                with open(py_file, 'r', encoding='utf-8') as f:
                    content = f.read()
                blocks = re.findall(r'["\']{3}(.*?)["\']{3}', content, re.DOTALL)
                phys = []
                for b in blocks:
                    for h, s in re.findall(r'>([^\n]+)\n([A-Z\n\s]+)', b):
                        p = get_physics(s)
                        if p: phys.append(p)
                
                if not phys: continue
                
                asym, c_diff, h_diff, complex_val, yield_val = 0, 0, 0, "Single-Chain", "N/A"
                l_logic = next((m['linker_desc'] for m in phys if m['linker_desc'] != "None"), "Standard")
                max_l_len = max([m['linker_len'] for m in phys])
                total_l_count = sum([m['linker_count'] for m in phys])

                if len(phys) >= 2:
                    asym = abs(phys[0]['mw'] - phys[1]['mw'])
                    c_diff = abs(phys[0]['charge'] - phys[1]['charge'])
                    h_diff = abs(phys[0]['hydro'] - phys[1]['hydro'])
                    
                    if folder == "Bispecific_mAb":
                        complex_val = "HIGH (4-Chain Problem)"
                    elif folder == "Bispecific_scFv":
                        complex_val = "LOW (Tethered)"
                    else:
                        complex_val = "MEDIUM"
                    
                    if asym < 2000 and c_diff < 1.0: 
                        yield_val = "LOW (<30%) - Mispairing Dominant"
                    elif asym > 5000 or c_diff > 3.0: 
                        yield_val = "HIGH (>85%) - Clean Assembly"
                    else: 
                        yield_val = "MODERATE (~60%)"

                results.append({
                    "Antibody": py_file.stem,
                    "Format": folder,
                    "Arm_Size_Asymmetry_Da": round(asym, 2),
                    "Charge_Difference": round(c_diff, 2),
                    "Hydro_Difference": round(h_diff, 3),
                    "Pairing_Complexity": complex_val,
                    "Linker_Constraints": l_logic,
                    "Max_Linker_Length": max_l_len,
                    "Total_Linker_Count": total_l_count,
                    "Predicted_Heterodimer_Yield": yield_val
                })
            except Exception as e:
                print(f"Error processing {py_file.name}: {e}")

    if results:
        df = pd.DataFrame(results).sort_values(by="Arm_Size_Asymmetry_Da", ascending=False)
        out = ROOT_DIR / "Format_Pairing_Prediction.xlsx"
        
        with pd.ExcelWriter(out, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Pairing_Predictions')
            wb, ws = writer.book, writer.sheets['Pairing_Predictions']
            
            # Format: Red for Risk, Green for High Yield
            risk_fmt = wb.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
            success_fmt = wb.add_format({'bg_color': '#C6EFCE', 'font_color': '#006100'})
            
            # 1. ADD AUTO-FILTER TO FIRST ROW
            ws.autofilter(0, 0, len(df), len(df.columns) - 1)
            
            # 2. FREEZE FIRST ROW
            ws.freeze_panes(1, 0)
            
            for i, col in enumerate(df.columns):
                ws.set_column(i, i, max(df[col].astype(str).map(len).max(), len(col)) + 4)
            
            ws.conditional_format(1, 5, len(df), 5, {'type': 'text', 'criteria': 'containing', 'value': 'HIGH', 'format': risk_fmt})
            ws.conditional_format(1, 6, len(df), 6, {'type': 'text', 'criteria': 'containing', 'value': 'SHORT', 'format': risk_fmt})
            ws.conditional_format(1, 9, len(df), 9, {'type': 'text', 'criteria': 'containing', 'value': 'LOW', 'format': risk_fmt})
            ws.conditional_format(1, 9, len(df), 9, {'type': 'text', 'criteria': 'containing', 'value': 'HIGH', 'format': success_fmt})

        print(f"SUCCESS: Pairing Predictor saved to {out.name}")

if __name__ == "__main__":
    run_predictor()