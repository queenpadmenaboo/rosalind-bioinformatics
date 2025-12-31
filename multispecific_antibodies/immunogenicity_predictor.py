import os
import re
import pandas as pd
from pathlib import Path

# ======================================================================================
# MODULE: IMMUNOGENICITY & FRAMEWORK (FR) RISK PREDICTOR
# ======================================================================================
# LITERATURE SOURCES & SCIENTIFIC RATIONALE:
#
# 1. T-CELL EPITOPES (MHC-II BINDING ANCHORS):
#    Source: Marshall et al. (2003) "MHC class II-binding conformations."
#    Logic: The 'Binding Groove' of MHC-II (HLA-DR) has a strict requirement for 
#    bulky hydrophobic 'anchor' residues at positions P1, P4, P6, and P9. 
#    Clusters of [F, W, Y, V, I, L] represent high-affinity binding potential.
#
# 2. FRAMEWORK (FR) RISK & HUMANIZATION:
#    Source: Harding et al. (2010) "The immunogenicity of therapeutic proteins."
#    Logic: Clinical data shows that "Non-Germline" mutations in the FR regions 
#    are 3x more likely to trigger Anti-Drug Antibodies (ADAs) than CDR mutations. 
#    Rare clusters of [C, P, G] disrupt the standard IgG scaffold, creating 
#    structural 'kinks' that the immune system flags as foreign.
#
# 3. CHEMICAL NEO-EPITOPES (DEAMIDATION):
#    Source: Vlasak et al. (2008) "Effects of deamidation on antibody structure."
#    Logic: The spontaneous conversion of Asparagine (N) to Aspartate (D) at 'NG' 
#    sites introduces a negative charge. This 'chemical scar' can create a 
#    neo-epitope, causing an immune response to an otherwise human protein.
# ======================================================================================

# --- Path Configuration ---
CANDIDATE_ROOTS = [
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p for p in CANDIDATE_ROOTS if p.exists()), Path("."))

CSV_SEARCH_PATHS = [Path(r"C:\Users\meeko"), Path(r"C:\Users\bunsr"), ROOT_DIR]

def get_master_csv():
    for search_path in CSV_SEARCH_PATHS:
        if search_path.exists():
            csv_files = list(search_path.glob("*TheraSAbDab*.csv"))
            if csv_files:
                return sorted(csv_files)[-1]
    return None

MASTER_CSV = get_master_csv()

# --- Define Search Patterns ---
RISK_PATTERNS = {
    'T_Cell_Epitope_Anchors': r'[F-WYVIL]{3,}', 
    'Framework_Structural_Kinks': r'[CPG]{3,}', 
    'PTM_NeoEpitope_Risk': r'NG'                
}

def analyze_immunogenicity(sequence):
    seq = str(sequence).upper().strip()
    results = {}
    total_risk_score = 0
    for label, pattern in RISK_PATTERNS.items():
        matches = len(re.findall(pattern, seq))
        results[label] = matches
        if label == 'T_Cell_Epitope_Anchors':
            total_risk_score += (matches * 4)
        else:
            total_risk_score += (matches * 2)
    results['Immuno_Risk_Density'] = round((total_risk_score / len(seq)) * 100, 2)
    return results

def run_immunogenicity_scan():
    print(f"--- STARTING IMMUNOGENICITY SCAN ---")
    
    target_data = []
    if MASTER_CSV:
        try:
            print(f"Found Master File: {MASTER_CSV}")
            df_master = pd.read_csv(MASTER_CSV, encoding='utf-8-sig')
            df_master.columns = df_master.columns.str.strip() 
            target_data = [
                (str(t).strip().lower(), str(g).strip()) 
                for t, g in zip(df_master['Therapeutic'], df_master['Target'])
                if pd.notnull(t)
            ]
            print(f"Successfully loaded {len(target_data)} therapeutics.")
        except Exception as e:
            print(f"CSV Error: {e}")

    results = []
    folders = ["Bispecific_mAb", "Bispecific_scFv", "Other_Formats", "Whole_mAb"]

    for folder in folders:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists(): continue
        
        for file_path in folder_path.glob("*.py"):
            fname_clean = file_path.stem.strip().lower()
            found_target = "Target Not Found"
            for thera_name, target_val in target_data:
                if thera_name == fname_clean or thera_name in fname_clean or fname_clean in thera_name:
                    found_target = target_val
                    break
            
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                all_blocks = re.findall(r'["\']{3}(.*?)["\']{3}', content, re.DOTALL)
                for block in all_blocks:
                    fasta_matches = re.findall(r'>([^\n]+)\n([A-Z\n\s]+)', block)
                    for header, seq_raw in fasta_matches:
                        clean_seq = re.sub(r'[^A-Z]', '', seq_raw.upper())
                        if len(clean_seq) < 20: continue
                        imm_data = analyze_immunogenicity(clean_seq)
                        results.append({
                            "Antibody": file_path.stem,
                            "Target": found_target,
                            "Chain": header.strip(),
                            "Format": folder,
                            **imm_data
                        })
            except Exception as e:
                print(f"Error in {file_path.name}: {e}")

    if results:
        df = pd.DataFrame(results).sort_values(by="Immuno_Risk_Density", ascending=False)
        output_file = ROOT_DIR / "Immunogenicity_Risk_Report.xlsx"
        
        # --- AUTO-WIDTH EXPORT LOGIC ---
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Risk_Report')
            worksheet = writer.sheets['Risk_Report']
            
            # Iterate through columns and find the max length for each
            for i, col in enumerate(df.columns):
                # Max length of string in column or column name
                max_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
                worksheet.set_column(i, i, max_len)
                
        print(f"SUCCESS: {len(results)} sequences processed into {output_file.name}")
    else:
        print("No sequences found.")

if __name__ == "__main__":
    run_immunogenicity_scan()