import os
# 1. GPU SECURITY BYPASS
os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] = "1"

import re
import torch
import warnings
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# Suppress noise
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

# --- CONFIG & EXCLUSIONS (STRICTLY PRESERVED) ---
BASE_DIR = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")
OUTPUT_DIR = BASE_DIR / "PDB_Output_Files_GPU"

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py', 'analyze_hotspots.py', 'aggregation_predictor.py', 'Aggregation_Risk_Report.xlsx',
    'antibody_diagnostic_tool.py', 'Antibody_Comparison_Report_2025.xlsx', 'build_pnas_shadow.py',
    'compare_hydrophobicity.py', 'developability_hotspots.xlsx', 'mAb_truth_engine.py',
    'mAb_Truth_Engine_Master.xlsx', 'MISSING_ANTIBODIES_LOG.xlsx', 'pnas_validator.py',
    'PNAS_VS_CALCULATIONS.xlsx', '3D_structure_builder_gpu.py', 'format_pairing_predictor.py',
    'immunogenicity_predictor.py', 'thermal_stability_predictor.py'
}

# Amino acids that are "Greasy" (Hydrophobic) for ML Predictor
HYDROPHOBIC_RESIDUES = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']

print("--- INITIALIZING BIOPHYSICAL ANALYSIS ENGINE ---")

def extract_chains_raw_text(file_path):
    """Original Regex Logic preserved to keep sequence data linked"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        blocks = re.findall(r'>(.*?)(?=>|\"\"\"|\'\'\'|\Z)', content, re.DOTALL)
        
        # FIX: Changed to lists to store ALL chains from multispecific antibodies
        h_chains = []
        l_chains = []
        
        for block in blocks:
            lines = block.strip().split('\n')
            if not lines: continue
            header = lines[0].lower()
            sequence = re.sub(r'[^A-Z]', '', "".join(lines[1:]).upper())
            if not sequence: continue
            
            # FIX: Append to list instead of overwriting
            if any(h in header for h in ["heavy", "vh", "vhh"]):
                h_chains.append(sequence)
            elif any(l in header for l in ["light", "vl"]):
                l_chains.append(sequence)
        
        # Return joined chains so Excel shows everything
        return {
            "H": " | ".join(h_chains) if h_chains else "N/A",
            "L": " | ".join(l_chains) if l_chains else "N/A"
        }
    except Exception:
        return None

def run_pipeline():
    if not OUTPUT_DIR.exists():
        print(f"ERROR: Output directory {OUTPUT_DIR} not found!")
        return

    # 1. FIND THE PY FILES
    antibody_files = [
        f for f in BASE_DIR.rglob("*.py") 
        if f.name not in EXCLUDE_FILES 
        and "shadow_benchmarks" not in f.parts 
        and not f.name.startswith("._")
    ]
    
    # 2. INDEX ALL PDBs IN FORMAT SUBFOLDERS
    print(f"Indexing structures in {OUTPUT_DIR.name}...")
    pdb_map = {f.stem.lower(): f for f in OUTPUT_DIR.rglob("*.pdb")}
    
    print(f"Found {len(antibody_files)} antibody scripts.")
    print(f"Found {len(pdb_map)} PDB files across all subfolders.")
    
    master_data = []
    parser = PDBParser(QUIET=True)
    sr = ShrakeRupley()

    # 3. ANALYZE
    for f_path in tqdm(antibody_files, desc="Calculating Biophysics"):
        ab_name = f_path.stem
        ab_name_lower = ab_name.lower()
        
        if ab_name_lower in pdb_map:
            pdb_path = pdb_map[ab_name_lower]
            
            try:
                struct = parser.get_structure(ab_name, str(pdb_path))
                
                # Global SASA
                sr.compute(struct, level="S")
                total_sasa = struct.sasa
                
                # Hydrophobic SASA
                sr.compute(struct, level="R")
                hydro_sasa = 0
                for res in struct.get_residues():
                    if res.get_resname() in HYDROPHOBIC_RESIDUES:
                        hydro_sasa += getattr(res, 'sasa', 0)

                seq_data = extract_chains_raw_text(f_path)
                
                master_data.append({
                    "Antibody_Name": ab_name,
                    "Format_Type": pdb_path.parent.name,
                    "Total_SASA": round(total_sasa, 2),
                    "Hydrophobic_SASA": round(hydro_sasa, 2),
                    "Hydro_Ratio": round(hydro_sasa / total_sasa, 3) if total_sasa > 0 else 0,
                    "Sequence_H": seq_data["H"] if seq_data else "N/A",
                    "Sequence_L": seq_data["L"] if seq_data else "N/A" # Added L column for completeness
                })

            except Exception as e:
                print(f"\n[SKIP] Error processing {ab_name}: {e}")

    # 4. SAVE FINAL DATASET (Excel with Formatting)
    if master_data:
        df = pd.DataFrame(master_data)
        excel_out = BASE_DIR / "all_antibody_sasa_chains.xlsx"
        
        # Use XlsxWriter engine for advanced formatting
        writer = pd.ExcelWriter(excel_out, engine='xlsxwriter')
        df.to_excel(writer, index=False, sheet_name='SASA_Data')
        
        workbook  = writer.book
        worksheet = writer.sheets['SASA_Data']

        # A. Freeze the first row
        worksheet.freeze_panes(1, 0)

        # B. Add filters to the first row
        worksheet.autofilter(0, 0, 0, len(df.columns) - 1)

        # C. Auto-adjust column width
        for i, col in enumerate(df.columns):
            column_len = max(df[col].astype(str).map(len).max(), len(col)) + 2
            # Cap width for sequence column to keep it readable
            if "Sequence" in col: column_len = 50 
            worksheet.set_column(i, i, column_len)

        writer.close()
        print(f"\n=== SUCCESS: Saved {len(df)} antibodies to {excel_out.name} ===")

        # 5. AUTOMATED CLEANUP
        print("\n--- CLEANING UP REDUNDANT .SASA.TXT FILES ---")
        deleted_count = 0
        for old_file in OUTPUT_DIR.rglob("*.sasa.txt"):
            try:
                old_file.unlink()
                deleted_count += 1
            except Exception as e:
                print(f"Failed to delete {old_file.name}: {e}")
        print(f"Cleanup complete. Removed {deleted_count} files.")

    else:
        print("\n[!] No data collected. Check if PDB filenames match .py filenames.")

if __name__ == "__main__":
    run_pipeline()