import os
from pathlib import Path
import pandas as pd

"""
ANTIBODY DIAGNOSTIC TOOL
========================
Use this script to find out WHY some antibodies are missing from your main report.
It scans the folders, reads the files, and checks for:
1. Missing '>' headers.
2. Missing 'heavy' or 'light' keywords in headers.
3. Empty sequences.
"""

# --- Path Configuration ---
CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), None)

CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

def parse_py_file_diagnostic(filepath: Path):
    """Deep scan of the file to see what headers exist."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        return ["FILE_READ_ERROR"], {}

    headers_found = []
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header:
                sequences[current_header] = ''.join(current_seq)
            current_header = line[1:]
            headers_found.append(current_header)
            current_seq = []
        # Identifying potential sequence lines
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            current_seq.append(line)
            
    if current_header:
        sequences[current_header] = ''.join(current_seq)
        
    return headers_found, sequences

def run_diagnostic():
    if not ROOT_DIR:
        print("!!! ERROR: ROOT_DIR not found. Check your CANDIDATE_ROOTS paths.")
        return

    audit_results = []
    print(f"--- STARTING AUDIT AT: {ROOT_DIR} ---")

    for folder in CATEGORY_FOLDERS:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists():
            print(f"Folder skipped (not found): {folder}")
            continue

        files = list(folder_path.glob("*.py"))
        for f in files:
            # Skip utility scripts
            if f.name in ['calculate_features.py', 'antibody_diagnostic_tool.py']:
                continue

            headers, seqs = parse_py_file_diagnostic(f)
            
            # Check for standard naming keywords
            heavy_match = [h for h in headers if any(k in h.lower() for k in ['heavy', 'vh', 'vhh'])]
            light_match = [h for h in headers if any(k in h.lower() for k in ['light', 'vl', 'vkappa', 'vlambda'])]
            
            status = "PASS"
            reason = ""
            
            # THE LOGIC THAT DROPS FILES:
            if not headers:
                status = "FAIL"
                reason = "No FASTA headers ('>') detected"
            elif folder == 'Whole_mAb' and (not heavy_match or not light_match):
                status = "FAIL"
                reason = f"Whole_mAb folder requires both H and L (Found H:{len(heavy_match)} L:{len(light_match)})"
            elif not heavy_match and not light_match:
                status = "FAIL"
                reason = "No chain labels found (Header needs 'heavy', 'light', 'VH', or 'VL')"
            
            audit_results.append({
                "Antibody": f.stem,
                "Folder": folder,
                "Status": status,
                "Reason": reason,
                "Raw_Headers": "|".join(headers)
            })

    df = pd.DataFrame(audit_results)
    failures = df[df['Status'] == "FAIL"]
    
    print("\n" + "="*50)
    print(f"AUDIT SUMMARY")
    print("="*50)
    print(f"Total Files Scanned: {len(df)}")
    print(f"Successful:          {len(df) - len(failures)}")
    print(f"Failed/Missing:      {len(failures)}")
    print("="*50)
    
    if not failures.empty:
        print("\nCOMMON FAILURE REASONS:")
        print(failures['Reason'].value_counts())
        
        log_file = ROOT_DIR / "MISSING_ANTIBODIES_LOG.xlsx"
        failures.to_excel(log_file, index=False)
        print(f"\n[!] Detailed log of failed files saved to: {log_file.name}")
    else:
        print("\nAll files are structurally correct for the feature calculator!")

if __name__ == "__main__":
    run_diagnostic()