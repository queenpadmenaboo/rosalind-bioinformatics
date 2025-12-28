"""
ANTIBODY DEVELOPABILITY & SEQUENCE HOTSPOT AUDIT
================================================
This script identifies Post-Translational Modification (PTM) hotspots and 
sequence motifs that impact Expression, Manufacturing, and Stability.

THEORY & BACKGROUND:
--------------------
1. N-Glycosylation (N-X-S/T): 
   Glycans are large sugar groups. If an "N-X-S/T" motif appears in the CDR (Binding region), 
   the sugar can physically block the antibody from hitting its target. 

2. Deamidation (NG, NS): 
   Asparagine (N) can chemically convert into Aspartic Acid (D). This adds a 
   NEGATIVE charge to the protein, causing "Charge Variants."

3. Isomerization (DG, DS, DT):
   Aspartic Acid (D) flips its chemical structure, often leading to potency loss.

4. Hydrophobic Clusters (Difficult Sequences):
   Clusters of V, I, L, F, W are "oily" patches that cause aggregation or 
   low expression yields in CHO cells.

5. Cleavage Risk (K, R):
   Lysine and Arginine are the primary "cutting sites" for proteases.
"""

from pathlib import Path
import pandas as pd
import re
import os

# --- Path Configuration ---
CANDIDATE_ROOTS = [
    Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"),
    Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies"),
]
ROOT_DIR = next((p.resolve() for p in CANDIDATE_ROOTS if p.exists()), None)
CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

# --- Regex Motifs for Hotspot Detection ---
HOTSPOT_MOTIFS = {
    'N_Glycosylation': r'N[^P][ST]',
    'Deamidation_NG': r'NG',
    'Deamidation_NS': r'NS',
    'Isomerization_D_motifs': r'D[GST]',
    'Hydrophobic_Clusters': r'[VILFW]{5,}'
}

def parse_py_file(filepath: Path) -> dict:
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    sequences = {}
    current_header, current_seq = None, []
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header: sequences[current_header] = ''.join(current_seq)
            current_header, current_seq = line[1:], []
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            if line.lower() not in ('na', 'n/a'): current_seq.append(line)
    if current_header: sequences[current_header] = ''.join(current_seq)
    return sequences

def scan_hotspots(sequence: str) -> dict:
    report = {}
    total_hotspots = 0
    for name, motif in HOTSPOT_MOTIFS.items():
        count = len(re.findall(motif, sequence.upper()))
        report[name] = count
        if name in ['N_Glycosylation', 'Hydrophobic_Clusters']:
            total_hotspots += (count * 2)
        else:
            total_hotspots += count
    
    k_r_count = sequence.upper().count('K') + sequence.upper().count('R')
    report['Cleavage_Risk_K_R'] = k_r_count
    report['Total_Hotspot_Score'] = total_hotspots + (1 if k_r_count > 40 else 0)
    return report

def main():
    print("=" * 60)
    print("RUNNING ANTIBODY HOTSPOT AUDIT (UPDATED COLUMN ORDER)")
    print("=" * 60)
    
    all_data = []
    for folder in CATEGORY_FOLDERS:
        folder_path = ROOT_DIR / folder
        if not folder_path.exists(): continue
        
        for py_file in folder_path.glob('*.py'):
            if py_file.name.startswith(('calculate', 'diagnostic', 'analyze')): continue
            
            try:
                seq_dict = parse_py_file(py_file)
                h_seqs = [s for h, s in seq_dict.items() if any(k in h.lower() for k in ['heavy', 'vh', '_h', 'chain a'])]
                l_seqs = [s for h, s in seq_dict.items() if any(k in h.lower() for k in ['light', 'vl', '_l', 'chain b'])]

                if not h_seqs and not l_seqs:
                    vals = list(seq_dict.values())
                    if len(vals) >= 2: h_seqs, l_seqs = [vals[0]], [vals[1]]
                    elif len(vals) == 1: h_seqs = [vals[0]]
                
                for chain_type, chain_list in [('Heavy', h_seqs), ('Light', l_seqs)]:
                    for i, seq in enumerate(chain_list):
                        res = scan_hotspots(seq)
                        # SWAPPED ORDER: Chain now comes before Folder
                        entry = {
                            'Antibody': py_file.stem,
                            'Chain': f"{chain_type}_{i+1}",
                            'Folder': folder,
                            'AA_Length': len(seq)
                        }
                        entry.update(res)
                        all_data.append(entry)
            except Exception: continue

    if all_data:
        df = pd.DataFrame(all_data)
        df = df.sort_values(by='Total_Hotspot_Score', ascending=False)
        
        output = ROOT_DIR / 'developability_hotspots.xlsx'
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
        df.to_excel(writer, index=False, sheet_name='Hotspot Audit')
        
        workbook = writer.book
        worksheet = writer.sheets['Hotspot Audit']
        danger_fmt = workbook.add_format({'bg_color': '#FFC7CE', 'font_color': '#9C0006'})
        warning_fmt = workbook.add_format({'bg_color': '#FFEB9C', 'font_color': '#9C5700'})
        
        worksheet.freeze_panes(1, 0)
        worksheet.autofilter(0, 0, len(df), len(df.columns) - 1)
        
        for i, col in enumerate(df.columns):
            worksheet.set_column(i, i, max(len(col), 12) + 2)

        score_idx = df.columns.get_loc('Total_Hotspot_Score')
        worksheet.conditional_format(1, score_idx, len(df), score_idx, 
                                     {'type': 'cell', 'criteria': '>=', 'value': 5, 'format': danger_fmt})
        worksheet.conditional_format(1, score_idx, len(df), score_idx, 
                                     {'type': 'cell', 'criteria': 'between', 'minimum': 3, 'maximum': 4, 'format': warning_fmt})
        
        writer.close()
        print(f"DONE! Report: {output.name}")
    print("=" * 60)

if __name__ == "__main__":
    main()