"""
Calculate sequence features for antibody files using BioPython.

FOLDER STRUCTURE:
    multispecific_antibodies/
    ├── Bispecific_mAb/
    ├── Bispecific_scFv/
    ├── Other_Formats/
    ├── Whole_mAb/
    ├── calculate_features.py
    ├── categorize_antibody_format.py
    ├── readme_count.py
    ├── README.md
    ├── sabdabconverter.py
    ├── selenium_antibody_scraper.py
    ├── thera_sabdab_scraper.py
    ├── validate_antibody_sequences.py
    └── validation_report.csv

PURPOSE:
    - Scans all 4 subfolders for antibody .py files
    - Parses FASTA sequences from each file
    - Calculates isoelectric point (pI) for each antibody molecule

NOTE ON pI:
    - This script calculates THEORETICAL pI based on amino acid sequence
    - Uses BioPython's ProteinAnalysis.isoelectric_point() method (uses 'Bjellqvist' pKa set by default)
    - Experimental pI (from icIEF, capillary isoelectric focusing) may differ
    - Differences due to post-translational modifications (glycosylation, etc.)

OUTPUT:
    - Console report with sample results
    - sequence_features.csv (created in multispecific_antibodies folder)
    - Columns: antibody, folder, heavy_chains, light_chains, total_length_aa, calculated_pI

USAGE:
    python calculate_features.py

REQUIRES:
    pip install biopython
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# Removed the manual IP and ProtParamData imports to avoid errors
import csv

# ============================
# CONFIGURATION
# ============================

ANTIBODY_FOLDER = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")

CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py',
    'categorize_antibody_format.py', 'categorize_simple.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py'
}

# ============================
# FUNCTIONS
# ============================

def parse_py_file(filepath: Path) -> dict:
    """
    Parses FASTA format sequences embedded within a .py file.
    Skips chains marked as "na" or "n/a" (not applicable)
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in content.split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_header:
                seq = ''.join(current_seq)
                if seq.upper() not in ('NA', 'N/A', ''):
                    sequences[current_header] = seq
            current_header = line[1:]
            current_seq = []
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            if line.lower() not in ('na', 'n/a'):
                current_seq.append(line)
    
    if current_header:
        seq = ''.join(current_seq)
        if seq.upper() not in ('NA', 'N/A', ''):
            sequences[current_header] = seq
    
    return sequences


def calculate_isoelectric_point(sequence: str) -> float:
    """
    Calculates the theoretical pI using Biopython's default method (Bjellqvist scale).
    """
    # Clean the sequence to only include standard amino acids
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in 'ACDEFGHIKLMNPQRSTVWY')
    if not clean_seq:
        return None
    
    # Use the reliable ProteinAnalysis method, which uses the default scale internally
    analysis = ProteinAnalysis(clean_seq)
    return round(analysis.isoelectric_point(), 2)


def get_all_antibody_files(folder: Path) -> list:
    """Finds all antibody .py files in category subfolders."""
    py_files = []
    for subfolder in CATEGORY_FOLDERS:
        subfolder_path = folder / subfolder
        if subfolder_path.exists():
            for f in subfolder_path.glob('*.py'):
                if f.is_file() and f.name not in EXCLUDE_FILES:
                    py_files.append(f)
    return py_files


def main():
    print("=" * 70)
    print("ANTIBODY FEATURE CALCULATION")
    print("=" * 70)
    
    py_files = get_all_antibody_files(ANTIBODY_FOLDER)
    print(f"Found {len(py_files)} antibody files\n")
    
    results = []
    
    for py_file in py_files:
        antibody_name = py_file.stem
        folder = py_file.parent.name
        
        try:
            sequences = parse_py_file(py_file)
        except Exception as e:
            print(f"ERROR parsing {antibody_name}: {e}")
            continue
        
        heavy_seqs = []
        light_seqs = []
        
        for header, seq in sequences.items():
            header_lower = header.lower()
            if 'heavy' in header_lower:
                heavy_seqs.append(seq)
            elif 'light' in header_lower:
                light_seqs.append(seq)
        
        # Build full molecule sequence (concatenates all chains present)
        # Whole_mAb format assumes 2 identical heavy + 2 identical light chains
        if folder == 'Whole_mAb':
            if len(heavy_seqs) == 0 or len(light_seqs) == 0:
                print(f"ERROR: {antibody_name} missing chains - heavy: {len(heavy_seqs)}, light: {len(light_seqs)}")
                continue
            # Duplicate the single chain sequence to represent the full molecule
            # Use [0] to access the string content of the first element in the list
            full_sequence = heavy_seqs[0] * 2 + light_seqs[0] * 2
            num_heavy = 2
            num_light = 2
        else:
            # For bispecifics/scFvs/etc, combine whatever unique chains are provided
            full_sequence = ''.join(heavy_seqs) + ''.join(light_seqs)
            num_heavy = len(heavy_seqs)
            num_light = len(light_seqs)
        
        total_length = len(full_sequence)
        # Calculate the single pI value for the entire concatenated sequence
        pI = calculate_isoelectric_point(full_sequence)
        
        results.append({
            'antibody': antibody_name,
            'folder': folder,
            'heavy_chains': num_heavy,
            'light_chains': num_light,
            'total_length_aa': total_length,
            'calculated_pI': pI
        })
    
    print("Sample results:")
    print("-" * 70)
    for r in results[:10]:
        print(f"  {r['antibody']}: pI = {r['calculated_pI']}, length = {r['total_length_aa']} aa")
    if len(results) > 10:
        print(f"  ... and {len(results) - 10} more\n")
    
    output_file = ANTIBODY_FOLDER / 'sequence_features.csv'
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['antibody', 'folder', 'heavy_chains', 'light_chains', 'total_length_aa', 'calculated_pI'])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results saved to: {output_file.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()