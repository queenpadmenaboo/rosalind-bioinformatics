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
    - Calculates additional features: GRAVY index, Instability Index, Aromaticity, and Amino Acid Composition (%)

NOTE ON pI:
    - This script calculates THEORETICAL pI based on amino acid sequence
    - Uses BioPython's ProteinAnalysis.isoelectric_point() method (uses 'Bjellqvist' pKa set by default)
    - Experimental pI (from icIEF, capillary isoelectric focusing) may differ
    - Differences due to post-translational modifications (glycosylation, etc.)

OUTPUT:
    - Console report with sample results
    - sequence_features.csv (created in multispecific_antibodies folder)
    - Columns: antibody, folder, heavy_chains, light_chains, total_length_aa, calculated_pI, gravy, instability_index, aromaticity, AA_A_percent, AA_C_percent, etc.

USAGE:
    python calculate_features.py

REQUIRES:
    pip install biopython
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import string # To get all standard amino acid letters

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

# Define the 20 standard amino acids for CSV header consistency
STANDARD_AAS = sorted("ACDEFGHIKLMNPQRSTVWY")

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

def analyze_sequence_features(sequence: str) -> dict:
    """
    Calculates various theoretical features using Biopython's ProteinAnalysis.
    """
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not clean_seq:
        return None
    
    analysis = ProteinAnalysis(clean_seq)

    # Calculate basic features
    features = {
        'calculated_pI': round(analysis.isoelectric_point(), 2),
        'gravy': round(analysis.gravy(), 2),
        'instability_index': round(analysis.instability_index(), 2),
        'aromaticity': round(analysis.aromaticity(), 4)
    }

    # Calculate Amino Acid Percentages (Composition) using the NEW attribute
    # Changed from: aa_percents = analysis.get_amino_acids_percent()
    aa_percents = analysis.amino_acids_percent 
    
    for aa in STANDARD_AAS:
         # Ensure all 20 AA columns exist, even if percentage is 0.0
         features[f'AA_{aa}_percent'] = round(aa_percents.get(aa, 0.0) * 100, 3) 

    return features


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
        
        # Build full molecule sequence
        if folder == 'Whole_mAb':
            if len(heavy_seqs) == 0 or len(light_seqs) == 0:
                print(f"ERROR: {antibody_name} missing chains - heavy: {len(heavy_seqs)}, light: {len(light_seqs)}")
                continue
            full_sequence = heavy_seqs[0] * 2 + light_seqs[0] * 2
            num_heavy = 2
            num_light = 2
        else:
            full_sequence = ''.join(heavy_seqs) + ''.join(light_seqs)
            num_heavy = len(heavy_seqs)
            num_light = len(light_seqs)
        
        total_length = len(full_sequence)
        
        # Calculate all features in one go
        features = analyze_sequence_features(full_sequence)
        
        if features:
            # Combine core data with calculated features
            antibody_data = {
                'antibody': antibody_name,
                'folder': folder,
                'heavy_chains': num_heavy,
                'light_chains': num_light,
                'total_length_aa': total_length,
            }
            antibody_data.update(features)
            results.append(antibody_data)
    
    print("Sample results (first 5):")
    print("-" * 70)
    for r in results[:5]:
        print(f"  {r['antibody']}: pI={r['calculated_pI']}, len={r['total_length_aa']}, gravy={r['gravy']}, instab={r['instability_index']}")
    if len(results) > 5:
        print(f"  ... and {len(results) - 5} more\n")
    
    output_file = ANTIBODY_FOLDER / 'sequence_features.csv'
    
    # Define all fieldnames for the CSV writer dynamically
    core_fieldnames = ['antibody', 'folder', 'heavy_chains', 'light_chains', 'total_length_aa']
    feature_fieldnames = ['calculated_pI', 'gravy', 'instability_index', 'aromaticity'] + [f'AA_{aa}_percent' for aa in STANDARD_AAS]
    all_fieldnames = core_fieldnames + feature_fieldnames
                  
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=all_fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results saved to: {output_file.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()
