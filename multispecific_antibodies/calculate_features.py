"""
Calculate sequence features for antibody files using BioPython.
... (omitted docstrings for brevity) ...
"""

from pathlib import Path
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import string 

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
    # ... (omitted parse_py_file function for brevity) ...
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
    ... (omitted docstrings for brevity) ...
    """
    # This function expects 'sequence' to be a single string, which the fix below ensures
    clean_seq = ''.join(aa for aa in sequence.upper() if aa in "ACDEFGHIKLMNPQRSTVWY")
    if not clean_seq:
        return None
    
    analysis = ProteinAnalysis(clean_seq)

    features = {
        'calculated_pI': round(analysis.isoelectric_point(), 2),
        'gravy': round(analysis.gravy(), 2),
        'instability_index': round(analysis.instability_index(), 2),
        'aromaticity': round(analysis.aromaticity(), 4)
    }

    aa_fractions = analysis.amino_acids_percent
    
    for aa in STANDARD_AAS:
         features[f'AA_{aa}_fraction'] = round(aa_fractions.get(aa, 0.0), 5) 

    return features


def get_all_antibody_files(folder: Path) -> list:
    # ... (omitted get_all_antibody_files function for brevity) ...
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
        if folder == 'Whole_mAb':
            if len(heavy_seqs) == 0 or len(light_seqs) == 0:
                print(f"ERROR: {antibody_name} missing chains - heavy: {len(heavy_seqs)}, light: {len(light_seqs)}")
                continue
            # FIX: Ensure we join strings correctly. Access the first string with [0]
            # and repeat the string content, then join them into one sequence string.
            full_sequence = "".join(heavy_seqs[0] * 2 + light_seqs[0] * 2)
            num_heavy = 2
            num_light = 2
        else:
            # This already works because ''.join() converts the list of sequences to one string
            full_sequence = ''.join(heavy_seqs) + ''.join(light_seqs)
            num_heavy = len(heavy_seqs)
            num_light = len(light_seqs)
        
        total_length = len(full_sequence)
        
        # Calculate all features in one go
        # This now receives a single string, fixing the AttributeError
        features = analyze_sequence_features(full_sequence)
        
        if features:
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
    
    core_fieldnames = ['antibody', 'folder', 'heavy_chains', 'light_chains', 'total_length_aa']
    feature_fieldnames = ['calculated_pI', 'gravy', 'instability_index', 'aromaticity'] + [f'AA_{aa}_fraction' for aa in STANDARD_AAS]
    all_fieldnames = core_fieldnames + feature_fieldnames
                  
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=all_fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    print(f"Results saved to: {output_file.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()