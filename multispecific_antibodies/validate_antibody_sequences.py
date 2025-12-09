"""
Validates antibody Python files against TheraSAbDab CSV database.

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
    - Scans root folder + all 4 subfolders for .py files
    - Checks if each antibody exists in CSV
    - Compares each chain sequence in .py file against CSV sequence
    - Validates sequences contain only valid amino acids (ACDEFGHIKLMNPQRSTVWY)
    - Counts chains per antibody
    - Counts files in each folder

OUTPUT:
    - Console report (matched, not in CSV, invalid)
    - validation_report.csv in multispecific_antibodies folder

USAGE:
    python validate_antibody_sequences.py
"""

import csv
from pathlib import Path
from typing import Dict, List, Set
from dataclasses import dataclass

# ============================
# CONFIGURATION
# ============================

CSV_PATH = Path(r"C:\Users\meeko\TheraSAbDab_SeqStruc_08Dec2025.csv")
ANTIBODY_FOLDER = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")

VALID_AA: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")

CATEGORY_FOLDERS = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']

EXCLUDE_FILES = {
    'readme_count.py',
    'sabdabconverter.py',
    'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py',
    'validate_antibody_sequences.py',
    'categorize_antibody_format.py',
    'categorize_simple.py',
    'therasabdab_analyze_formats.py',
    'calculate_features.py',
    'fix_sequences.py',
    'validation_report.csv'
}

# ============================
# DATA STRUCTURES
# ============================

@dataclass
class AntibodyData:
    # Stores antibody info from CSV
    name: str       # Therapeutic name
    heavy1: str     # Heavy chain sequence (primary)
    light1: str     # Light chain sequence (primary)
    heavy2: str     # Heavy chain 2 (bispecifics only)
    light2: str     # Light chain 2 (bispecifics only)
    format: str     # Antibody format (Whole mAb, Bispecific, etc.)
    status: str     # Clinical status (Active, Discontinued, etc.)

@dataclass
class ValidationResult:
    # Stores validation results for reporting
    matched: List[Dict]           # Antibodies that passed validation
    not_in_csv: List[Dict]        # Files not found in CSV
    invalid_sequences: List[Dict] # Files with sequence errors

# ============================
# FUNCTIONS
# ============================

def load_csv_data(csv_file: Path) -> Dict[str, AntibodyData]:
    # Reads TheraSAbDab CSV and returns dict of antibody data
    # Key = antibody name (lowercase, underscores)
    antibodies = {}
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            therapeutic_raw = row['Therapeutic'].strip()
            key = therapeutic_raw.lower().replace(' ', '_').replace('-', '_')
            antibodies[key] = AntibodyData(
                name=therapeutic_raw,
                heavy1=row['HeavySequence'].strip(),
                light1=row['LightSequence'].strip(),
                heavy2=row['HeavySequence(ifbispec)'].strip(),
                light2=row['LightSequence(ifbispec)'].strip(),
                format=row['Format'],
                status=row['Est. Status']
            )
    return antibodies


def parse_py_file(filepath: Path) -> Dict[str, str]:
    # Parses antibody .py file and extracts FASTA sequences
    # Returns dict: {header: sequence}
    # Handles 'na' as missing chains (valid for some formats)
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
                sequences[current_header] = 'NA' if seq.lower() in ('na', 'n/a', '') else seq
            current_header = line[1:]
            current_seq = []
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            if line.lower() in ('na', 'n/a'):
                current_seq.append('NA')
            else:
                current_seq.append(line)
    
    if current_header:
        seq = ''.join(current_seq)
        sequences[current_header] = 'NA' if seq.lower() in ('na', 'n/a', '') else seq
    
    return sequences


def validate_sequence(seq: str) -> bool:
    # Checks if sequence contains only valid amino acids
    # Returns True if valid, False if invalid characters found
    # 'NA' or empty = valid (means missing chain)
    seq_clean = seq.strip().upper()
    if seq_clean in ('NA', 'N/A', ''):
        return True
    return all(aa in VALID_AA for aa in seq_clean)


def count_chains(sequences: Dict[str, str]) -> int:
    # Counts non-NA sequences
    # Helps identify bispecifics (4 chains) vs whole mAb (2 chains)
    return sum(1 for seq in sequences.values() if seq not in ('NA', ''))


def get_all_antibody_files(py_folder: Path) -> List[Path]:
    # Gets all antibody .py files from root and 4 subfolders
    # Excludes utility scripts (validate, categorize, etc.)
    py_files = []
    for f in py_folder.glob('*.py'):
        if f.is_file() and f.name not in EXCLUDE_FILES:
            py_files.append(f)
    for subfolder in CATEGORY_FOLDERS:
        subfolder_path = py_folder / subfolder
        if subfolder_path.exists():
            for f in subfolder_path.glob('*.py'):
                if f.is_file() and f.name not in EXCLUDE_FILES:
                    py_files.append(f)
    return py_files


def validate_antibody_files(py_folder: Path, csv_antibodies: Dict[str, AntibodyData]) -> ValidationResult:
    # Main validation logic
    # Checks: 1) antibody exists in CSV
    #         2) sequences match CSV exactly
    #         3) sequences contain valid amino acids
    matched = []
    not_in_csv = []
    invalid_sequences = []
    
    py_files = get_all_antibody_files(py_folder)
    
    print(f"Found {len(py_files)} antibody files to validate\n")
    
    for py_file in py_files:
        antibody_name = py_file.stem
        
        try:
            py_sequences = parse_py_file(py_file)
        except Exception as e:
            invalid_sequences.append({'name': antibody_name, 'error': f"Parse error: {str(e)}"})
            continue
        
        if not py_sequences:
            invalid_sequences.append({'name': antibody_name, 'error': "No sequences found"})
            continue
        
        if antibody_name not in csv_antibodies:
            not_in_csv.append({'name': antibody_name, 'chain_count': count_chains(py_sequences), 'path': str(py_file)})
            continue
        
        csv_data = csv_antibodies[antibody_name]
        
        # Check sequences match CSV
        sequence_mismatches = []
        for header, seq in py_sequences.items():
            if seq == 'NA':
                continue
            header_lower = header.lower()
            
            # Determine which CSV sequence to compare against
            csv_seq = None
            if 'heavy' in header_lower:
                if '_2' in header or 'chain_2' in header_lower:
                    csv_seq = csv_data.heavy2
                else:
                    csv_seq = csv_data.heavy1
            elif 'light' in header_lower:
                if '_2' in header or 'chain_2' in header_lower:
                    csv_seq = csv_data.light2
                else:
                    csv_seq = csv_data.light1
            
            # Compare sequences
            if csv_seq and seq != csv_seq:
                sequence_mismatches.append({
                    'chain': header,
                    'py_len': len(seq),
                    'csv_len': len(csv_seq)
                })
        
        if sequence_mismatches:
            invalid_sequences.append({
                'name': antibody_name,
                'error': f"Sequence mismatch: {len(sequence_mismatches)} chain(s) differ from CSV",
                'mismatches': sequence_mismatches
            })
            continue
        
        # Validate amino acids
        invalid_chains = []
        for header, seq in py_sequences.items():
            if not validate_sequence(seq):
                invalid_chains.append({'header': header, 'sequence': seq[:50]})
        
        if invalid_chains:
            invalid_sequences.append({'name': antibody_name, 'invalid_chains': invalid_chains})
            continue
        
        matched.append({
            'name': antibody_name,
            'chain_count': count_chains(py_sequences),
            'format': csv_data.format,
            'status': csv_data.status
        })
    
    return ValidationResult(matched=matched, not_in_csv=not_in_csv, invalid_sequences=invalid_sequences)


def print_report(results: ValidationResult):
    # Prints detailed validation report to console
    print("=" * 70)
    print("VALIDATION REPORT")
    print("=" * 70)
    
    print(f"\nMATCHED: {len(results.matched)} antibodies validated")
    if results.matched:
        for item in results.matched[:5]:
            print(f"  - {item['name']} ({item['chain_count']} chains, {item['format']})")
        if len(results.matched) > 5:
            print(f"  ... and {len(results.matched) - 5} more")
    
    print(f"\nNOT IN CSV: {len(results.not_in_csv)} antibodies")
    if results.not_in_csv:
        for item in results.not_in_csv:
            print(f"  - {item['name']} (Path: {item['path']})")
    
    print(f"\nINVALID SEQUENCES: {len(results.invalid_sequences)} antibodies")
    if results.invalid_sequences:
        for item in results.invalid_sequences:
            print(f"  - {item['name']}: {item.get('error', 'Invalid amino acids')}")
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Validated: {len(results.matched)}")
    print(f"  Not in CSV: {len(results.not_in_csv)}")
    print(f"  Invalid: {len(results.invalid_sequences)}")
    print(f"  Total issues: {len(results.not_in_csv) + len(results.invalid_sequences)}")
    
    print("\n  Files by folder:")
    total = 0
    for subfolder in CATEGORY_FOLDERS:
        folder_path = ANTIBODY_FOLDER / subfolder
        if folder_path.exists():
            count = len([f for f in folder_path.glob('*.py') if f.name not in EXCLUDE_FILES])
            total += count
            print(f"    - {subfolder}: {count}")
    print(f"    ----------------------")
    print(f"    - Total: {total}")


def save_report(results: ValidationResult, output_file: Path):
    # Saves validation results to CSV file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Antibody', 'Status', 'Chain_Count', 'Format', 'Clinical_Status', 'Notes'])
        for item in results.matched:
            writer.writerow([item['name'], 'VALIDATED', item['chain_count'], item['format'], item['status'], ''])
        for item in results.not_in_csv:
            writer.writerow([item['name'], 'NOT_IN_CSV', item['chain_count'], '', '', item.get('path', '')])
        for item in results.invalid_sequences:
            notes = item.get('error', 'Invalid amino acids')
            writer.writerow([item['name'], 'INVALID', 0, '', '', notes])
    print(f"\nReport saved to: {output_file.name}")


def main():
    # Main entry point - runs validation workflow
    print("=" * 70)
    print("ANTIBODY SEQUENCE VALIDATION")
    print("=" * 70)
    print(f"CSV: {CSV_PATH.name}")
    print(f"Folder: {ANTIBODY_FOLDER.name}\n")
    
    if not CSV_PATH.exists():
        print(f"ERROR: CSV not found: {CSV_PATH}")
        return
    if not ANTIBODY_FOLDER.exists():
        print(f"ERROR: Folder not found: {ANTIBODY_FOLDER}")
        return
    
    csv_antibodies = load_csv_data(CSV_PATH)
    print(f"Loaded {len(csv_antibodies)} therapeutics from CSV")
    
    results = validate_antibody_files(ANTIBODY_FOLDER, csv_antibodies)
    print_report(results)
    save_report(results, ANTIBODY_FOLDER / 'validation_report.csv')
    
    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()