"""
Validates antibody Python files against TheraSAbDab CSV database.
"""

import csv
from pathlib import Path
from typing import Dict, List, Set
from dataclasses import dataclass

# ============================
# CONFIGURATION
# ============================

CSV_PATH = Path(r"C:\Users\bunsr\TheraSAbDab_SeqStruc_07Dec2025.csv")
ANTIBODY_FOLDER = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")

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
    'therasabdab_analyze_formats.py'
}

# ============================
# DATA STRUCTURES
# ============================

@dataclass
class AntibodyData:
    name: str
    heavy1: str
    light1: str
    heavy2: str
    light2: str
    format: str
    status: str

@dataclass
class ValidationResult:
    matched: List[Dict]
    not_in_csv: List[Dict]
    invalid_sequences: List[Dict]

# ============================
# FUNCTIONS
# ============================

def load_csv_data(csv_file: Path) -> Dict[str, AntibodyData]:
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
    seq_clean = seq.strip().upper()
    if seq_clean in ('NA', 'N/A', ''):
        return True
    return all(aa in VALID_AA for aa in seq_clean)


def count_chains(sequences: Dict[str, str]) -> int:
    return sum(1 for seq in sequences.values() if seq not in ('NA', ''))


def get_all_antibody_files(py_folder: Path) -> List[Path]:
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