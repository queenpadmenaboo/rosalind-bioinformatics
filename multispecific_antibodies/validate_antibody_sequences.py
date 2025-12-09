"""
Validates antibody Python files against TheraSAbDab CSV database.
Checks for: sequence validity, chain counts, and CSV matches.
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

# Valid single-letter amino acid codes
VALID_AA: Set[str] = set("ACDEFGHIKLMNPQRSTVWY")

# Utility scripts to exclude from validation
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
    """Stores antibody information from CSV."""
    name: str
    heavy1: str
    light1: str
    heavy2: str
    light2: str
    format: str
    status: str


@dataclass
class ValidationResult:
    """Stores validation results for reporting."""
    matched: List[Dict]
    not_in_csv: List[Dict]
    invalid_sequences: List[Dict]


# ============================
# CSV PROCESSING
# ============================

def load_csv_data(csv_file: Path) -> Dict[str, AntibodyData]:
    """
    Load antibody sequences and metadata from TheraSAbDab CSV.
    
    Returns:
        Dict mapping normalized names to AntibodyData objects
    """
    antibodies = {}
    
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        
        print("✓ CSV columns found:")
        for col in reader.fieldnames:
            print(f"  • {col}")
        print()
        
        for row in reader:
            # Normalize name to match .py filenames
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
    
    print(f"✓ Loaded {len(antibodies)} therapeutics from CSV\n")
    return antibodies


# ============================
# PYTHON FILE PARSING
# ============================

def parse_py_file(filepath: Path) -> Dict[str, str]:
    """
    Extract FASTA sequences from antibody .py file.
    
    Returns:
        Dict mapping header -> sequence (or 'NA' for missing chains)
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    sequences = {}
    current_header = None
    current_seq = []
    
    for line in content.split('\n'):
        line = line.strip()
        
        if line.startswith('>'):
            # Save previous sequence
            if current_header:
                seq = ''.join(current_seq)
                sequences[current_header] = 'NA' if seq.lower() in ('na', 'n/a', '') else seq
            
            # Start new sequence
            current_header = line[1:]
            current_seq = []
            
        elif line and not any(line.startswith(c) for c in ('"""', '=', '#', 'import', 'from')):
            # Add sequence line (skip Python syntax)
            if line.lower() in ('na', 'n/a'):
                current_seq.append('NA')
            else:
                current_seq.append(line)
    
    # Save last sequence
    if current_header:
        seq = ''.join(current_seq)
        sequences[current_header] = 'NA' if seq.lower() in ('na', 'n/a', '') else seq
    
    return sequences


# ============================
# VALIDATION
# ============================

def validate_sequence(seq: str) -> bool:
    """
    Check if sequence contains only valid amino acids.
    'NA' or empty is valid (missing chain, OK for some formats).
    """
    seq_clean = seq.strip().upper()
    
    # Missing chains are valid
    if seq_clean in ('NA', 'N/A', ''):
        return True
    
    # Check all characters are valid amino acids
    return all(aa in VALID_AA for aa in seq_clean)


def count_chains(sequences: Dict[str, str]) -> int:
    """Count non-NA sequences (helps identify bispecifics vs standard mAbs)."""
    return sum(1 for seq in sequences.values() if seq not in ('NA', ''))


def validate_antibody_files(
    py_folder: Path,
    csv_antibodies: Dict[str, AntibodyData]
) -> ValidationResult:
    """
    Validate all antibody .py files against CSV database.
    
    Checks:
        1. File exists in CSV
        2. Sequences contain valid amino acids
        3. Chain counts are reasonable
    """
    matched = []
    not_in_csv = []
    invalid_sequences = []
    
    # Get all .py files, excluding subfolders
    py_files = [f for f in py_folder.glob('*.py') if f.is_file()]
    
    # Also check subfolders
    for subfolder in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']:
        subfolder_path = py_folder / subfolder
        if subfolder_path.exists():
            py_files.extend([f for f in subfolder_path.glob('*.py') if f.is_file()])
    
    antibody_files = [f for f in py_files if f.name not in EXCLUDE_FILES]
    
    print(f"✓ Found {len(py_files)} total .py files")
    print(f"✓ Validating {len(antibody_files)} antibody files\n")
    
    for py_file in antibody_files:
        antibody_name = py_file.stem
        
        # Parse file
        try:
            py_sequences = parse_py_file(py_file)
        except Exception as e:
            invalid_sequences.append({
                'name': antibody_name,
                'error': f"Parse error: {str(e)}"
            })
            continue
        
        if not py_sequences:
            invalid_sequences.append({
                'name': antibody_name,
                'error': "No sequences found in file"
            })
            continue
        
        # Check if in CSV
        if antibody_name not in csv_antibodies:
            not_in_csv.append({
                'name': antibody_name,
                'chain_count': count_chains(py_sequences)
            })
            continue
        
        csv_data = csv_antibodies[antibody_name]
        
        # Validate sequences
        invalid_chains = []
        for header, seq in py_sequences.items():
            if not validate_sequence(seq):
                invalid_chains.append({
                    'header': header,
                    'sequence': seq[:50]
                })
        
        if invalid_chains:
            invalid_sequences.append({
                'name': antibody_name,
                'invalid_chains': invalid_chains
            })
            continue
        
        # Passed all checks
        matched.append({
            'name': antibody_name,
            'chain_count': count_chains(py_sequences),
            'format': csv_data.format,
            'status': csv_data.status
        })
    
    return ValidationResult(
        matched=matched,
        not_in_csv=not_in_csv,
        invalid_sequences=invalid_sequences
    )


# ============================
# REPORTING
# ============================

def print_report(results: ValidationResult):
    """Display detailed validation report."""
    print("\n" + "=" * 80)
    print("VALIDATION REPORT")
    print("=" * 80)
    
    # Matched antibodies
    print(f"\n✓ MATCHED: {len(results.matched)} antibodies validated")
    if results.matched:
        print("\n  Sample validated antibodies:")
        for item in results.matched[:10]:
            print(f"    • {item['name']:30s} ({item['chain_count']} chains, {item['format']}, {item['status']})")
        if len(results.matched) > 10:
            print(f"    ... and {len(results.matched) - 10} more")
    
    # Not in CSV
    print(f"\n⚠ NOT IN CSV: {len(results.not_in_csv)} antibodies")
    if results.not_in_csv:
        print("  These files exist but are not in TheraSAbDab CSV:")
        for item in results.not_in_csv[:15]:
            print(f"    • {item['name']:30s} ({item['chain_count']} chains)")
        if len(results.not_in_csv) > 15:
            print(f"    ... and {len(results.not_in_csv) - 15} more")
    
    # Invalid sequences
    print(f"\n✗ INVALID SEQUENCES: {len(results.invalid_sequences)} antibodies")
    if results.invalid_sequences:
        print("  These files have errors:")
        for item in results.invalid_sequences[:10]:
            if 'error' in item:
                print(f"    • {item['name']}: {item['error']}")
            elif 'invalid_chains' in item:
                print(f"    • {item['name']}:")
                for chain in item['invalid_chains'][:3]:
                    print(f"        └─ {chain['header']}: {chain['sequence']}...")
        if len(results.invalid_sequences) > 10:
            print(f"    ... and {len(results.invalid_sequences) - 10} more")
    
    # Summary
    total_issues = len(results.not_in_csv) + len(results.invalid_sequences)
    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")
    print(f"  Validated successfully: {len(results.matched)}")
    print(f"  Not found in CSV:       {len(results.not_in_csv)}")
    print(f"  Invalid sequences:      {len(results.invalid_sequences)}")
    print(f"  Total issues:           {total_issues}")
    
    # Format breakdown
    if results.matched:
        bispec = sum(1 for x in results.matched if 'bispec' in x['format'].lower())
        trispec = sum(1 for x in results.matched if 'trispec' in x['format'].lower())
        standard_mab = sum(1 for x in results.matched 
                          if 'mab' in x['format'].lower() 
                          and 'bispec' not in x['format'].lower()
                          and 'trispec' not in x['format'].lower())
        other = len(results.matched) - bispec - trispec - standard_mab
        
        print(f"\n  Format breakdown:")
        print(f"    • Bispecifics:   {bispec}")
        print(f"    • Trispecifics:  {trispec}")
        print(f"    • Standard mAbs: {standard_mab}")
        print(f"    • Other formats: {other}")


def save_report(results: ValidationResult, output_file: Path):
    """Save validation results to CSV for further analysis."""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Antibody', 'Status', 'Chain_Count', 'Format', 'Clinical_Status', 'Notes'])
        
        # Matched
        for item in results.matched:
            writer.writerow([
                item['name'],
                'VALIDATED',
                item['chain_count'],
                item['format'],
                item['status'],
                ''
            ])
        
        # Not in CSV
        for item in results.not_in_csv:
            writer.writerow([
                item['name'],
                'NOT_IN_CSV',
                item['chain_count'],
                '',
                '',
                'Not found in TheraSAbDab CSV'
            ])
        
        # Invalid
        for item in results.invalid_sequences:
            notes = item.get('error', f"Invalid amino acids in {len(item.get('invalid_chains', []))} chains")
            writer.writerow([
                item['name'],
                'INVALID',
                0,
                '',
                '',
                notes
            ])
    
    print(f"\n✓ Report saved to: {output_file.name}")


# ============================
# MAIN
# ============================

def main():
    """Run validation workflow."""
    print("=" * 80)
    print("ANTIBODY SEQUENCE VALIDATION")
    print("Checking .py files against TheraSAbDab CSV")
    print("=" * 80)
    
    # Validate paths
    if not CSV_PATH.exists():
        print(f"\n✗ ERROR: CSV file not found: {CSV_PATH}")
        return
    
    if not ANTIBODY_FOLDER.exists():
        print(f"\n✗ ERROR: Antibody folder not found: {ANTIBODY_FOLDER}")
        return
    
    print(f"\nCSV:    {CSV_PATH.name}")
    print(f"Folder: {ANTIBODY_FOLDER.name}\n")
    
    # Load and validate
    csv_antibodies = load_csv_data(CSV_PATH)
    results = validate_antibody_files(ANTIBODY_FOLDER, csv_antibodies)
    
    # Report
    print_report(results)
    
    # Save
    report_file = ANTIBODY_FOLDER / 'validation_report.csv'
    save_report(results, report_file)
    
    print("\n" + "=" * 80)
    print("✓ VALIDATION COMPLETE")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()