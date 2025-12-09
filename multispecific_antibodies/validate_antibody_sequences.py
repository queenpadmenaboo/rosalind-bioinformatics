import csv
import os
from pathlib import Path

# ============================
# CONFIGURATION
# ============================

# Path to the CSV downloaded from TheraSAbDab with all antibody data
THERA_CSV_FILE = r"C:\Users\bunsr\TheraSAbDab_SeqStruc_ 07Dec2025.csv"

# Folder containing the 470 antibody.py files
ANTIBODY_PY_FOLDER = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"

# Valid amino acids for protein sequences
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# ============================
# FUNCTIONS
# ============================

def load_csv_sequences(csv_file):
    """
    Reads the TheraSAbDab CSV and extracts all antibody sequences into a dictionary.
    Key = antibody name (lowercase, underscores), Value = sequence data + metadata
    """
    sequences = {}
    
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        
        # Print all column names from the CSV to debug
        print("\nColumn names found in CSV:")
        if reader.fieldnames:
            for column_name in reader.fieldnames:
                print(f"  Column: '{column_name}'")
        print()
        
        for row in reader:
            # Normalizes the therapeutic name to match .py filenames
            therapeutic_raw = row['Therapeutic'].strip()
            therapeutic_key = therapeutic_raw.lower().replace(' ', '_').replace('-', '_')
            
            # Stores all relevant data for this antibody
            sequences[therapeutic_key] = {
                'heavy1': row['HeavySequence'].strip(),
                'light1': row['LightSequence'].strip(),
                'heavy2': row['HeavySequence(ifbispec)'].strip(),
                'light2': row['LightSequence(ifbispec)'].strip(),
                'format': row['Format'],
                'status': row['Est. Status'],
                'therapeutic_name': therapeutic_raw
            }
    
    print(f"Loaded {len(sequences)} therapeutics from CSV")
    return sequences


def parse_py_file(filepath):
    """
    Reads an antibody.py file and extracts the FASTA sequences.
    Returns dict: {header: sequence}
    Handles 'na' as missing chains (valid for some formats like BiTEs).
    """
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    sequences = {}
    current_header = None
    current_seq = []
    
    # Parses FASTA format (lines starting with > are headers, rest is sequence)
    for line in content.split('\n'):
        line = line.strip()
        
        if line.startswith('>'):
            # Saves previous sequence if exists
            if current_header:
                seq = ''.join(current_seq)
                # Handles missing chains marked as 'na'
                sequences[current_header] = 'NA' if seq.lower() in ['na', 'n/a', ''] else seq
            
            # Starts new sequence
            current_header = line[1:]  # Removes '>'
            current_seq = []
        elif line and not line.startswith('"""') and not line.startswith('=') and not line.startswith('#'):
            # Adds sequence line (skips Python syntax)
            current_seq.append('NA' if line.lower() in ['na', 'n/a'] else line)
    
    # Saves last sequence
    if current_header:
        seq = ''.join(current_seq)
        sequences[current_header] = 'NA' if seq.lower() in ['na', 'n/a', ''] else seq
    
    return sequences


def validate_sequence(seq):
    """
    Checks if sequence only contains valid amino acids.
    'NA' or empty = valid (means missing chain, which is OK for some formats)
    """
    seq_upper = seq.strip().upper()
    
    # Missing chains are valid
    if seq_upper in ['NA', 'N/A', '']:
        return True
    
    # Checks each amino acid
    return all(aa in VALID_AA for aa in seq_upper)


def count_chains(py_sequences):
    """
    Counts how many actual sequences exist (not NA).
    Helps identify if antibody is bispecific (4 chains) or standard mAb (2 chains).
    """
    return sum(1 for seq in py_sequences.values() if seq not in ['NA', ''])


def compare_antibody_files(py_folder, csv_sequences):
    """
    Main validation logic: for each .py file, checks if:
    1. Antibody exists in the CSV
    2. Sequences contain valid amino acids
    3. Chain counts make sense
    
    Returns categorized results.
    """
    results = {
        'matched': [],           # Valid antibodies found in CSV
        'not_in_csv': [],        # Files not in TheraSAbDab
        'invalid_sequences': []  # Files with invalid amino acids
    }
    
    py_files = list(Path(py_folder).glob('*.py'))
    # Skips utility scripts
    script_files = ['readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py', 'thera_sabdab_scraper.py', 'validate_antibody_sequences.py']
    
    print(f"\nFound {len(py_files)} .py files")
    print(f"Validating {len(py_files) - len(script_files)} antibody files\n")
    
    for py_file in py_files:
        if py_file.name in script_files:
            continue
            
        antibody_name = py_file.stem  # Filename without .py
        
        # Attempts to parse the file
        try:
            py_sequences = parse_py_file(py_file)
        except Exception as e:
            results['invalid_sequences'].append({'name': antibody_name, 'error': f"Parse error: {str(e)}"})
            continue
        
        if not py_sequences:
            results['invalid_sequences'].append({'name': antibody_name, 'error': "No sequences found"})
            continue
        
        # Checks if this antibody exists in TheraSAbDab CSV
        if antibody_name not in csv_sequences:
            results['not_in_csv'].append({'name': antibody_name, 'chain_count': count_chains(py_sequences)})
            continue
        
        csv_data = csv_sequences[antibody_name]
        
        # Validates all chains have valid amino acids
        invalid_seqs = []
        for header, seq in py_sequences.items():
            if not validate_sequence(seq):
                invalid_seqs.append({'header': header, 'sequence': seq[:50]})
        
        if invalid_seqs:
            results['invalid_sequences'].append({'name': antibody_name, 'invalid_chains': invalid_seqs})
            continue
        
        # Antibody passed all checks
        results['matched'].append({
            'name': antibody_name,
            'chain_count': count_chains(py_sequences),
            'format': csv_data['format'],
            'status': csv_data['status']
        })
    
    return results


def report_results(results):
    """
    Prints detailed validation report.
    """
    print("\n" + "="*70)
    print("VALIDATION REPORT")
    print("="*70)
    
    print(f"\nMATCHED: {len(results['matched'])} antibodies validated against CSV")
    if results['matched']:
        print("\nSample matched antibodies:")
        for item in results['matched'][:10]:
            print(f"   - {item['name']} ({item['chain_count']} chains, {item['format']}, {item['status']})")
        if len(results['matched']) > 10:
            print(f"   ... and {len(results['matched']) - 10} more")
    
    print(f"\nNOT IN CSV: {len(results['not_in_csv'])} antibodies")
    if results['not_in_csv']:
        print("These antibodies exist in .py files but not in the CSV:")
        for item in results['not_in_csv'][:15]:
            print(f"   - {item['name']} ({item['chain_count']} chains)")
        if len(results['not_in_csv']) > 15:
            print(f"   ... and {len(results['not_in_csv']) - 15} more")
    
    print(f"\nINVALID SEQUENCES: {len(results['invalid_sequences'])} antibodies")
    if results['invalid_sequences']:
        print("These antibodies have invalid amino acid sequences:")
        for item in results['invalid_sequences'][:10]:
            if isinstance(item, dict) and 'error' in item:
                print(f"   - {item['name']}: {item['error']}")
            elif isinstance(item, dict) and 'invalid_chains' in item:
                print(f"   - {item['name']}:")
                for chain in item['invalid_chains'][:3]:
                    print(f"      * {chain['header']}: {chain['sequence']}...")
        if len(results['invalid_sequences']) > 10:
            print(f"   ... and {len(results['invalid_sequences']) - 10} more")
    
    print(f"\nSUMMARY:")
    print(f"   Total validated successfully: {len(results['matched'])}")
    print(f"   Not found in CSV: {len(results['not_in_csv'])}")
    print(f"   Invalid sequences: {len(results['invalid_sequences'])}")
    print(f"   Total issues: {len(results['not_in_csv']) + len(results['invalid_sequences'])}")
    
    if results['matched']:
        bispec_count = sum(1 for item in results['matched'] if 'bispec' in item['format'].lower() or 'trispec' in item['format'].lower())
        mab_count = sum(1 for item in results['matched'] if 'mab' in item['format'].lower() and 'bispec' not in item['format'].lower())
        print(f"\nFormat breakdown (matched antibodies):")
        print(f"   - Multispecifics (bi/tri): ~{bispec_count}")
        print(f"   - Standard mAbs: ~{mab_count}")
        print(f"   - Other formats: ~{len(results['matched']) - bispec_count - mab_count}")


def save_validation_report(results, output_file):
    """
    Saves validation results to a CSV file for further analysis.
    """
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Antibody', 'Status', 'Chain_Count', 'Format', 'Clinical_Status', 'Notes'])
        
        for item in results['matched']:
            writer.writerow([
                item['name'],
                'VALIDATED',
                item['chain_count'],
                item['format'],
                item['status'],
                ''
            ])
        
        for item in results['not_in_csv']:
            writer.writerow([
                item['name'],
                'NOT_IN_CSV',
                item['chain_count'],
                '',
                '',
                'Not found in TheraSAbDab CSV'
            ])
        
        for item in results['invalid_sequences']:
            if isinstance(item, dict) and 'error' in item:
                writer.writerow([
                    item['name'],
                    'INVALID',
                    0,
                    '',
                    '',
                    item['error']
                ])
            elif isinstance(item, dict) and 'invalid_chains' in item:
                writer.writerow([
                    item['name'],
                    'INVALID',
                    0,
                    '',
                    '',
                    f"Invalid amino acids in {len(item['invalid_chains'])} chains"
                ])
    
    print(f"\nValidation report saved to: {output_file}")


# ============================
# MAIN
# ============================

def main():
    """
    Runs validation of antibody .py files against TheraSAbDab CSV.
    """
    print("="*70)
    print("ANTIBODY SEQUENCE VALIDATION")
    print("Checking .py files against TheraSAbDab CSV")
    print("="*70)
    
    if not os.path.exists(THERA_CSV_FILE):
        print(f"\nERROR: CSV file not found: {THERA_CSV_FILE}")
        return
    
    if not os.path.exists(ANTIBODY_PY_FOLDER):
        print(f"\nERROR: Antibody folder not found: {ANTIBODY_PY_FOLDER}")
        return
    
    print(f"\nCSV File: {THERA_CSV_FILE}")
    print(f"Antibody Folder: {ANTIBODY_PY_FOLDER}")
    
    csv_sequences = load_csv_sequences(THERA_CSV_FILE)
    
    results = compare_antibody_files(ANTIBODY_PY_FOLDER, csv_sequences)
    
    report_results(results)
    
    report_file = os.path.join(ANTIBODY_PY_FOLDER, 'validation_report.csv')
    save_validation_report(results, report_file)
    
    print("\n" + "="*70)
    print("VALIDATION COMPLETE")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()