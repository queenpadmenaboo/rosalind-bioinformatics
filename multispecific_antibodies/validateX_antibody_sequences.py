import csv
import os

# ============================
# CONFIGURATION
# ============================

# Thera-SAbDab CSV file path
THERA_CSV_FILE = r"C:\Users\bunsr\TheraSAbDab_SeqStruc_ 07Dec2025.csv"

# Folder where the individual FASTA-style .py files are saved
OUTPUT_FOLDER = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\thera_fastas"

# Valid amino acids for antibody sequences
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")  # Standard 20 amino acids

# ============================
# FUNCTIONS
# ============================

def load_therasabdab_csv(csv_file):
    """
    Load the Thera-SAbDab CSV file and return a list of rows (as dicts).
    """
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")

    rows = []

    with open(csv_file, newline='', encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(row)

    print(f"✅ Loaded {len(rows)} rows from CSV: {csv_file}")
    return rows


def validate_sequence(seq):
    """
    Check if the amino acid sequence contains only valid amino acids.
    Returns True if valid, False otherwise.
    """
    seq = seq.upper().strip()
    for aa in seq:
        if aa not in VALID_AA:
            return False
    return True


def validate_antibody_sequences(rows):
    """
    Validate the Heavy and Light antibody sequences in the CSV rows.
    Returns a list of tuples: (row_index, row, heavy_valid, light_valid)
    """
    results = []

    for index, row in enumerate(rows):
        # Strip BOM from headers if present
        row = {k.replace("\ufeff", ""): v for k, v in row.items()}

        # Get sequences
        heavy_seq = row.get("HeavySequence", "")
        light_seq = row.get("LightSequence", "")

        # Validate sequences
        heavy_valid = validate_sequence(heavy_seq) if heavy_seq else False
        light_valid = validate_sequence(light_seq) if light_seq else False

        results.append((index, row, heavy_valid, light_valid))

    return results


def report_results(results):
    """
    Print a summary of sequence validation results.
    """
    total = len(results)
    heavy_valid_count = sum(1 for r in results if r[2])
    light_valid_count = sum(1 for r in results if r[3])

    print("\n>> Validation Summary:")
    print(f"Total sequences checked: {total}")
    print(f"Valid heavy sequences: {heavy_valid_count}")
    print(f"Valid light sequences: {light_valid_count}")
    print(f"Invalid heavy sequences: {total - heavy_valid_count}")
    print(f"Invalid light sequences: {total - light_valid_count}")

    # Print first 5 invalid heavy sequences
    print("\n>>> First 5 invalid heavy sequences:")
    count = 0
    for index, row, heavy_valid, light_valid in results:
        if not heavy_valid:
            print(f"Row {index + 1}: {row.get('HeavySequence', '')}")
            count += 1
            if count >= 5:
                break

    # Print first 5 invalid light sequences
    print("\n>>> First 5 invalid light sequences:")
    count = 0
    for index, row, heavy_valid, light_valid in results:
        if not light_valid:
            print(f"Row {index + 1}: {row.get('LightSequence', '')}")
            count += 1
            if count >= 5:
                break


def save_sequences_as_py(results, output_folder):
    """
    Save each valid antibody sequence to a separate .py file
    using HEAVY_CHAIN and LIGHT_CHAIN nomenclature.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for index, row, heavy_valid, light_valid in results:
        if not (heavy_valid and light_valid):
            # Skip sequences that are invalid
            continue

        # Strip BOM for safety
        row = {k.replace("\ufeff", ""): v for k, v in row.items()}

        therapeutic_name = row.get("Therapeutic", f"thera_{index + 1}")
        # Clean name for filename: lowercase, underscores instead of spaces
        filename = therapeutic_name.strip().replace(" ", "_").replace("-", "_").lower()
        filepath = os.path.join(output_folder, f"{filename}.py")

        heavy_seq = row.get("HeavySequence", "")
        light_seq = row.get("LightSequence", "")

        # Write to .py file
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(f'HEAVY_CHAIN = """{heavy_seq}"""\n')
            f.write(f'LIGHT_CHAIN = """{light_seq}"""\n')

    print(f"\n✅ Saved valid sequences to .py files in folder: {output_folder}")


# ============================
# MAIN
# ============================

def main():
    """
    Main function to run the antibody sequence validation and save FASTA-style .py files.
    """
    print(">> Starting Thera-SAbDab antibody sequence validation script...\n")

    # Load CSV
    thera_rows = load_therasabdab_csv(THERA_CSV_FILE)

    # Preview first 3 rows
    print(">>> First 3 rows preview:")
    for i, row in enumerate(thera_rows):
        if i >= 3:
            break
        row_clean = {k.replace("\ufeff", ""): v for k, v in row.items()}
        print(row_clean)

    # Validate sequences
    validation_results = validate_antibody_sequences(thera_rows)

    # Report results
    report_results(validation_results)

    # Save valid sequences as .py files
    save_sequences_as_py(validation_results, OUTPUT_FOLDER)

    print("\n>> Script completed.")


# ============================
# ENTRY POINT
# ============================

if __name__ == "__main__":
    main()
