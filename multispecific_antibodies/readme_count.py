import os
import re
from datetime import datetime
from collections import defaultdict

EXCLUDE_FILES = {
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'validate_antibody_sequences.py', 'validation_report.csv',
    'categorize_antibody_format.py', 'fix_sequences.py',
    'therasabdab_analyze_formats.py', 'calculate_features.py',
    'sequence_features.csv', 'sequence_features.xlsx', 'ml_sasa_predictor.py',
    'all_antibody_sasa_features.csv', 'ml_sasa_predictor_chains.py', 'all_antibody_sasa_chains.csv',
    '3D_structure_builder.py'
}
def update_readme_count():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    print("--------------------------------------------------")
    print(f"Running readme_count.py from: {script_dir}")
    print("Scanning ONLY existing antibody folders.")
    print("--------------------------------------------------")

    # Antibody folders to scan
    antibody_folders = [
        "Bispecific_mAb",
        "Bispecific_scFv",
        "Main",
        "Other_Formats",
        "Whole_mAb"
    ]

    # Track files by folder
    folder_counts = defaultdict(list)

    # Scan folders and count .py files
    for folder in antibody_folders:
        folder_path = os.path.join(script_dir, folder)
        if not os.path.exists(folder_path):
            print(f" Folder not found, skipping: {folder}")
            continue

        for filename in os.listdir(folder_path):
            if filename.endswith(".py") and filename not in EXCLUDE_FILES:
                folder_counts[folder].append(filename)

    # Display breakdown by folder
    print("\n" + "=" * 60)
    print("ANTIBODY FILES BY FOLDER")
    print("=" * 60)

    total = 0
    for folder in sorted(folder_counts.keys()):
        count = len(folder_counts[folder])
        total += count
        print(f"\n{folder}: {count} files")
        for fname in sorted(folder_counts[folder])[:3]:
            print(f"  • {fname}")
        if count > 3:
            print(f"  ... and {count - 3} more")

    print("\n" + "=" * 60)
    print(f"TOTAL ANTIBODIES: {total}")
    print("=" * 60)

    # Update README.md
    readme_path = os.path.join(script_dir, "README.md")
    if not os.path.exists(readme_path):
        print("\n ERROR: README.md not found!")
        return

    with open(readme_path, "r", encoding="utf-8") as f:
        content = f.read()

    current_date = datetime.now().strftime("%B %Y")

    updated = re.sub(
        r"- Data last updated:.*\n- Total antibodies:.*",
        f"- Data last updated: {current_date}\n- Total antibodies: {total}",
        content
    )

    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(updated)

    print("\n README.md updated successfully!")
    print(f"New total: {total}")
    print(f"Updated date: {current_date}")
    print("--------------------------------------------------")

if __name__ == "__main__":
    update_readme_count()
