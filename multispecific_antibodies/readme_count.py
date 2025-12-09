import os
import re
from datetime import datetime

def update_readme_count():
    # ------------------------------------------------
    # Ensure script runs relative to its own directory
    # ------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    print("--------------------------------------------------")
    print(f"Running readme_count.py from: {script_dir}")
    print("Recursively scanning ALL subfolders under this directory.")
    print("--------------------------------------------------")

    # ------------------------------------------------
    # Exclude utility / helper scripts
    # ------------------------------------------------
    excluded = {
        "readme_count.py",
        "sabdabconverter.py",
        "selenium_antibody_scraper.py",
        "thera_sabdab_scraper.py",
        "validate_antibody_sequences.py",
        "categorize_antibody_format.py",
        "validate_antibody_sequences.py",
        "validation_report.csv"
    }

    antibody_files = []

    # ------------------------------------------------
    # Walk the entire directory tree starting at "."
    # ------------------------------------------------
    for root, dirs, files in os.walk(script_dir):
        for filename in files:
            if filename.endswith(".py") and filename not in excluded:
                full_path = os.path.join(root, filename)
                antibody_files.append(full_path)

    # ------------------------------------------------
    # Display summary
    # ------------------------------------------------
    print("\nDetected antibody .py files:")
    for path in antibody_files:
        print(" -", path)

    total = len(antibody_files)
    print(f"\nTotal antibodies counted: {total}")

    # ------------------------------------------------
    # Update README.md
    # ------------------------------------------------
    readme_path = os.path.join(script_dir, "README.md")

    if not os.path.exists(readme_path):
        print("\n❌ ERROR: README.md not found!")
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

    print("\nREADME.md updated successfully!")
    print(f"New total: {total}")
    print(f"Updated date: {current_date}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    update_readme_count()