import os
import re
from datetime import datetime
from collections import defaultdict

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
        "validation_report.csv"
    }

    # Track files by folder
    folder_counts = defaultdict(list)

    # ------------------------------------------------
    # Walk the entire directory tree starting at "."
    # ------------------------------------------------
    for root, dirs, files in os.walk(script_dir):
        for filename in files:
            if filename.endswith(".py") and filename not in excluded:
                full_path = os.path.join(root, filename)
                
                # Determine folder category
                rel_path = os.path.relpath(root, script_dir)
                if rel_path == ".":
                    folder_name = "Main"
                else:
                    folder_name = rel_path.split(os.sep)[0]
                
                folder_counts[folder_name].append(filename)

    # ------------------------------------------------
    # Display breakdown by folder
    # ------------------------------------------------
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

    print("\n✓ README.md updated successfully!")
    print(f"New total: {total}")
    print(f"Updated date: {current_date}")
    print("--------------------------------------------------")


if __name__ == "__main__":
    update_readme_count()