"""
Author: Bunsree Patel
Date: November 3, 2025
"""
``
import os
import re
import time

# Path to your Rosalind repo folder
rosalind_folder = r"C:\Users\bunsr\rosalind-bioinformatics\bioinformatics_stronghold"

# Output Markdown file inside the repo
output_file = os.path.join(rosalind_folder, "rosalind_table.md")

# Regex to detect human-readable or numeric dates
date_pattern = re.compile(r'([A-Za-z]+\s\d{1,2},\s\d{4}|\d{4}-\d{2}-\d{2}|\d{2}/\d{2}/\d{4})')

def generate_table():
    """Scans .py files and updates the Markdown tracker."""
    results = []
    for file_name in os.listdir(rosalind_folder):
        if file_name.endswith(".py") and not file_name.startswith("__"):
            file_path = os.path.join(rosalind_folder, file_name)
            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()[:10]
                date_found = None
                for line in lines:
                    match = date_pattern.search(line)
                    if match:
                        date_found = match.group()
                        break
                if not date_found:
                    date_found = "Unknown"
                problem_name = file_name.replace(".py", "")
                results.append((date_found, problem_name))

    # Sort by month names / human-readable dates first, then numeric
    def date_sort_key(item):
        date = item[0]
        if date == "Unknown":
            return (9999, 99, 99)
        try:
            return time.strptime(date, "%B %d, %Y")
        except ValueError:
            return time.strptime("December 31, 9999", "%B %d, %Y")

    results.sort(key=date_sort_key)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("# Rosalind Problem Tracker\n\n")
        f.write("| Date | Problem | Completed |\n")
        f.write("|------|---------|-----------|\n")
        for date, problem in results:
            f.write(f"| {date} | {problem} | [x] |\n")

    print(f"✅ Rosalind table updated: {output_file}")

# --- Auto-refresh loop ---
print("🔁 Watching for file changes... Press Ctrl+C to stop.\n")
prev_state = {}
while True:
    current_state = {
        f: os.path.getmtime(os.path.join(rosalind_folder, f))
        for f in os.listdir(rosalind_folder) if f.endswith(".py")
    }
    if current_state != prev_state:
        generate_table()
        prev_state = current_state
    time.sleep(3)