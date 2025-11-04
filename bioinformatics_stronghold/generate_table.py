"""
Author: Bunsree Patel
Date: November 3, 2025
"""

import os
import re

# Path to local Rosalind repo folder
rosalind_folder = r"C:\Users\bunsr\rosalind-bioinformatics\bioinformatics_stronghold"

# Output Markdown file inside the repo
output_file = os.path.join(rosalind_folder, "rosalind_table.md")

# Regex to match dates like "October 24, 2025" or "Nov 3, 2025"
date_pattern = re.compile(r'([A-Za-z]+\s+\d{1,2},\s+\d{4})')

results = []

# Loop through all Python files in the folder
for file_name in os.listdir(rosalind_folder):
    if file_name.endswith(".py") and not file_name.startswith("__") and file_name not in ["generate_table.py", "main.py", "update_readme.py"]:
        file_path = os.path.join(rosalind_folder, file_name)
        with open(file_path, "r", encoding="utf-8") as f:
            # Check first 10 lines for a date
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

# Helper: turn "October 24, 2025" into sortable tuple (year, month, day)
def sort_key(item):
    date_str, _ = item
    months = {
        'January': 1, 'February': 2, 'March': 3, 'April': 4, 'May': 5,
        'June': 6, 'July': 7, 'August': 8, 'September': 9,
        'October': 10, 'November': 11, 'December': 12
    }
    if date_str == "Unknown":
        return (9999, 99, 99)
    try:
        month_name, day_year = date_str.split(" ", 1)
        day, year = day_year.replace(",", "").split()
        return (int(year), months.get(month_name, 0), int(day))
    except Exception:
        return (9999, 99, 99)

# Sort ascending
results.sort(key=sort_key)

# Write Markdown
with open(output_file, "w", encoding="utf-8") as f:
    f.write("# Rosalind Problem Tracker\n\n")
    f.write("| Date | Problem | Completed |\n")
    f.write("|------|---------|-----------|\n")
    for date, problem in results:
        f.write(f"| {date} | {problem} | [x] |\n")

print(f"\n✅ Rosalind table updated: {output_file}")