"""
Author: Bunsree Patel
Date: November 3, 2025
"""

import os
import re

# Paths to Rosalind problem sections
folders_to_scan = [
    r"C:\Users\bunsr\rosalind-bioinformatics\python_village",
    r"C:\Users\bunsr\rosalind-bioinformatics\bioinformatics_stronghold"
]

# Output Markdown file inside the repo
output_file = os.path.join(folders_to_scan[1], "rosalind_table.md")

# Regex to match dates like "October 24, 2025" or "Nov 3, 2025"
date_pattern = re.compile(r'([A-Za-z]+\s+\d{1,2},\s+\d{4})')

results = []

# Loop through all Python files in the folder
for folder in folders_to_scan:
    for file_name in os.listdir(folder):
        if file_name.endswith(".py") and not file_name.startswith("__") and file_name not in [
            "generate_table.py", "main.py", "update_readme.py",
            "DNAToolkit.py", "FibonacciNumbers.py", "ProteinToolkit.py"
        
        ]:
            file_path = os.path.join(folder, file_name)
                

            with open(file_path, "r", encoding="utf-8") as f:
                lines = f.readlines()[:10]          # Check first 10 lines for a date and location

            # Defaults
            date_found = "Unknown"
            location_found = "Unknown"
            full_name = file_name.replace("py", "")

            # Extract information from header
            for line in lines:
                line = line.strip()

                # Find date
                match_date = date_pattern.search(line)
                if match_date:
                    date_found = match_date.group()

                # Find location information from header
                if "Location:" in line:
                    location_found = line.split("Location:")[1].strip()

                # Full descriptive problem name
                if "Rosalind Problem:" in line:
                    parts = line.split(":", 1)[1].strip() 
                    if " - " in parts:
                        full_name = parts.split(" - ", 1)[1].strip()
                    else:
                        full_name = parts

            results.append((date_found, full_name, location_found))

# Helper: convert "October 24, 2025" into sortable tuple (year, month, day)
def sort_key(item):
    date_str, _, _ = item       # unpack (date, problem, location)
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

# Group results by location
grouped_results = {}
for date, problem, location in results:
    if location not in grouped_results:
        grouped_results[location] = []
    grouped_results[location].append((date, problem, location))

# Sort each location by date
for loc in grouped_results:
    grouped_results[loc].sort(key=sort_key)      # Python’s list.sort() defaults to ascending order.

# Write Markdown table grouped by location
with open(output_file, "w", encoding="utf-8") as f:
    f.write("# Rosalind Problem Tracker\n\n")
    f.write("| Date | Problem | Location | Completed |\n")
    f.write("|------|---------|----------|-----------|\n")
    for location in sorted(grouped_results.keys()):
        for date, problem, loc in grouped_results[location]:
            f.write(f"| {date} | {problem} | {location} | [✅] |\n")

print(f"\n✅ Rosalind table updated (grouped by location): {output_file}")