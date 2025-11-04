"""
Author: Bunsree Patel
Date: November 3, 2025
"""

import os
import re
from datetime import datetime

# Path to your Rosalind repo folder
rosalind_folder = r"C:\Users\bunsr\rosalind-bioinformatics\bioinformatics_stronghold"

# Output Markdown file inside that folder
output_file = os.path.join(rosalind_folder, "rosalind_table.md")

# Regex that matches:
# - October 24, 2025
# - Oct 24, 2025
# - 2025-10-24
# - 10/24/2025
date_pattern = re.compile(
    r'([A-Za-z]+ \d{1,2}, \d{4}|\d{4}-\d{2}-\d{2}|\d{2}/\d{2}/\d{4})'
)

results = []

def normalize_date(date_str):
    """Convert any supported date string to a sortable datetime object."""
    for fmt in ("%B %d, %Y", "%b %d, %Y", "%Y-%m-%d", "%m/%d/%Y"):
        try:
            return datetime.strptime(date_str, fmt)
        except ValueError:
            continue
    return None

# Loop through all Python files in the repo folder
for file_name in os.listdir(rosalind_folder):
    if file_name.endswith(".py") and not file_name.startswith("__"):
        file_path = os.path.join(rosalind_folder, file_name)
        with open(file_path, "r", encoding="utf-8") as f:
            lines = f.readlines()[:10]
            date_found = "Unknown"
            for line in lines:
                match = date_pattern.search(line)
                if match:
                    date_found = match.group()
                    break

            norm_date = normalize_date(date_found)
            results.append((norm_date, date_found, file_name.replace(".py", "")))

# Sort oldest to newest, with Unknown last
results.sort(key=lambda x: (x[0] is None, x[0]))

# Write Markdown table
with open(output_file, "w", encoding="utf-8") as f:
    f.write("# Rosalind Problem Tracker\n\n")
    f.write("| Date | Problem | Completed |\n")
    f.write("|------|---------|-----------|\n")
    for norm_date, date_str, problem in results:
        f.write(f"| {date_str} | {problem} | [x] |\n")

print(f"\n✅ Table generated and written to {output_file}")