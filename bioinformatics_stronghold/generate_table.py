import os
import re

# Path to your local Rosalind repo folder
rosalind_folder = r"C:\Users\bunsr\rosalind-bioinformatics"

# Path to output Markdown file (anywhere you like)
output_file = r"C:\Users\bunsr\Documents\rosalind_tracker.md"
# ===========================

# Regex to match dates in YYYY-MM-DD or MM/DD/YYYY
date_pattern = re.compile(r'(\d{4}-\d{2}-\d{2}|\d{2}/\d{2}/\d{4})')

results = []

# Loop through all Python files in the Rosalind folder
for file_name in os.listdir(rosalind_folder):
    if file_name.endswith(".py"):
        file_path = os.path.join(rosalind_folder, file_name)
        with open(file_path, "r", encoding="utf-8") as f:
            # Check first 10 lines for date
            lines = f.readlines()[:10]
            date_found = None
            for line in lines:
                match = date_pattern.search(line)
                if match:
                    date_found = match.group()
                    break
            if not date_found:
                date_found = "Unknown"
            # Extract problem name from filename (remove .py)
            problem_name = file_name.replace(".py", "")
            results.append((date_found, problem_name))

# Sort results by date (Unknown dates go last)
results.sort(key=lambda x: x[0])

# Write to Markdown file
with open(output_file, "w", encoding="utf-8") as f:
    f.write("# Rosalind Problem Tracker\n\n")
    f.write("| Date | Problem | Completed |\n")
    f.write("|------|---------|-----------|\n")
    for date, problem in results:
        f.write(f"| {date} | {problem} | [x] |\n")  # [x] = completed
        
print(f"\n✅ Table generated and written to {output_file}")