import os

# Paths
repo_folder = r"C:\Users\bunsr\rosalind-bioinformatics\bioinformatics_stronghold"
table_file = os.path.join(repo_folder, "rosalind_table.md")
readme_file = os.path.join(repo_folder, "README.md")

# Read the table
with open(table_file, "r", encoding="utf-8") as f:
    table_content = f.read()

# Static README header and progress tracker
readme_header = """# 🧬 Rosalind-Bioinformatics

Tracking my solutions to Rosalind’s challenges across multiple tracks:

- **Python Village** – Basic programming problems with Python  
- **Bioinformatics Stronghold** – Core bioinformatics problems  
- **Bioinformatics Armory** – Advanced bioinformatics problems  

---

## 📊 Progress Tracker

| Status | Description |
|--------|-------------|
| ✅ | Completed problems automatically added via `generate_table.py` |
| 📅 | Each script includes the completion date and problem title |
| ⚙️ | Markdown table sorted chronologically |

---
"""

# Combine header + table
full_readme = f"{readme_header}\n{table_content}"

# Write README.md
with open(readme_file, "w", encoding="utf-8") as f:
    f.write(full_readme)

print(f"✅ README.md updated dynamically with {table_file}")