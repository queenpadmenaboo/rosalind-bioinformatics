import os       # Provides a way to interact with operating system. Makes script portable
import sys      # Gives access to system-specific parameters and functions. Provides proper termination when error occurs.

# Paths
repo_root = r"C:\Users\bunsr\rosalind-bioinformatics"                   # Root of the repository
stronghold_folder = os.path.join(repo_root, "bioinformatics_stronghold")
table_file = os.path.join(stronghold_folder, "rosalind_table.md")       # Table source
readme_file = os.path.join(repo_root, "README.md")                      # Output README at root level

# Check if table file exists
if not os.path.exists(table_file):
    print(f"❌ Error: {table_file} not found!")
    sys.exit(1)

# Read the table
with open(table_file, "r", encoding="utf-8") as f:
    table_content = f.read()

# Static README header and progress tracker
readme_header = """# 🧬 Rosalind-Bioinformatics

My solutions to Rosalind bioinformatics problems:

- **Python Village** - Python programming fundamentals
- **Bioinformatics Stronghold** - Core bioinformatics algorithms
- **Bioinformatics Armory** - Advanced bioinformatics techniques

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
full_readme = f"{readme_header}\n{table_content}\n"

# Write README.md
with open(readme_file, "w", encoding="utf-8") as f:
    f.write(full_readme)

print(f"✅ README.md updated at {readme_file}")
print(f"📄 Table source: {table_file}")