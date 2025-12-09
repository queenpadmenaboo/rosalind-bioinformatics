"""
Categorize antibody Python files into 4 subdirectories based on format.

FOLDER STRUCTURE:
    multispecific_antibodies/
    ├── Bispecific_mAb/
    ├── Bispecific_scFv/
    ├── Other_Formats/
    ├── Whole_mAb/
    ├── calculate_features.py
    ├── categorize_antibody_format.py
    ├── readme_count.py
    ├── README.md
    ├── sabdabconverter.py
    ├── selenium_antibody_scraper.py
    ├── thera_sabdab_scraper.py
    ├── validate_antibody_sequences.py
    └── validation_report.csv

PURPOSE:
    - Reads format from TheraSAbDab CSV for each antibody file in root folder
    - Categorizes into: Whole_mAb, Bispecific_mAb, Bispecific_scFv, or Other_Formats
    - If not in CSV, file is flagged and NOT moved
    - Moves files to appropriate subfolder after confirmation

CATEGORY LOGIC:
    - Whole_mAb: format contains "whole mab" or "whole ab"
    - Bispecific_mAb: format contains "bispecific mab", "bispecific igg", etc.
    - Bispecific_scFv: format contains "t-cell engager", "bite", "scfv", etc.
    - Other_Formats: everything else (trispecifics, tetraspecifics, nanobodies, fabs)

OUTPUT:
    - Console report showing categorization
    - Files moved to subfolders after user confirms

USAGE:
    python categorize_antibody_format.py
"""

import pandas as pd
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ============================
# CONFIGURATION
# ============================

CSV_PATH = Path(r"C:\Users\bunsr\TheraSAbDab_SeqStruc_07Dec2025.csv")
ANTIBODY_FOLDER = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies")

CATEGORY_KEYWORDS = {
    'Whole_mAb': ['whole mab', 'whole ab'],
    'Bispecific_mAb': [
        'bispecific mab', 'bispecific mixed mab', 'bispecific igg',
        'bispecific dual variable domain', 'bispecific mixed format',
        'bispecific (h-gamma1', 'bispecific (vl-vh',
        'bispecific (l-kappa-h-gamma1', 'bispecific (half'
    ],
    'Bispecific_scFv': [
        'bispecific t-cell engager', 'bite', 'tce',
        'bispecific single domains', 'bispecific mixed domains',
        'single domain variable fragment'
    ]
}

EXCLUDE_FILES = {
    '__init__.py', 'categorize_antibody_format.py', 'categorize_simple.py',
    'therasabdab_analyze_formats.py', 'validate_antibody_sequences.py',
    'readme_count.py', 'sabdabconverter.py', 'selenium_antibody_scraper.py',
    'thera_sabdab_scraper.py', 'calculate_features.py'
}

# ============================
# FUNCTIONS
# ============================

def get_category(format_string: Optional[str]) -> str:
    if not format_string or pd.isna(format_string):
        return 'Other_Formats'
    
    format_lower = format_string.lower()
    
    for keyword in CATEGORY_KEYWORDS['Whole_mAb']:
        if keyword in format_lower:
            return 'Whole_mAb'
    
    for keyword in CATEGORY_KEYWORDS['Bispecific_mAb']:
        if keyword in format_lower:
            return 'Bispecific_mAb'
    
    for keyword in CATEGORY_KEYWORDS['Bispecific_scFv']:
        if keyword in format_lower:
            return 'Bispecific_scFv'
    
    return 'Other_Formats'


def load_csv(csv_path: Path) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_csv(csv_path)
        if 'Therapeutic' not in df.columns or 'Format' not in df.columns:
            print("ERROR: CSV missing Therapeutic or Format columns")
            return None
        print(f"Loaded {len(df)} entries from CSV")
        return df
    except Exception as e:
        print(f"ERROR reading CSV: {e}")
        return None


def get_antibody_files(folder: Path) -> List[str]:
    py_files = [f.name for f in folder.glob('*.py') if f.is_file() and f.name not in EXCLUDE_FILES]
    print(f"Found {len(py_files)} antibody files in root folder")
    if py_files:
        print(f"  Files: {', '.join(py_files)}")
    return py_files


def categorize_files(py_files: List[str], df: pd.DataFrame) -> Dict[str, List[Tuple[str, str, str]]]:
    results = defaultdict(list)
    
    for py_file in py_files:
        name = py_file.replace('.py', '')
        match = df[df['Therapeutic'].str.lower().str.strip() == name.lower().strip()]
        
        if not match.empty:
            format_type = match.iloc[0]['Format']
            category = get_category(format_type)
            results[category].append((name, format_type, py_file))
        else:
            results['NOT_IN_CSV'].append((name, 'Not in CSV', py_file))
    
    return dict(results)


def print_results(results: Dict[str, List[Tuple[str, str, str]]]):
    print("\n" + "=" * 70)
    print("CATEGORIZATION RESULTS")
    print("=" * 70)
    
    for category in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']:
        if category in results:
            items = results[category]
            print(f"\n{category}: {len(items)} files")
            for name, fmt, _ in sorted(items):
                print(f"  - {name}: {fmt}")
    
    if 'NOT_IN_CSV' in results and results['NOT_IN_CSV']:
        print(f"\nNOT IN CSV (will not be moved): {len(results['NOT_IN_CSV'])} files")
        for name, _, _ in results['NOT_IN_CSV']:
            print(f"  - {name}")


def organize_files(py_folder: Path, results: Dict[str, List[Tuple[str, str, str]]]):
    print("\n" + "=" * 70)
    print("MOVING FILES")
    print("=" * 70)
    
    for category in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']:
        if category not in results:
            continue
        
        category_folder = py_folder / category
        category_folder.mkdir(parents=True, exist_ok=True)
        
        for name, fmt, py_file in results[category]:
            src = py_folder / py_file
            dst = category_folder / py_file
            try:
                shutil.move(str(src), str(dst))
                print(f"  Moved {py_file} -> {category}/")
            except Exception as e:
                print(f"  ERROR moving {py_file}: {e}")


def main():
    print("=" * 70)
    print("ANTIBODY FILE CATEGORIZATION")
    print("=" * 70)
    print(f"CSV: {CSV_PATH.name}")
    print(f"Folder: {ANTIBODY_FOLDER.name}\n")
    
    df = load_csv(CSV_PATH)
    if df is None:
        return
    
    py_files = get_antibody_files(ANTIBODY_FOLDER)
    if not py_files:
        print("\nNo files to categorize in root folder.")
        return
    
    results = categorize_files(py_files, df)
    print_results(results)
    
    movable = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']
    files_to_move = sum(len(results.get(cat, [])) for cat in movable)
    
    if files_to_move == 0:
        print("\n" + "=" * 70)
        print("No files to move.")
        print("=" * 70)
        return
    
    print("\n" + "=" * 70)
    response = input("Proceed with moving files? (yes/no): ").strip().lower()
    if response != 'yes':
        print("Cancelled.")
        return
    
    organize_files(ANTIBODY_FOLDER, results)
    
    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()