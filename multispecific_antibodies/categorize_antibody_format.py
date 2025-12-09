import pandas as pd
import os
import shutil
from collections import defaultdict
from pathlib import Path

# 4-Category mapping based on your data
# ============================================================================
# EDIT 1: Added keywords to catch unusual bispecific formats
# ============================================================================
CATEGORY_MAPPING = {
    'Whole_mAb': [
        'whole mab', 'whole ab'
    ],
    'Bispecific_mAb': [
        'bispecific mixed mab', 'bispecific mab', 'bispecific igg',
        'bispecific dual variable domain ig', 'bispecific mixed format',
        'bispecific (h-gamma1',           # ← ADDED: bexatamig
        'bispecific (vl-vh',              # ← ADDED: other unusual formats
        'bispecific (l-kappa-h-gamma1',   # ← ADDED: fanastomig
        'bispecific (half'                # ← ADDED: rezetamig
    ],
    'Bispecific_scFv': [
        'bispecific t-cell engager', 'bite', 'tce',
        'bispecific single domains', 'bispecific mixed domains',
        'single domain variable fragment'
    ],
    'Other_Formats': [
        'fab', 'nanobody', 'fv fusion',
        'bispecific (vh-vk'  # ← ADDED: acimtamig (VH-VK homodimer)
    ]
}

def get_category(format_string):
    """Get category from format string"""
    if not format_string or pd.isna(format_string):
        return 'Unclassified'
    
    format_lower = format_string.lower()
    
    for category, keywords in CATEGORY_MAPPING.items():
        for keyword in keywords:
            if keyword.lower() in format_lower:
                return category
    
    return 'Unclassified'

def main():
    csv_path = Path(r"C:\Users\meeko\TheraSAbDab_SeqStruc_08Dec2025.csv")
    py_folder = Path(r"C:\Users\meeko\rosalind-bioinformatics\multispecific_antibodies")
    
    print("=" * 80)
    print("CATEGORIZING ANTIBODY FILES INTO 4 SUBFOLDERS")
    print("=" * 80)
    print(f"CSV: {csv_path}")
    print(f"Folder: {py_folder}")
    print()
    
    # Load and validate CSV
    try:
        df = pd.read_csv(csv_path)
        print(f"Loaded {len(df)} entries from CSV")
    except Exception as e:
        print(f"ERROR reading CSV: {e}")
        return
    
    # ============================================================================
    # EDIT 2: Changed 'Antibody' to 'Therapeutic' (column name fix)
    # ============================================================================
    # Validate required columns - FIXED: Use 'Therapeutic' instead of 'Antibody'
    required_cols = ['Therapeutic', 'Format']
    if not all(col in df.columns for col in required_cols):
        print(f"ERROR: CSV missing required columns: {required_cols}")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Get all Python files
    try:
        all_files = os.listdir(py_folder)
        print(f"Total files in folder: {len(all_files)}")
    except Exception as e:
        print(f"ERROR reading folder: {e}")
        return
    
    # Filter for .py files (exclude utility scripts)
    exclude_files = {'__init__.py', 'categorize_antibody_format.py', 'categorize_simple.py', 
                     'therasabdab_analyze_formats.py'}
    py_files = [f for f in all_files if f.endswith('.py') and f not in exclude_files 
                and not os.path.isdir(os.path.join(py_folder, f))]
    
    print(f"Found {len(py_files)} antibody .py files to categorize")
    if py_files:
        print(f"Sample: {py_files[:3]}")
    print()
    
    if len(py_files) == 0:
        print("WARNING: No antibody .py files found!")
        return
    
    # Categorize files
    results = defaultdict(list)
    
    for py_file in py_files:
        name = py_file.replace('.py', '')
        
        # ============================================================================
        # EDIT 3: Changed df['Antibody'] to df['Therapeutic'] (column name fix)
        # ============================================================================
        # Find in CSV (case-insensitive, with whitespace handling) - FIXED: Use 'Therapeutic'
        match = df[df['Therapeutic'].str.lower().str.strip() == name.lower().strip()]
        
        if not match.empty:
            format_type = match.iloc[0]['Format']
            category = get_category(format_type)
            results[category].append((name, format_type, py_file))
        else:
            results['Unclassified'].append((name, 'N/A (not in CSV)', py_file))
    
    # Print categorization results
    print("=" * 80)
    print("CATEGORIZATION RESULTS")
    print("=" * 80)
    
    categories_order = ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats', 'Unclassified']
    
    for category in categories_order:
        if category in results:
            count = len(results[category])
            print(f"\n{category}: {count} files")
            for name, fmt, _ in sorted(results[category])[:3]:
                print(f"  • {name}")
                if fmt != 'N/A (not in CSV)':
                    print(f"    Format: {fmt}")
            if count > 3:
                print(f"  ... and {count - 3} more")
    
    # Confirmation before organizing
    print("\n" + "=" * 80)
    print("READY TO ORGANIZE")
    print("=" * 80)
    print("\nThis will create 4 folders and move files:")
    for cat in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']:
        if cat in results:
            print(f"  • {cat}/ ({len(results[cat])} files)")
    
    response = input("\nContinue? (yes/no): ").strip().lower()
    
    if response != 'yes':
        print("Cancelled.")
        return
    
    # Create folders and move files
    print("\n" + "=" * 80)
    print("ORGANIZING FILES")
    print("=" * 80)
    
    for category in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Other_Formats']:
        category_folder = py_folder / category
        
        if not category_folder.exists():
            category_folder.mkdir(parents=True)
            print(f"Created: {category}/")
        
        if category in results:
            for name, fmt, py_file in results[category]:
                src = py_folder / py_file
                dst = category_folder / py_file
                try:
                    shutil.move(str(src), str(dst))
                except Exception as e:
                    print(f"ERROR moving {py_file}: {e}")
            print(f"✓ Moved {len(results[category])} files to {category}/")
    
    # Handle unclassified files
    if 'Unclassified' in results and results['Unclassified']:
        print(f"\n⚠ {len(results['Unclassified'])} files could not be classified:")
        for name, fmt, _ in results['Unclassified']:
            print(f"  • {name}")
    
    print("\n" + "=" * 80)
    print("Done!")
    print("=" * 80)

if __name__ == "__main__":
    main()