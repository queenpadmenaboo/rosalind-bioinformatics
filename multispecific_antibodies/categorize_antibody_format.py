import pandas as pd
import os
import shutil
from collections import defaultdict

# 4-Category mapping
CATEGORY_MAPPING = {
    'Whole_mAb': ['Whole mAb', 'Whole Ab'],
    'Bispecific_mAb': ['Bispecific mAb', 'Bispecific IgG', 'Bispecific Whole mAb'],
    'Bispecific_scFv': ['scFv', 'BiTE', 'TCE', 'Tandem'],
    'Multispecific_Advanced': ['Trispecific', 'Tetraspecific']
}

def get_category(format_string):
    """Get category from format string"""
    for category, keywords in CATEGORY_MAPPING.items():
        for keyword in keywords:
            if keyword.lower() in format_string.lower():
                return category
    return 'Other'

def main():
    csv_path = r"C:\Users\bunsr\TheraSAbDab_SeqStruc_ 07Dec2025.csv"
    py_folder = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"
    
    print("=" * 80)
    print("CATEGORIZING ANTIBODY FILES")
    print("=" * 80)
    print(f"Looking in folder: {py_folder}")
    print()
    
    # Load CSV
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} entries from CSV")
    
    # Get all files
    try:
        all_files = os.listdir(py_folder)
        print(f"Total files in folder: {len(all_files)}")
    except Exception as e:
        print(f"ERROR reading folder: {e}")
        return
    
    # Filter for .py files
    py_files = []
    for f in all_files:
        if f.endswith('.py'):
            if f not in ['__init__.py', 'categorize_antibody_format.py', 'categorize_simple.py', 'therasabdab_analyze_formats.py']:
                py_files.append(f)
    
    print(f"Found {len(py_files)} antibody .py files")
    print(f"First 5: {py_files[:5]}")
    print()
    
    if len(py_files) == 0:
        print("ERROR: No antibody .py files found!")
        return
    
    # Categorize
    results = defaultdict(list)
    
    for py_file in py_files:
        name = py_file.replace('.py', '')
        
        # Find in CSV
        match = df[df['Therapeutic'].str.lower() == name.lower()]
        
        if not match.empty:
            format_type = match.iloc[0]['Format']
            category = get_category(format_type)
            results[category].append((name, format_type, py_file))
        else:
            results['Not_in_CSV'].append((name, 'N/A', py_file))
    
    # Print results
    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    
    for category in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Multispecific_Advanced', 'Other', 'Not_in_CSV']:
        if category in results:
            count = len(results[category])
            print(f"\n{category}: {count} files")
            for name, fmt, _ in sorted(results[category])[:3]:
                print(f"  - {name}: {fmt}")
            if count > 3:
                print(f"  ... and {count - 3} more")
    
    # Ask before moving
    print("\n" + "=" * 80)
    print("Ready to organize files into folders.")
    print("This will create 4 folders and move files into them.")
    response = input("Continue? (yes/no): ")
    
    if response.lower() != 'yes':
        print("Cancelled.")
        return
    
    # Create folders and move files
    print("\n" + "=" * 80)
    print("ORGANIZING FILES")
    print("=" * 80)
    
    for category in ['Whole_mAb', 'Bispecific_mAb', 'Bispecific_scFv', 'Multispecific_Advanced']:
        category_folder = os.path.join(py_folder, category)
        
        if not os.path.exists(category_folder):
            os.makedirs(category_folder)
            print(f"Created: {category}/")
        
        if category in results:
            for name, fmt, py_file in results[category]:
                src = os.path.join(py_folder, py_file)
                dst = os.path.join(category_folder, py_file)
                try:
                    shutil.move(src, dst)
                except Exception as e:
                    print(f"ERROR moving {py_file}: {e}")
            print(f"Moved {len(results[category])} files to {category}/")
    
    print("\nDone!")

if __name__ == "__main__":
    main()