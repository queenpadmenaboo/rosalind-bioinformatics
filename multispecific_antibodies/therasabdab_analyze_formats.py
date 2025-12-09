import pandas as pd
import os
from collections import defaultdict, Counter

def analyze_antibody_formats(csv_path, antibody_folder):
    """
    Analyze TheraSAbDab CSV to extract format categories and counts.
    Suggests optimal grouping for machine learning.
    """
    
    print("=" * 80)
    print("ANTIBODY FORMAT ANALYSIS")
    print("=" * 80)
    
    # Load CSV
    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} therapeutics from CSV\n")
    
    # Get all unique formats
    format_counts = df['Format'].value_counts()
    print("=" * 80)
    print("ALL FORMATS IN DATABASE")
    print("=" * 80)
    for format_name, count in format_counts.items():
        print(f"  {format_name}: {count} antibodies")
    
    print(f"\nTotal unique formats: {len(format_counts)}")
    
    # Identify multispecific formats
    print("\n" + "=" * 80)
    print("MULTISPECIFIC vs MONOSPECIFIC BREAKDOWN")
    print("=" * 80)
    
    multispecific_keywords = ['Bispecific', 'Trispecific', 'Tetraspecific', 'bispecific', 'trispecific']
    
    multispecific_formats = []
    monospecific_formats = []
    
    for format_name in format_counts.index:
        is_multi = any(keyword in format_name for keyword in multispecific_keywords)
        if is_multi:
            multispecific_formats.append(format_name)
        else:
            monospecific_formats.append(format_name)
    
    print(f"\nMULTISPECIFIC FORMATS ({len(multispecific_formats)} types):")
    multispecific_total = 0
    for format_name in multispecific_formats:
        count = format_counts[format_name]
        multispecific_total += count
        print(f"  {format_name}: {count}")
    
    print(f"\nMONOSPECIFIC/OTHER FORMATS ({len(monospecific_formats)} types):")
    monospecific_total = 0
    for format_name in monospecific_formats:
        count = format_counts[format_name]
        monospecific_total += count
        print(f"  {format_name}: {count}")
    
    print(f"\nTOTAL MULTISPECIFIC: {multispecific_total}")
    print(f"TOTAL MONOSPECIFIC/OTHER: {monospecific_total}")
    
    # Check which antibodies exist in your .py files
    print("\n" + "=" * 80)
    print("CHECKING YOUR ANTIBODY .PY FILES")
    print("=" * 80)
    
    py_files = [f.replace('.py', '') for f in os.listdir(antibody_folder) 
                if f.endswith('.py') and f != '__init__.py']
    print(f"Found {len(py_files)} .py files\n")
    
    # Match py files to CSV entries
    df['therapeutic_lower'] = df['Therapeutic'].str.lower()
    
    format_distribution_in_files = defaultdict(int)
    antibodies_by_format = defaultdict(list)
    
    for py_name in py_files:
        match = df[df['therapeutic_lower'] == py_name.lower()]
        if not match.empty:
            format_name = match.iloc[0]['Format']
            format_distribution_in_files[format_name] += 1
            antibodies_by_format[format_name].append(py_name)
    
    print("FORMAT DISTRIBUTION IN YOUR .PY FILES:")
    for format_name, count in sorted(format_distribution_in_files.items(), 
                                     key=lambda x: x[1], reverse=True):
        print(f"  {format_name}: {count} files")
    
    # Suggest ML categories - using EXACT TheraSAbDab format names
    print("\n" + "=" * 80)
    print("SUGGESTED ML CATEGORIES (Using exact TheraSAbDab format names)")
    print("=" * 80)
    
    print("\nFor a multispecific predictor, suggested grouping:\n")
    
    # Group 1: Whole mAb (using exact CSV term)
    whole_mab = [f for f in format_distribution_in_files.keys() 
                 if 'Whole mAb' in f or f == 'Whole mAb']
    print("GROUP 1: Whole mAb")
    for fmt in whole_mab:
        print(f"  - {fmt}: {format_distribution_in_files[fmt]} files")
    
    # Group 2: Bispecific mAb
    bispecific_mab = [f for f in format_distribution_in_files.keys() 
                      if 'Bispecific mAb' in f or 'Bispecific IgG' in f]
    print("\nGROUP 2: Bispecific mAb")
    for fmt in bispecific_mab:
        print(f"  - {fmt}: {format_distribution_in_files[fmt]} files")
    
    # Group 3: Bispecific scFv
    bispecific_scfv = [f for f in format_distribution_in_files.keys() 
                       if 'scFv' in f or 'TCE' in f or 'tandem' in f.lower()]
    print("\nGROUP 3: Bispecific scFv")
    for fmt in bispecific_scfv:
        print(f"  - {fmt}: {format_distribution_in_files[fmt]} files")
    
    # Group 4: Trispecific
    trispecific = [f for f in format_distribution_in_files.keys() 
                   if 'Trispecific' in f]
    print("\nGROUP 4: Trispecific")
    for fmt in trispecific:
        print(f"  - {fmt}: {format_distribution_in_files[fmt]} files")
    
    # Group 5: Other formats
    categorized = set(whole_mab + bispecific_mab + bispecific_scfv + trispecific)
    other = [f for f in format_distribution_in_files.keys() if f not in categorized]
    print("\nGROUP 5: OTHER FORMATS")
    for fmt in other:
        print(f"  - {fmt}: {format_distribution_in_files[fmt]} files")
    
    # Generate Python dictionary for metadata
    print("\n" + "=" * 80)
    print("GENERATING ANTIBODY METADATA DICTIONARY")
    print("=" * 80)
    
    metadata_dict = {}
    
    for py_name in py_files:
        match = df[df['therapeutic_lower'] == py_name.lower()]
        if not match.empty:
            row = match.iloc[0]
            
            # Determine if multispecific
            is_multi = any(keyword in str(row['Format']) for keyword in multispecific_keywords)
            
            # Count chains (rough estimate based on format)
            if 'Bispecific' in str(row['Format']):
                chains = 4
            elif 'Trispecific' in str(row['Format']):
                chains = 6
            elif 'Fab' in str(row['Format']):
                chains = 2
            elif 'Whole mAb' in str(row['Format']):
                chains = 2
            else:
                chains = 'Unknown'
            
            metadata_dict[py_name] = {
                'format': str(row['Format']),
                'chains': chains,
                'is_multispecific': is_multi,
                'clinical_status': str(row['Est. Status']) if pd.notna(row['Est. Status']) else 'Unknown',
                'target': str(row['Target']) if pd.notna(row['Target']) else 'Unknown',
                'companies': str(row['Companies']) if pd.notna(row['Companies']) else 'Unknown'
            }
    
    print(f"\nCreated metadata for {len(metadata_dict)} antibodies")
    
    # Save to Python file
    output_file = os.path.join(antibody_folder, 'antibody_metadata.py')
    
    with open(output_file, 'w') as f:
        f.write('"""\n')
        f.write('Antibody Metadata\n')
        f.write('Auto-generated from TheraSAbDab CSV\n')
        f.write('Maps antibody name -> format and attributes\n')
        f.write('"""\n\n')
        
        f.write('ANTIBODY_FORMATS = {\n')
        for ab_name in sorted(metadata_dict.keys()):
            metadata = metadata_dict[ab_name]
            f.write(f"    '{ab_name}': {{\n")
            f.write(f"        'format': '{metadata['format']}',\n")
            f.write(f"        'chains': {metadata['chains']},\n")
            f.write(f"        'is_multispecific': {metadata['is_multispecific']},\n")
            f.write(f"        'clinical_status': '{metadata['clinical_status']}',\n")
            f.write(f"        'target': '{metadata['target']}',\n")
            f.write(f"        'companies': '{metadata['companies']}'\n")
            f.write(f"    }},\n")
        f.write('}\n\n')
        
        # Create reverse index by format
        f.write('# Reverse index: format -> list of antibodies\n')
        f.write('FORMAT_GROUPS = {\n')
        for format_name in sorted(antibodies_by_format.keys()):
            f.write(f"    '{format_name}': [\n")
            for ab_name in sorted(antibodies_by_format[format_name]):
                f.write(f"        '{ab_name}',\n")
            f.write(f"    ],\n")
        f.write('}\n\n')
        
        # Add helper functions
        f.write('def is_multispecific(antibody_name):\n')
        f.write('    """Check if antibody is multispecific"""\n')
        f.write('    if antibody_name not in ANTIBODY_FORMATS:\n')
        f.write('        return None\n')
        f.write('    return ANTIBODY_FORMATS[antibody_name]["is_multispecific"]\n\n')
        
        f.write('def get_format(antibody_name):\n')
        f.write('    """Get format type for antibody"""\n')
        f.write('    if antibody_name not in ANTIBODY_FORMATS:\n')
        f.write('        return None\n')
        f.write('    return ANTIBODY_FORMATS[antibody_name]["format"]\n\n')
        
        f.write('def get_antibodies_by_format(format_name):\n')
        f.write('    """Get all antibodies of a specific format"""\n')
        f.write('    return FORMAT_GROUPS.get(format_name, [])\n')
    
    print(f"\nMetadata file saved to: {output_file}")
    
    # Final recommendations
    print("\n" + "=" * 80)
    print("RECOMMENDATIONS FOR ML MODEL")
    print("=" * 80)
    
    print("\nBased on your data distribution, here are recommendations:\n")
    
    print("1. PRIMARY CATEGORIES (for 'format' feature):")
    print("   - Use 5 categories as shown in GROUP 1-5 above")
    print("   - This balances granularity with sample size\n")
    
    print("2. BINARY CLASSIFICATION option:")
    print("   - Multispecific (all bi/tri formats) vs Monospecific")
    print("   - Simpler model, more samples per class\n")
    
    print("3. SAMPLE SIZE CONSIDERATIONS:")
    total_multi = sum(format_distribution_in_files[f] for f in format_distribution_in_files 
                      if any(kw in f for kw in multispecific_keywords))
    print(f"   - Total multispecifics in files: {total_multi}")
    print(f"   - For ML, aim for 30+ samples per category")
    print(f"   - Consider grouping rare formats into 'Other'\n")
    
    print("4. FEATURE ENGINEERING:")
    print("   - Bispecific mAbs: Focus on Fc engineering, chain pairing")
    print("   - scFv formats: Linker length/composition, aggregation risk")
    print("   - All formats: Standard developability (pI, hydrophobicity, etc)\n")
    
    return metadata_dict, format_counts, antibodies_by_format


if __name__ == "__main__":
    csv_path = r"C:\Users\bunsr\TheraSAbDab_SeqStruc_ 07Dec2025.csv"
    antibody_folder = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"
    
    analyze_antibody_formats(csv_path, antibody_folder)