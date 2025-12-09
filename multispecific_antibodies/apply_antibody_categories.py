"""
Add 4-category classification to antibody_metadata.py
Reads existing metadata and adds 'ml_category' field
"""

import os
import sys

# Import the category mapping
from antibody_categories import CATEGORY_MAPPING, get_category_from_format, get_category_statistics


def add_categories_to_metadata(metadata_file_path, output_file_path=None):
    """
    Read antibody_metadata.py and add ml_category field to each antibody
    
    Args:
        metadata_file_path: Path to antibody_metadata.py
        output_file_path: Path to save updated file (if None, overwrites original)
    """
    
    if output_file_path is None:
        output_file_path = metadata_file_path
    
    # Import the metadata
    metadata_dir = os.path.dirname(metadata_file_path)
    sys.path.insert(0, metadata_dir)
    
    try:
        import antibody_metadata
        ANTIBODY_FORMATS = antibody_metadata.ANTIBODY_FORMATS
        FORMAT_GROUPS = antibody_metadata.FORMAT_GROUPS
    except ImportError:
        print(f"ERROR: Could not import antibody_metadata from {metadata_file_path}")
        return
    
    print("=" * 80)
    print("ADDING ML CATEGORIES TO ANTIBODY METADATA")
    print("=" * 80)
    print(f"Input file: {metadata_file_path}")
    print(f"Output file: {output_file_path}")
    print()
    
    # Add ml_category to each antibody
    updated_count = 0
    unmapped_formats = set()
    
    for antibody_name, metadata in ANTIBODY_FORMATS.items():
        format_name = metadata['format']
        category = get_category_from_format(format_name)
        
        if category:
            metadata['ml_category'] = category
            updated_count += 1
        else:
            metadata['ml_category'] = 'UNMAPPED'
            unmapped_formats.add(format_name)
    
    print(f"Updated {updated_count} antibodies with ML categories")
    
    if unmapped_formats:
        print(f"\nWARNING: {len(unmapped_formats)} format(s) not mapped to any category:")
        for fmt in sorted(unmapped_formats):
            print(f"  - {fmt}")
    
    # Get category statistics
    category_stats = {}
    for antibody_name, metadata in ANTIBODY_FORMATS.items():
        cat = metadata['ml_category']
        category_stats[cat] = category_stats.get(cat, 0) + 1
    
    print()
    print("=" * 80)
    print("CATEGORY DISTRIBUTION")
    print("=" * 80)
    for category in sorted(category_stats.keys()):
        count = category_stats[category]
        print(f"{category}: {count} antibodies")
    
    # Create reverse index by ml_category
    CATEGORY_GROUPS = {}
    for category in CATEGORY_MAPPING.keys():
        CATEGORY_GROUPS[category] = []
    CATEGORY_GROUPS['UNMAPPED'] = []
    
    for antibody_name, metadata in ANTIBODY_FORMATS.items():
        cat = metadata['ml_category']
        CATEGORY_GROUPS[cat].append(antibody_name)
    
    # Write updated metadata file
    print()
    print("=" * 80)
    print("WRITING UPDATED METADATA FILE")
    print("=" * 80)
    
    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write('"""\n')
        f.write('Antibody Metadata with ML Categories\n')
        f.write('Auto-generated from TheraSAbDab CSV\n')
        f.write('Maps antibody name -> format, attributes, and ML category\n')
        f.write('"""\n\n')
        
        # Write ANTIBODY_FORMATS with ml_category
        f.write('ANTIBODY_FORMATS = {\n')
        for ab_name in sorted(ANTIBODY_FORMATS.keys()):
            metadata = ANTIBODY_FORMATS[ab_name]
            f.write(f"    '{ab_name}': {{\n")
            f.write(f"        'format': '{metadata['format']}',\n")
            f.write(f"        'chains': {metadata['chains']},\n")
            f.write(f"        'is_multispecific': {metadata['is_multispecific']},\n")
            f.write(f"        'clinical_status': '{metadata['clinical_status']}',\n")
            f.write(f"        'target': '{metadata['target']}',\n")
            f.write(f"        'companies': '{metadata['companies']}',\n")
            f.write(f"        'ml_category': '{metadata['ml_category']}'\n")
            f.write(f"    }},\n")
        f.write('}\n\n')
        
        # Write FORMAT_GROUPS
        f.write('# Reverse index: format -> list of antibodies\n')
        f.write('FORMAT_GROUPS = {\n')
        for format_name in sorted(FORMAT_GROUPS.keys()):
            f.write(f"    '{format_name}': [\n")
            for ab_name in sorted(FORMAT_GROUPS[format_name]):
                f.write(f"        '{ab_name}',\n")
            f.write(f"    ],\n")
        f.write('}\n\n')
        
        # Write CATEGORY_GROUPS
        f.write('# ML category index: category -> list of antibodies\n')
        f.write('CATEGORY_GROUPS = {\n')
        for category in sorted(CATEGORY_GROUPS.keys()):
            f.write(f"    '{category}': [\n")
            for ab_name in sorted(CATEGORY_GROUPS[category]):
                f.write(f"        '{ab_name}',\n")
            f.write(f"    ],\n")
        f.write('}\n\n')
        
        # Write helper functions
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
        
        f.write('def get_ml_category(antibody_name):\n')
        f.write('    """Get ML category for antibody"""\n')
        f.write('    if antibody_name not in ANTIBODY_FORMATS:\n')
        f.write('        return None\n')
        f.write('    return ANTIBODY_FORMATS[antibody_name]["ml_category"]\n\n')
        
        f.write('def get_antibodies_by_format(format_name):\n')
        f.write('    """Get all antibodies of a specific format"""\n')
        f.write('    return FORMAT_GROUPS.get(format_name, [])\n\n')
        
        f.write('def get_antibodies_by_category(category):\n')
        f.write('    """Get all antibodies in a specific ML category"""\n')
        f.write('    return CATEGORY_GROUPS.get(category, [])\n\n')
        
        f.write('def get_category_counts():\n')
        f.write('    """Get count of antibodies in each category"""\n')
        f.write('    return {cat: len(abs) for cat, abs in CATEGORY_GROUPS.items()}\n')
    
    print(f"Metadata file updated: {output_file_path}")
    print()
    print("=" * 80)
    print("COMPLETE")
    print("=" * 80)
    print()
    print("You can now use:")
    print("  - get_ml_category('amivantamab') -> 'Bispecific_mAb'")
    print("  - get_antibodies_by_category('Bispecific_scFv') -> list of antibodies")
    print("  - get_category_counts() -> {'Whole_mAb': 308, ...}")


if __name__ == "__main__":
    # Update the metadata file
    antibody_folder = r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies"
    metadata_file = os.path.join(antibody_folder, "antibody_metadata.py")
    
    if not os.path.exists(metadata_file):
        print(f"ERROR: antibody_metadata.py not found at {metadata_file}")
        print("Please run therasabdab_analyze_formats.py first to generate metadata file")
        sys.exit(1)
    
    add_categories_to_metadata(metadata_file)