# Multispecific Antibodies Sequence Database

A collection of therapeutic multispecific and monospecific antibody sequences in FASTA format.

## Data Source
Antibody sequences obtained from Thera-SAbDab (http://opig.stats.ox.ac.uk/webapps/therasabdab/search/)

Citation: Raybould et al. (2020) Thera-SAbDab: the Therapeutic Structural Antibody Database. Nucleic Acids Research.

## Repository Contents
- Individual `.py` files containing FASTA sequences for each therapeutic antibody
- `sabdabconverter.py` - Tool for converting Thera-SAbDab data to FASTA format
- `readme_count.py` - Script to update antibody count in README

## Usage
Each antibody file contains sequences in FASTA format stored as a string variable:
```python
from odronextamab import odronextamab
print(odronextamab)
```

Or parse the FASTA:
```python
with open('odronextamab.py', 'r') as f:
    content = f.read()
    # Extract sequences from triple-quoted string
```

## License
- **Code** (sabdabconverter.py, readme_count.py): MIT License
- **Data** (antibody sequences): Sourced from Thera-SAbDab (publicly available data)

## Disclaimer
These sequences are from publicly available therapeutic antibodies as curated by Thera-SAbDab. For commercial use, please verify current patent status and regulatory information.

## Statistics
- Data last updated: December 2025
- Total antibodies: 62