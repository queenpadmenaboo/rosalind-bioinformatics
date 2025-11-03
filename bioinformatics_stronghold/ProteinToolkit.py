# Protein Toolkit file
"""
Author: Bunsree Patel
Date: November 3, 2025
"""

def clean_fasta(fasta_text):
    lines = fasta_text.strip().split("\n")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence

# Remove FASTA headers and join sequence lines
# Input: fasta_text (string)
# Output: cleaned sequence (string)
