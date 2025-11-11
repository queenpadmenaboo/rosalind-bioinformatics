"""
Rosalind Problem: SPLC - RNA Splicing
Problem URL: http://rosalind.info/problems/splc/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 11, 2025
"""

"""
'RNA polymerase' (RNAP) - breaks the bonds joining complementary bases of DNA
'complementary bases'
'precursor mRNA' (pre-mRNA)
'template strand'
'coding strand'
replacement of Thymine with Uracil by RNAP
'double helix'

'introns'
'exons'
'splicing'
'spliceosome'
exons from a gene are known as gene's 'coding region'
"""

"""Sample Dataset Problem
After identifying the exons and introns of an RNA string, 
we only need to delete the introns and concatenate the exons to form a new string ready for translation.

    Given: A DNA string s (of length at most 1 kbp) and a collection of substrings of s acting as introns. All strings are given in FASTA format.
    Return: A protein string resulting from transcribing and translating the exons of s. 
    (Note: Only one solution will exist for the dataset provided.)
"""

from DNAToolkit import solve_splc

fasta_text = """
>Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT
"""

# In FASTA format, the first sequence is always the main gene, and all subsequent sequences are the introns to remove.
# In real genes: You start with the full gene (exons + introns mixed together)
