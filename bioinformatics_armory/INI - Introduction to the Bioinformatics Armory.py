"""
Rosalind Problem: INI - Introduction to the Bioinformatics Armory
Problem URL: http://rosalind.info/problems/ini/
Location: Bioinformatics Armory
Author: Bunsree Patel
Date: November 27, 2025
"""

"""Sample Dataset Problem

This initial problem is aimed at familiarizing you with Rosalind's task-solving pipeline.
To solve it, you merely have to take a given DNA sequence and find its nucleotide counts; this problem is equivalent to “Counting DNA Nucleotides” in the Stronghold.

Of the many tools for DNA sequence analysis, one of the most popular is the Sequence Manipulation Suite.
Commonly known as SMS 2, it comprises a collection of programs for generating, formatting, and analyzing short strands of DNA and polypeptides.

One of the simplest SMS 2 programs, called 'DNA stats', counts the number of occurrences of each nucleotide in a given strand of DNA.
An online interface for DNA stats can be found at https://www.bioinformatics.org/sms2/dna_stats.html

    Given: A DNA string s of length at most 1000 bp.

    Return: Four integers (separated by spaces) representing the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s.
    
    Note: You must provide your answer in the format shown in the sample output below."""


dna_seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
def nucleotide_count(dna_seq):
    TmpFreqDict = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for nuc in dna_seq:
        if nuc in TmpFreqDict:
            TmpFreqDict[nuc] += 1

    return TmpFreqDict

print(nucleotide_count(dna_seq))
"""Output:
{'A': 20, 'C': 12, 'G': 17, 'T': 21}"""

# Abbreviated Method without TmpFreqDict
def nuc_count(dna_seq):
    print(
        "A:", dna_seq.count('A'),
        "C:", dna_seq.count('C'),
        "G:", dna_seq.count('G'),
        "T:", dna_seq.count('T')
    )

dna_seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
nuc_count(dna_seq)
"""Output:
A: 20 C: 12 G: 17 T: 21"""

"""Actual Dataset Problem"""

def nuc_counttwo(dna_seq):
    print(
        "A:", dna_seq.count('A'),
        "C:", dna_seq.count('C'),
        "G:", dna_seq.count('G'),
        "T:", dna_seq.count('T')
    )
dna_seq ="GGTCAGATCGTGGGGGCAATTGTACCGCGACTACTAACCCACGCCGAAACCCATCGCGTGTACGGATTACGCCGGCGAGGTGGAATACGCAAATTGATATGTACTAAGTAAATATGCGCTCAACATTTCTTGGCCTGTCTTAAGTTACTTGATGACATAATCCGCAATCGGAATTCTATGAAACGTCGAACTGGGTATTTCGCAGTTGTAATGTTATGAAGTCACCCAGAAACGAGTTTTAATGGCGGCCTGGACGTATCGTAGCGTGCAGAGGAGGTGGGCAGCTCTAGTACCCGAGCTACCCGTCGTATAGTCGTGAATGTCGGTCATACTGACGTCAAGAGCTCGATAGATGTAGGTACTGGTGTAGTTCCTTTGCCATAGCATCTTGGAGGAGCGAATGTTACGATTGAAGGTCATGCCCGCAATGTCCATAATACTATTCCGGCTAGGCCATTGATGAATTCCGGATCGTCGAAGCGTGAAATTGTACGGGACAGATATTGTGTGTATATCTACTGTCCAGAAACATTTATGTCGGTATGCGGGACGCAACTAAGTAACATACCGTTAGTCCTCGCAAGCTGGATCCTAGGATTCCCAAGTATGAACGAACCTCTTCAGCTCTCACTTCAATAAATTCCAGGTCGCCCGAACCACTGTTCGCTACAAGTATAGGGCCACGAGACATGTATTAATCAAGATTATAAACCTTTCACCTTCTTATACAGACCCTCATTCGTTTGCAATCCGCAATTATCGGTGAAAATTGACTCGTGAGCGCTTATCCACTAGGGTGACGGCCAATAGTGGCCGTACCTAACGACCCTGGCATATTCTTCTTTATCCATATTAG"
nuc_counttwo(dna_seq)

"""Output:
A: 227 C: 197 G: 201 T: 231"""

from Bio.Seq import Seq

my_seq = Seq("AGTACACTGGGGTTTTCCCAACTACTTGGAGA")
print(
    my_seq.count('A'),
    my_seq.count('C'),
    my_seq.count('G'),
    my_seq.count('T')
)
"""Output:
8 7 8 9"""
