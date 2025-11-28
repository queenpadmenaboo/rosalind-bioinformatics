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

def nucleotide_count(dna_seq):
    