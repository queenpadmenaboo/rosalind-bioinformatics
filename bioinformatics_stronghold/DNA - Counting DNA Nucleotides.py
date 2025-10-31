"""
Rosalind Problem: DNA - Counting DNA Nucleotides
Problem URL: http://rosalind.info/problems/dna/

Author: Bunsree Patel
Date: October 24, 2025
"""

"""DNA (deoxyribonucleic acid) is a molecule that carries the genetic instructions used in the growth, development, functioning, and reproduction of all known living organisms and many viruses. 
DNA is composed of four types of nucleotides, adenine (A), cytosine (C), guanine (G), and thymine (T). In RNA, uracil (U) replaces thymine.

The nucleic acid monomer is called a nucleotide (abbreaviated to nt), and is used as a unit of strand length. 
Each nucleotide consists of a sugar molecule, a negatively charged phosphate group, and a nitrogenous base (nucleobase)
A sugar-phosphate backbone is formed by covalent bonds between the sugar of one nucleotide and the phosphate of the next nucleotide in the chain.
Nucleotides of a specific type of nucleic acid (DNA or RNA) differ only in the nitrogenous base they contain.

The ordering of bases defines a nucleic acid's primary structure.

DNA is found in all living organisms on Earth, including bacteria, and many viruses (which are often considered to be nonliving).
Genome refers to the sum total of the DNA in an organism's chromosomes."""


"""Sample Problem:
A string is an ordered collection of symbols selected from some alphabet and formed into a word; the length of a string is the number of symbols that it contains.
An example of a length 21 DNA string (whose alphabet contains the symbols 'A', 'C', 'G', and 'T') is "ATGCTTCAGAAAGGTCTTACG."

Given: A DNA string s of length at most 1000 nt.
Return: Four integers (separated by spaces) counting the respective number of times that the symbols 'A', 'C', 'G', and 'T' occur in s."""

"""Sample Dataset Problem"""

from DNAToolkit import validateSeq, countDNANucFrequency
dna_seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"

print(validateSeq(dna_seq))
print(countDNANucFrequency(dna_seq))
"""Sample Output:
20 12 17 21
Explanation:
In the sample dataset, there are 20 occurrences of 'A', 12 occurrences of 'C', 17 occurrences of 'G', and 21 occurrences of 'T'."""

"""Actual Dataset Problem"""

from DNAToolkit import validateSeq, countDNANucFrequency
dna_seq = "CGCCAGAGCATCGGATACCTAGGATAAAACGAAGAGGGCATCCTCTCCATGATTCGGACAGGATAACTCATGCGCAAAATTTTATGTAACCGCGAGGGCACCGGGGACACAATTTACGTCTGAAGCTAAGCAACTAACGAAGTTGGATCAGAACTACTTTCTCATTATCGCTCACAATACATAAGCTTCCATCATCAACTTGAGGCATAATCGGATATTCGTGGGGGGGTATACCGTCGGTACATCAGCCCCTCTAATTAGTCCGTCTCCATATGTTAGACCTTTCTGGAACGTGACTGGGCTGTTATGAAGATGATTTCACTAACTTAGCCCTAGGTCTACAGAATGTGGCTAAGCCAATGTTACAATAATACAGGTGCTCGTCGTGATTCCCCACTTCTTGATAATGGATTCGTGTTACATGAGATCGTGGGCCGCGTCTCACACAGCGATGGCGTCTACGTGAGCTAATGCTGACATATCACGCTATAGGTAAGGATGGTACGGTAGGCTCGTATGCAGCTGCCTATTAGTGGCAGGATCGGCCATTTACCTCCAGCACATGGATCGCGGTCTGGCAGGTCTTCGTTATATCACGTCGCATTGGCCACTTTGCGCTCAGTTGGTATCGCCGGAAGTCAGCGTTACTCGTGTAAAAAGTCGCCTCATCGAGAGGAACCTCAGCTCTTATCCTCGCGAGGGCAAACGCGTTCACGCTATAAACAGCCCGGCTCGGTTGTACATACGTCGCCGTTTATAGCTGAGAGCCTGTCTCTAGACAGACAAGTAACATCCGATCTAGCTATCGATCCAGCATCACGGATACTCCGGTAAACTATTATACGCTAGCCCTCCTGTGCTTGTTTTAACATAAAACCCG"

print(validateSeq(dna_seq))
print(countDNANucFrequency(dna_seq))
# Expected Output: {'A': 225, 'C': 221, 'G': 206, 'T': 228}
# Expected Output: "225 221 206 228"