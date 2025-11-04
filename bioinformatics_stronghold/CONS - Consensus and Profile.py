"""
Rosalind Problem: CONS - Consensus and Profile
Problem URL: http://rosalind.info/problems/cons/

Author: Bunsree Patel
Date: November 4, 2025
"""

"""A set of DNA strings of equal length is considered. Each position across
the strings forms a column. For each column, the occurrences of the four
nucleotides (A, C, G, T) are counted. These counts form a 4 x n profile
matrix, where n is the length of the DNA strings. The rows of the profile
matrix correspond to nucleotide counts in the fixed order: A, C, G, T.

The consensus string is formed by selecting, at each column, the nucleotide
with the highest frequency in that column according to the profile matrix.
If multiple nucleotides share the highest frequency at a position, any of
them may be chosen.

Example:

DNA Strings:
A T C C A G C T
G G G C A A C T
A T G G A T C T
A A G C A A C C
T T G G A A C T
A T G C C A T T
A T G G C A C T

Profile Matrix:
A: 5 1 0 0 5 5 0 0
C: 0 0 1 4 2 0 6 1
G: 1 1 6 3 0 1 0 0
T: 1 5 0 0 0 1 1 6

Consensus String:
A T G C A A C T

Input Format:
A collection of at most 10 DNA strings (each up to 1000 bp) in FASTA format.

Output Format:
1) A consensus string.
2) The profile matrix, printed with rows for A, C, G, and T in that order.
"""

"""Sample Dataset Problem
    Given: A collection of at most 10 DNA strings of equal length (at most 1 kbp) in FASTA format.
    Return: A consensus string and profile matrix for the collection. (If several possible consensus strings exist, then you may return any one of them.)"""

from ProteinToolkit import parse_fasta

# Compute consensus string and profile matrix from sequences.
def consensus_and_profile(fasta_dict):
    sequences = list(fasta_dict.values())
    length = len(sequences[0])
    nucleotides = "ACGT"

    # Initialize profile dictionary
    profile = {nuc: [0]*length for nuc in nucleotides}

    # Fill profile matrix
    for i in range(length):
        for seq in sequences:
            profile[seq[i]][i] += 1

    # Build consensus string
    consensus = ""
    for i in range(length):
        # Find nucleotide with max count in this column
        max_nuc = max(nucleotides, key=lambda n: profile[n][i])
        consensus += max_nuc

    return consensus, profile


fasta_text = """
>Rosalind_1
ATCCAGCT
>Rosalind_2
GGGCAACT
>Rosalind_3
ATGGATCT
>Rosalind_4
AAGCAACC
>Rosalind_5
TTGGAACT
>Rosalind_6
ATGCCATT
>Rosalind_7
ATGGCACT
"""

fasta_dict = parse_fasta(fasta_text)
consensus, profile = consensus_and_profile(fasta_dict)

print(consensus)
for nuc in "ACGT":
    print(f"{nuc}: {''.join(map(str, profile[nuc]))}")

