"""
Rosalind Problem: SUBS - Finding a Motif in DNA
Problem URL: http://rosalind.info/problems/subs/

Author: Bunsree Patel
Date: October 30, 2025
"""

"""Finding the same interval of DNA in the genomes of two different organisms (often taken from different species) is highly suggestive that the interval has the same function in both organisms.
A motif is defined as a commonly shared interval of DNA.
A common task in molecular biology is to search an organism's genome for a known motif.

Genomes have intervals of DNA that occur multiple times (possibly with slight modifications), called repeats.
These repeats occur far more often than would be dictated by random chance, indicating that genomes
are anything but random --> shows that language of DNA must be very powerful.

The most common repeat in humans is the Alu repeat, which is approximately 300 bp long and recurs around a million times throughout every human genome.
Alu, however, has not been found to serve a positive purpose, and appears in fact to be parasitic: when a new Alu repeat is inserted into a genome, it frequently causes genetic disorders."""

"""Sample Dataset Problem
Given two strings s and t, t is a substring of s if t is contained as a contiguous collection of symbols in s (as a result, t must be no longer than s).

The position of a symbol in a string is the total number of symbols found to its left, including itself (e.g., the positions of all occurrences of 'U' in "AUGCUUCAGAAAGGUCUUACG" are 2, 5, 6, 15, 17, and 18). The symbol at position i of s is denoted by s[i].

A substring of s can be represented as s[j:k], where j and k represent the starting and ending positions of the substring in s; for example, if s = "AUGCUUCAGAAAGGUCUUACG", then s[2:5] = "UGCU".

The location of a substring s[j:k] is its beginning position j; note that t will have multiple locations in s if it occurs more than once as a substring of s (see the Sample below).

    Given: Two DNA strings s and t (each of length at most 1 kbp).
    Return: All locations of t as a substring of s.
"""

from DNAToolkit import find_substring_locations

# Define strings
s = "GATATATGCATATACTT"
t = "ATAT"

# Call the function
result = find_substring_locations(s,t)

# Format and print output
print(*result)
"""Output: 2 4 10"""