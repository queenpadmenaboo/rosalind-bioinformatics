"""
Rosalind Problem: PERM - Enumerating Gene Orders
Problem URL: http://rosalind.info/problems/perm/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 8, 2025
"""

"""
Genome Rearrangement Concepts:

'Point mutations' are small-scale changes to individual DNA nucleotides, like changing an 'A' to a 'T'. 
    These single letter changes in the genetic code don't affect the order or structure of genes on chromosomes.
'Genome rearrangements' are large-scale structural changes where entire chunks of DNA get moved, flipped, duplicated, or deleted. 
    Unlike point mutations, these rearrangements change the order and position of genes along chromosomes.
'Synteny blocks' are regions of DNA that contain the same genes in the same order across different species. 
    Comparing genomes between species reveals these conserved blocks, which help understand evolutionary relationships. Even though the exact DNA sequence within these blocks may differ slightly due to point mutations, the blocks are considered equivalent because the same genes appear in the same order.
'Chromosomes' are the long DNA molecules that store genetic information. 
    Different species have different numbers of chromosomes, and over evolutionary time, these chromosomes can undergo rearrangements that change how genes are organized.
Studying genome rearrangements involves labeling each synteny block with a positive integer. 
    This approach ignores the small differences within blocks and focuses only on the order of blocks. 
    When comparing genomes, only the order of numbered synteny blocks matters.
    For example, if one species has blocks in order [1, 2, 3, 4] and another has [1, 3, 2, 4], blocks 2 and 3 have been swapped, representing a rearrangement event."""


"""Sample Dataset Problem
A permutation of length n is an ordering of the positive integers {1,2,...,n}.
For example π = (5,3,2,1,4) is a permutation of length 5.

    Given: A positive integer ≤ 7.
    Return: The total number of permutations of length n, followed by a list of all such permutations (in any order)."""

def permutation(arr):
    # Base case: single element has only one permutations
    if len(arr) == 1:
        return [arr]
    
    result = []
    # For each element in an array
    for i in range(len(arr)):
        # Take element at index i
        current = arr[i]
        # Get remaining elements (exclude current)
        """ arr[:i] = slices everything before index i (does not include i)
            arr[i+1:] = slices everything after index i (from i+1 to end)
            code below combines both parts = everything except index i"""
        remaining = arr [:i] + arr[i+1:]
        # Recursively permute remaining elements
        for p in permutation(remaining):
            # Add current element to front of each permutation
            result.append([current] + p)

    return result

# Input Sample Data
n = 3

# Generate all permutations of [1, 2, ..., n]
perms = permutation(list(range(1, n+1)))

# Output: total count
print(len(perms))

# Output: each permutation as space-seperated integers
for p in perms:
    print(' '.join(map(str, p)))
"""Output: 
6
1 2 3
1 3 2
2 1 3
2 3 1
3 1 2
3 2 1
"""

"""Actual Dataset Problem"""
n = 6
perms = permutation(list(range(1, n+1)))
print(len(perms))
for p in perms:
    print(' '.join(map(str, p)))
"""Output:
720
1 2 3 4 5 6
1 2 3 4 6 5
1 2 3 5 4 6
1 2 3 5 6 4
1 2 3 6 4 5
1 2 3 6 5 4
1 2 4 3 5 6
1 2 4 3 6 5
1 2 4 5 3 6
1 2 4 5 6 3
1 2 4 6 3 5
1 2 4 6 5 3
1 2 5 3 4 6
1 2 5 3 6 4
1 2 5 4 3 6
1 2 5 4 6 3
1 2 5 6 3 4
1 2 5 6 4 3
1 2 6 3 4 5
1 2 6 3 5 4
1 2 6 4 3 5
1 2 6 4 5 3
1 2 6 5 3 4
1 2 6 5 4 3
1 3 2 4 5 6
1 3 2 4 6 5
1 3 2 5 4 6
1 3 2 5 6 4
1 3 2 6 4 5
1 3 2 6 5 4
1 3 4 2 5 6
1 3 4 2 6 5
1 3 4 5 2 6
1 3 4 5 6 2
1 3 4 6 2 5
1 3 4 6 5 2
1 3 5 2 4 6
1 3 5 2 6 4
1 3 5 4 2 6
1 3 5 4 6 2
1 3 5 6 2 4
1 3 5 6 4 2
1 3 6 2 4 5
1 3 6 2 5 4
1 3 6 4 2 5
1 3 6 4 5 2
1 3 6 5 2 4
1 3 6 5 4 2
1 4 2 3 5 6
1 4 2 3 6 5
1 4 2 5 3 6
1 4 2 5 6 3
1 4 2 6 3 5
1 4 2 6 5 3
1 4 3 2 5 6
1 4 3 2 6 5
1 4 3 5 2 6
1 4 3 5 6 2
1 4 3 6 2 5
1 4 3 6 5 2
1 4 5 2 3 6
1 4 5 2 6 3
1 4 5 3 2 6
1 4 5 3 6 2
1 4 5 6 2 3
1 4 5 6 3 2
1 4 6 2 3 5
1 4 6 2 5 3
1 4 6 3 2 5
1 4 6 3 5 2
1 4 6 5 2 3
1 4 6 5 3 2
1 5 2 3 4 6
1 5 2 3 6 4
1 5 2 4 3 6
1 5 2 4 6 3
1 5 2 6 3 4
1 5 2 6 4 3
1 5 3 2 4 6
1 5 3 2 6 4
1 5 3 4 2 6
1 5 3 4 6 2
1 5 3 6 2 4
1 5 3 6 4 2
1 5 4 2 3 6
1 5 4 2 6 3
1 5 4 3 2 6
1 5 4 3 6 2
1 5 4 6 2 3
1 5 4 6 3 2
1 5 6 2 3 4
1 5 6 2 4 3
1 5 6 3 2 4
1 5 6 3 4 2
1 5 6 4 2 3
1 5 6 4 3 2
1 6 2 3 4 5
1 6 2 3 5 4
1 6 2 4 3 5
1 6 2 4 5 3
1 6 2 5 3 4
1 6 2 5 4 3
1 6 3 2 4 5
1 6 3 2 5 4
1 6 3 4 2 5
1 6 3 4 5 2
1 6 3 5 2 4
1 6 3 5 4 2
1 6 4 2 3 5
1 6 4 2 5 3
1 6 4 3 2 5
1 6 4 3 5 2
1 6 4 5 2 3
1 6 4 5 3 2
1 6 5 2 3 4
1 6 5 2 4 3
1 6 5 3 2 4
1 6 5 3 4 2
1 6 5 4 2 3
1 6 5 4 3 2
2 1 3 4 5 6
2 1 3 4 6 5
2 1 3 5 4 6
2 1 3 5 6 4
2 1 3 6 4 5
2 1 3 6 5 4
2 1 4 3 5 6
2 1 4 3 6 5
2 1 4 5 3 6
2 1 4 5 6 3
2 1 4 6 3 5
2 1 4 6 5 3
2 1 5 3 4 6
2 1 5 3 6 4
2 1 5 4 3 6
2 1 5 4 6 3
2 1 5 6 3 4
2 1 5 6 4 3
2 1 6 3 4 5
2 1 6 3 5 4
2 1 6 4 3 5
2 1 6 4 5 3
2 1 6 5 3 4
2 1 6 5 4 3
2 3 1 4 5 6
2 3 1 4 6 5
2 3 1 5 4 6
2 3 1 5 6 4
2 3 1 6 4 5
2 3 1 6 5 4
2 3 4 1 5 6
2 3 4 1 6 5
2 3 4 5 1 6
2 3 4 5 6 1
2 3 4 6 1 5
2 3 4 6 5 1
2 3 5 1 4 6
2 3 5 1 6 4
2 3 5 4 1 6
2 3 5 4 6 1
2 3 5 6 1 4
2 3 5 6 4 1
2 3 6 1 4 5
2 3 6 1 5 4
2 3 6 4 1 5
2 3 6 4 5 1
2 3 6 5 1 4
2 3 6 5 4 1
2 4 1 3 5 6
2 4 1 3 6 5
2 4 1 5 3 6
2 4 1 5 6 3
2 4 1 6 3 5
2 4 1 6 5 3
2 4 3 1 5 6
2 4 3 1 6 5
2 4 3 5 1 6
2 4 3 5 6 1
2 4 3 6 1 5
2 4 3 6 5 1
2 4 5 1 3 6
2 4 5 1 6 3
2 4 5 3 1 6
2 4 5 3 6 1
2 4 5 6 1 3
2 4 5 6 3 1
2 4 6 1 3 5
2 4 6 1 5 3
2 4 6 3 1 5
2 4 6 3 5 1
2 4 6 5 1 3
2 4 6 5 3 1
2 5 1 3 4 6
2 5 1 3 6 4
2 5 1 4 3 6
2 5 1 4 6 3
2 5 1 6 3 4
2 5 1 6 4 3
2 5 3 1 4 6
2 5 3 1 6 4
2 5 3 4 1 6
2 5 3 4 6 1
2 5 3 6 1 4
2 5 3 6 4 1
2 5 4 1 3 6
2 5 4 1 6 3
2 5 4 3 1 6
2 5 4 3 6 1
2 5 4 6 1 3
2 5 4 6 3 1
2 5 6 1 3 4
2 5 6 1 4 3
2 5 6 3 1 4
2 5 6 3 4 1
2 5 6 4 1 3
2 5 6 4 3 1
2 6 1 3 4 5
2 6 1 3 5 4
2 6 1 4 3 5
2 6 1 4 5 3
2 6 1 5 3 4
2 6 1 5 4 3
2 6 3 1 4 5
2 6 3 1 5 4
2 6 3 4 1 5
2 6 3 4 5 1
2 6 3 5 1 4
2 6 3 5 4 1
2 6 4 1 3 5
2 6 4 1 5 3
2 6 4 3 1 5
2 6 4 3 5 1
2 6 4 5 1 3
2 6 4 5 3 1
2 6 5 1 3 4
2 6 5 1 4 3
2 6 5 3 1 4
2 6 5 3 4 1
2 6 5 4 1 3
2 6 5 4 3 1
3 1 2 4 5 6
3 1 2 4 6 5
3 1 2 5 4 6
3 1 2 5 6 4
3 1 2 6 4 5
3 1 2 6 5 4
3 1 4 2 5 6
3 1 4 2 6 5
3 1 4 5 2 6
3 1 4 5 6 2
3 1 4 6 2 5
3 1 4 6 5 2
3 1 5 2 4 6
3 1 5 2 6 4
3 1 5 4 2 6
3 1 5 4 6 2
3 1 5 6 2 4
3 1 5 6 4 2
3 1 6 2 4 5
3 1 6 2 5 4
3 1 6 4 2 5
3 1 6 4 5 2
3 1 6 5 2 4
3 1 6 5 4 2
3 2 1 4 5 6
3 2 1 4 6 5
3 2 1 5 4 6
3 2 1 5 6 4
3 2 1 6 4 5
3 2 1 6 5 4
3 2 4 1 5 6
3 2 4 1 6 5
3 2 4 5 1 6
3 2 4 5 6 1
3 2 4 6 1 5
3 2 4 6 5 1
3 2 5 1 4 6
3 2 5 1 6 4
3 2 5 4 1 6
3 2 5 4 6 1
3 2 5 6 1 4
3 2 5 6 4 1
3 2 6 1 4 5
3 2 6 1 5 4
3 2 6 4 1 5
3 2 6 4 5 1
3 2 6 5 1 4
3 2 6 5 4 1
3 4 1 2 5 6
3 4 1 2 6 5
3 4 1 5 2 6
3 4 1 5 6 2
3 4 1 6 2 5
3 4 1 6 5 2
3 4 2 1 5 6
3 4 2 1 6 5
3 4 2 5 1 6
3 4 2 5 6 1
3 4 2 6 1 5
3 4 2 6 5 1
3 4 5 1 2 6
3 4 5 1 6 2
3 4 5 2 1 6
3 4 5 2 6 1
3 4 5 6 1 2
3 4 5 6 2 1
3 4 6 1 2 5
3 4 6 1 5 2
3 4 6 2 1 5
3 4 6 2 5 1
3 4 6 5 1 2
3 4 6 5 2 1
3 5 1 2 4 6
3 5 1 2 6 4
3 5 1 4 2 6
3 5 1 4 6 2
3 5 1 6 2 4
3 5 1 6 4 2
3 5 2 1 4 6
3 5 2 1 6 4
3 5 2 4 1 6
3 5 2 4 6 1
3 5 2 6 1 4
3 5 2 6 4 1
3 5 4 1 2 6
3 5 4 1 6 2
3 5 4 2 1 6
3 5 4 2 6 1
3 5 4 6 1 2
3 5 4 6 2 1
3 5 6 1 2 4
3 5 6 1 4 2
3 5 6 2 1 4
3 5 6 2 4 1
3 5 6 4 1 2
3 5 6 4 2 1
3 6 1 2 4 5
3 6 1 2 5 4
3 6 1 4 2 5
3 6 1 4 5 2
3 6 1 5 2 4
3 6 1 5 4 2
3 6 2 1 4 5
3 6 2 1 5 4
3 6 2 4 1 5
3 6 2 4 5 1
3 6 2 5 1 4
3 6 2 5 4 1
3 6 4 1 2 5
3 6 4 1 5 2
3 6 4 2 1 5
3 6 4 2 5 1
3 6 4 5 1 2
3 6 4 5 2 1
3 6 5 1 2 4
3 6 5 1 4 2
3 6 5 2 1 4
3 6 5 2 4 1
3 6 5 4 1 2
3 6 5 4 2 1
4 1 2 3 5 6
4 1 2 3 6 5
4 1 2 5 3 6
4 1 2 5 6 3
4 1 2 6 3 5
4 1 2 6 5 3
4 1 3 2 5 6
4 1 3 2 6 5
4 1 3 5 2 6
4 1 3 5 6 2
4 1 3 6 2 5
4 1 3 6 5 2
4 1 5 2 3 6
4 1 5 2 6 3
4 1 5 3 2 6
4 1 5 3 6 2
4 1 5 6 2 3
4 1 5 6 3 2
4 1 6 2 3 5
4 1 6 2 5 3
4 1 6 3 2 5
4 1 6 3 5 2
4 1 6 5 2 3
4 1 6 5 3 2
4 2 1 3 5 6
4 2 1 3 6 5
4 2 1 5 3 6
4 2 1 5 6 3
4 2 1 6 3 5
4 2 1 6 5 3
4 2 3 1 5 6
4 2 3 1 6 5
4 2 3 5 1 6
4 2 3 5 6 1
4 2 3 6 1 5
4 2 3 6 5 1
4 2 5 1 3 6
4 2 5 1 6 3
4 2 5 3 1 6
4 2 5 3 6 1
4 2 5 6 1 3
4 2 5 6 3 1
4 2 6 1 3 5
4 2 6 1 5 3
4 2 6 3 1 5
4 2 6 3 5 1
4 2 6 5 1 3
4 2 6 5 3 1
4 3 1 2 5 6
4 3 1 2 6 5
4 3 1 5 2 6
4 3 1 5 6 2
4 3 1 6 2 5
4 3 1 6 5 2
4 3 2 1 5 6
4 3 2 1 6 5
4 3 2 5 1 6
4 3 2 5 6 1
4 3 2 6 1 5
4 3 2 6 5 1
4 3 5 1 2 6
4 3 5 1 6 2
4 3 5 2 1 6
4 3 5 2 6 1
4 3 5 6 1 2
4 3 5 6 2 1
4 3 6 1 2 5
4 3 6 1 5 2
4 3 6 2 1 5
4 3 6 2 5 1
4 3 6 5 1 2
4 3 6 5 2 1
4 5 1 2 3 6
4 5 1 2 6 3
4 5 1 3 2 6
4 5 1 3 6 2
4 5 1 6 2 3
4 5 1 6 3 2
4 5 2 1 3 6
4 5 2 1 6 3
4 5 2 3 1 6
4 5 2 3 6 1
4 5 2 6 1 3
4 5 2 6 3 1
4 5 3 1 2 6
4 5 3 1 6 2
4 5 3 2 1 6
4 5 3 2 6 1
4 5 3 6 1 2
4 5 3 6 2 1
4 5 6 1 2 3
4 5 6 1 3 2
4 5 6 2 1 3
4 5 6 2 3 1
4 5 6 3 1 2
4 5 6 3 2 1
4 6 1 2 3 5
4 6 1 2 5 3
4 6 1 3 2 5
4 6 1 3 5 2
4 6 1 5 2 3
4 6 1 5 3 2
4 6 2 1 3 5
4 6 2 1 5 3
4 6 2 3 1 5
4 6 2 3 5 1
4 6 2 5 1 3
4 6 2 5 3 1
4 6 3 1 2 5
4 6 3 1 5 2
4 6 3 2 1 5
4 6 3 2 5 1
4 6 3 5 1 2
4 6 3 5 2 1
4 6 5 1 2 3
4 6 5 1 3 2
4 6 5 2 1 3
4 6 5 2 3 1
4 6 5 3 1 2
4 6 5 3 2 1
5 1 2 3 4 6
5 1 2 3 6 4
5 1 2 4 3 6
5 1 2 4 6 3
5 1 2 6 3 4
5 1 2 6 4 3
5 1 3 2 4 6
5 1 3 2 6 4
5 1 3 4 2 6
5 1 3 4 6 2
5 1 3 6 2 4
5 1 3 6 4 2
5 1 4 2 3 6
5 1 4 2 6 3
5 1 4 3 2 6
5 1 4 3 6 2
5 1 4 6 2 3
5 1 4 6 3 2
5 1 6 2 3 4
5 1 6 2 4 3
5 1 6 3 2 4
5 1 6 3 4 2
5 1 6 4 2 3
5 1 6 4 3 2
5 2 1 3 4 6
5 2 1 3 6 4
5 2 1 4 3 6
5 2 1 4 6 3
5 2 1 6 3 4
5 2 1 6 4 3
5 2 3 1 4 6
5 2 3 1 6 4
5 2 3 4 1 6
5 2 3 4 6 1
5 2 3 6 1 4
5 2 3 6 4 1
5 2 4 1 3 6
5 2 4 1 6 3
5 2 4 3 1 6
5 2 4 3 6 1
5 2 4 6 1 3
5 2 4 6 3 1
5 2 6 1 3 4
5 2 6 1 4 3
5 2 6 3 1 4
5 2 6 3 4 1
5 2 6 4 1 3
5 2 6 4 3 1
5 3 1 2 4 6
5 3 1 2 6 4
5 3 1 4 2 6
5 3 1 4 6 2
5 3 1 6 2 4
5 3 1 6 4 2
5 3 2 1 4 6
5 3 2 1 6 4
5 3 2 4 1 6
5 3 2 4 6 1
5 3 2 6 1 4
5 3 2 6 4 1
5 3 4 1 2 6
5 3 4 1 6 2
5 3 4 2 1 6
5 3 4 2 6 1
5 3 4 6 1 2
5 3 4 6 2 1
5 3 6 1 2 4
5 3 6 1 4 2
5 3 6 2 1 4
5 3 6 2 4 1
5 3 6 4 1 2
5 3 6 4 2 1
5 4 1 2 3 6
5 4 1 2 6 3
5 4 1 3 2 6
5 4 1 3 6 2
5 4 1 6 2 3
5 4 1 6 3 2
5 4 2 1 3 6
5 4 2 1 6 3
5 4 2 3 1 6
5 4 2 3 6 1
5 4 2 6 1 3
5 4 2 6 3 1
5 4 3 1 2 6
5 4 3 1 6 2
5 4 3 2 1 6
5 4 3 2 6 1
5 4 3 6 1 2
5 4 3 6 2 1
5 4 6 1 2 3
5 4 6 1 3 2
5 4 6 2 1 3
5 4 6 2 3 1
5 4 6 3 1 2
5 4 6 3 2 1
5 6 1 2 3 4
5 6 1 2 4 3
5 6 1 3 2 4
5 6 1 3 4 2
5 6 1 4 2 3
5 6 1 4 3 2
5 6 2 1 3 4
5 6 2 1 4 3
5 6 2 3 1 4
5 6 2 3 4 1
5 6 2 4 1 3
5 6 2 4 3 1
5 6 3 1 2 4
5 6 3 1 4 2
5 6 3 2 1 4
5 6 3 2 4 1
5 6 3 4 1 2
5 6 3 4 2 1
5 6 4 1 2 3
5 6 4 1 3 2
5 6 4 2 1 3
5 6 4 2 3 1
5 6 4 3 1 2
5 6 4 3 2 1
6 1 2 3 4 5
6 1 2 3 5 4
6 1 2 4 3 5
6 1 2 4 5 3
6 1 2 5 3 4
6 1 2 5 4 3
6 1 3 2 4 5
6 1 3 2 5 4
6 1 3 4 2 5
6 1 3 4 5 2
6 1 3 5 2 4
6 1 3 5 4 2
6 1 4 2 3 5
6 1 4 2 5 3
6 1 4 3 2 5
6 1 4 3 5 2
6 1 4 5 2 3
6 1 4 5 3 2
6 1 5 2 3 4
6 1 5 2 4 3
6 1 5 3 2 4
6 1 5 3 4 2
6 1 5 4 2 3
6 1 5 4 3 2
6 2 1 3 4 5
6 2 1 3 5 4
6 2 1 4 3 5
6 2 1 4 5 3
6 2 1 5 3 4
6 2 1 5 4 3
6 2 3 1 4 5
6 2 3 1 5 4
6 2 3 4 1 5
6 2 3 4 5 1
6 2 3 5 1 4
6 2 3 5 4 1
6 2 4 1 3 5
6 2 4 1 5 3
6 2 4 3 1 5
6 2 4 3 5 1
6 2 4 5 1 3
6 2 4 5 3 1
6 2 5 1 3 4
6 2 5 1 4 3
6 2 5 3 1 4
6 2 5 3 4 1
6 2 5 4 1 3
6 2 5 4 3 1
6 3 1 2 4 5
6 3 1 2 5 4
6 3 1 4 2 5
6 3 1 4 5 2
6 3 1 5 2 4
6 3 1 5 4 2
6 3 2 1 4 5
6 3 2 1 5 4
6 3 2 4 1 5
6 3 2 4 5 1
6 3 2 5 1 4
6 3 2 5 4 1
6 3 4 1 2 5
6 3 4 1 5 2
6 3 4 2 1 5
6 3 4 2 5 1
6 3 4 5 1 2
6 3 4 5 2 1
6 3 5 1 2 4
6 3 5 1 4 2
6 3 5 2 1 4
6 3 5 2 4 1
6 3 5 4 1 2
6 3 5 4 2 1
6 4 1 2 3 5
6 4 1 2 5 3
6 4 1 3 2 5
6 4 1 3 5 2
6 4 1 5 2 3
6 4 1 5 3 2
6 4 2 1 3 5
6 4 2 1 5 3
6 4 2 3 1 5
6 4 2 3 5 1
6 4 2 5 1 3
6 4 2 5 3 1
6 4 3 1 2 5
6 4 3 1 5 2
6 4 3 2 1 5
6 4 3 2 5 1
6 4 3 5 1 2
6 4 3 5 2 1
6 4 5 1 2 3
6 4 5 1 3 2
6 4 5 2 1 3
6 4 5 2 3 1
6 4 5 3 1 2
6 4 5 3 2 1
6 5 1 2 3 4
6 5 1 2 4 3
6 5 1 3 2 4
6 5 1 3 4 2
6 5 1 4 2 3
6 5 1 4 3 2
6 5 2 1 3 4
6 5 2 1 4 3
6 5 2 3 1 4
6 5 2 3 4 1
6 5 2 4 1 3
6 5 2 4 3 1
6 5 3 1 2 4
6 5 3 1 4 2
6 5 3 2 1 4
6 5 3 2 4 1
6 5 3 4 1 2
6 5 3 4 2 1
6 5 4 1 2 3
6 5 4 1 3 2
6 5 4 2 1 3
6 5 4 2 3 1
6 5 4 3 1 2
6 5 4 3 2 1
"""