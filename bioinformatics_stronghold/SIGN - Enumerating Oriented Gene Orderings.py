"""
Rosalind Problem: SIGN - Enumerating Oriented Gene Orderings
Problem URL: http://rosalind.info/problems/sign/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 8, 2025
"""


"""
Genome Rearrangment Orientation Concepts: 

Orientation refers to the direction in which DNA is read or transcribed. 
    RNA transcription only occurs in one direction along the DNA strand (5' to 3' direction). 
    This directionality is critical because genes have a "forward" and "reverse" orientation on chromosomes.

'Synteny block orientation' becomes important when studying genome rearrangements because a block of genes can be flipped (inverted) during evolution. 
    A block might contain the same genes in the same order, but they could be read in the opposite direction.

'Adding orientations to synteny blocks' expands the concept of permutations. 
    Instead of just tracking the order of numbered blocks like [1, 2, 3, 4], each block now carries a sign indicating its orientation. 
    A positive sign represents forward orientation, while a negative sign represents reverse orientation.

'Signed permutation' is an ordering of positive integers {1, 2, ..., n} where each integer is assigned either a positive or negative sign. 
    The sign indicates the orientation of that element. 
    For example, π = (5, -3, -2, 1, 4) is a signed permutation of length 5, where blocks 3 and 2 are in reverse orientation while blocks 5, 1, and 4 are in forward orientation.
    Example: [1, -2, 3, 4] means block 2 has been inverted compared to the reference genome. The genes within block 2 are still together, but they're now read in the opposite direction.

Comparing signed permutations between species reveals not just which blocks moved, but also which blocks got flipped during evolutionary rearrangements."""

"""Sample Dataset Problem
A 'signed permutation' of length n is some ordering of the positive integers {1,2,...,n} 
in which each integer is then provided with either a positive or negative sign (for the sake of simplicity, we omit the positive sign).
For example, π=(5,-3,-2,1,4) is a signed permutation of length 5.

    Given: A positive integer n ≤ 6.
    Return: The total number of signed permutations of length n, followed by a list of all such permutations (you may list the signed perms. in any order)."""

def signed_permutation(arr):
    # Base case: single element
    if len(arr) == 1:
        return [[arr[0]], [-arr[0]]]    # positive and negative
    
    result = []
    # For each element in an array
    for i in range(len(arr)):
        # Take element at index i
        current = arr[i]
        # Get remaining elements (exclude current)
        remaining = arr [:i] + arr[i+1:]
        # Recursively permute remaining elements
        for p in signed_permutation(remaining):
            # Add current element with positive sign
            result.append([current] + p)
            # Add current element with negative sign
            result.append([-current] + p)

    return result

# Input Sample Data
n = 2

# Generate all permutations of [1, 2, ..., n]
perms = signed_permutation(list(range(1, n+1)))

# Output: total count
print(len(perms))

# Output: each permutation as space-seperated integers
for p in perms:
    print(' '.join(map(str, p)))
"""Output:
8
1 2
-1 2
1 -2
-1 -2
2 1
-2 1
2 -1
-2 -1
"""

"""Actual Dataset Problem"""
# Input Sample Data
n = 3
perms = signed_permutation(list(range(1, n+1)))
print(len(perms))
for p in perms:
    print(' '.join(map(str, p)))
"""Output:
48
1 2 3
-1 2 3
1 -2 3
-1 -2 3
1 2 -3
-1 2 -3
1 -2 -3
-1 -2 -3
1 3 2
-1 3 2
1 -3 2
-1 -3 2
1 3 -2
-1 3 -2
1 -3 -2
-1 -3 -2
2 1 3
-2 1 3
2 -1 3
-2 -1 3
2 1 -3
-2 1 -3
2 -1 -3
-2 -1 -3
2 3 1
-2 3 1
2 -3 1
-2 -3 1
2 3 -1
-2 3 -1
2 -3 -1
-2 -3 -1
3 1 2
-3 1 2
3 -1 2
-3 -1 2
3 1 -2
-3 1 -2
3 -1 -2
-3 -1 -2
3 2 1
-3 2 1
3 -2 1
-3 -2 1
3 2 -1
-3 2 -1
3 -2 -1
-3 -2 -1
"""