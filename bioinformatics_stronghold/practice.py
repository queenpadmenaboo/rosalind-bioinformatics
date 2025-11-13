
# Write a function that counts how many times 'A' appears in a DNA sequence.
def count_adenine(dna_seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna_seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict["A"]

dna_seq = "AATCGA"
print(count_adenine(dna_seq))
"""Output: 3"""

# Alternate method without dictionary
def count_cytosine(dna_seq):
    count = 0
    for nuc in dna_seq:
       if nuc == 'C':
           count += 1
    return count

dna_seq = "AATCGA"
print(count_cytosine(dna_seq))
"""Output: 1"""


# Even simpler method
def adenine_count(dna_seq):
    return dna_seq.count('A')

dna_seq = "AATCGA"
print(adenine_count(dna_seq))
"""Output: 3"""

dna_seq = "TTCCGG"
print(adenine_count(dna_seq))
"""Output: 0"""

dna_seq = "AAAA"
print(adenine_count(dna_seq))
"""Output: 4"""