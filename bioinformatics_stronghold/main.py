# DNA Toolset/Code testing file
"""
Author: Bunsree Patel
Date: October 24, 2025
"""

"""Testing validateSeq function from DNAToolkit.py"""
from DNAToolkit import validateSeq
rndDNAStr = "ATTTCGT"
print(validateSeq(rndDNAStr))
# Expected Output: ATTTCGT

"""Testing validateSeq function from DNAToolkit.py, invalid character 'X' in string."""
from DNAToolkit import validateSeq
rndDNAStr = "ATTTCGTX"
print(validateSeq(rndDNAStr))
# Expected Output: False

"""Mixed case (upper and lower case) letters should give string output in all upper case."""
from DNAToolkit import validateSeq
rndDNAStr = "AtCCgGGtGGt"
print(validateSeq(rndDNAStr))
# Expected Output: ATCCGGGTGGT

"""All 3 testing scenarios passed successfully."""


from DNAToolkit import *
import random

# Creating a random DNA string of length 20
randDNAStr = ''.join([random.choice(Nucleotides)
                     for nuc in range(50)])
print(validateSeq(randDNAStr))

print(validateSeq(randDNAStr))
print(countNucFrequency(randDNAStr))