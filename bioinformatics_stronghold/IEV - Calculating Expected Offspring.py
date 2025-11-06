"""
Rosalind Problem: IEV - Calculating Expected Offspring
Problem URL: http://rosalind.info/problems/iev/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 6, 2025
"""

""" 'alleles'
    'random variable'
    'expected value'
    'uniform random variable'
    'genotype'
    'factor'

"""

"""Sample Dataset Problem
    Given: Six nonnegative integers, each of which does not exceed 20,000. The integers correspond to the  number of couples in a population possessing each genotype pairing for a given factor.
    In order, the six given integers represent the number of couples having the following genotypes:
        1. AA-AA
        2. AA-Aa
        3. AA-aa
        4. Aa-Aa
        5. Aa-aa
        6. aa-aa
    Return: The expected number of offspring displaying the dominant phenotype in the next generation, under the assumption that every couple has exactly two offspring.
    """

from ProteinToolkit import calculate_expected_dominant_offspring

def iev():
    input_line = input().strip()
    couples = list(map(int, input_line.split()))

    # Calculated expected value
    result = calculate_expected_dominant_offspring(couples)

    # Print result
    print(result)

sample_input = [1, 0, 0, 1, 0, 1]
print(f"Sample input: {sample_input}")
print(f"Expected dominant offspring: {calculate_expected_dominant_offspring(sample_input)}")
"""Output: 
Sample input: [1, 0, 0, 1, 0, 1]
Expected dominant offspring: 3.5"""



"""Actual Dataset Problem"""

sample_input = [16148, 17391, 16865, 19555, 16960, 18881]
print(f"Sample input: {sample_input}")
print(f"Expected dominant offspring: {calculate_expected_dominant_offspring(sample_input)}")
"""Output:
Sample input: [16148, 17391, 16865, 19555, 16960, 18881]
Expected dominant offspring: 147100.5"""