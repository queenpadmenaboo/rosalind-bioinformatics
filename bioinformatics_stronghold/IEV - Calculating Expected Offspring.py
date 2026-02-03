"""
Rosalind Problem: IEV - Calculating Expected Offspring
Problem URL: http://rosalind.info/problems/iev/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 6, 2025
"""

"""
Key Concepts:
    - A 'random variable' is a value that changes based on chance - it's the results of something random, like rolling a die.
    - An 'expected value' is the average value that is expected if a random process is repeated many times.
    - A 'uniform random variable' is a random value where every possible outcome is equally likely.
"""

"""Sample Dataset Problem
    Given: Six nonnegative integers, each of which does not exceed 20,000.
    The integers correspond to the number of couples in a population possessing each genotype pairing for a given factor.
    
    In order, the six given integers represent the number of couples having the following genotypes:
        1. AA-AA
        2. AA-Aa
        3. AA-aa
        4. Aa-Aa
        5. Aa-aa
        6. aa-aa
    
    Return: The expected number of offspring displaying the dominant phenotype in the next generation, 
    under the assumption that every couple has exactly two offspring.
"""

from ProteinToolkit import calculate_expected_dominant_offspring

"""Algorithm: 'calculate_expected_dominant_offspring' function uses:
    1. Define a list of probabilities representing the chance that each genotype pairing produces a dominant-phenotype offspring:
        [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    2. Set offspring_per_couple = 2 (each couple produces exactly 2 offspring).
    3. Initialize 'expected_value = 0' to store the running total of expected dominant offspring.
    4. Loop through each of the 6 genotype pairings using their index i:
        - Multiply inline:
            
            couples[i]  * offspring_per_couple * probabilities[i]
            
            This calculates the expected number of dominant offspring from the i-th genotype pairing.
        - Add the result directly to 'expected_value':

            expected_value += couples [i] * offspring_per_couple * probabilities[i]

    5. After the loop finishes, return expected_value as the total expected number of dominant phenotype offspring."""

def iev():
    """
    Reads input of six integers representing genotype pairings;
    calculates the expected number of offspring with the dominant phenotype, and prints the results.
    """

    input_line = input().strip()
    couples = list(map(int, input_line.split()))

    # Calculated expected value (dominant phenotype offspring)
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