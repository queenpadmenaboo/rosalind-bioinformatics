"""
Rosalind Problem: FIBD - Mortal Fibonacci Rabbits
Problem URL: http://rosalind.info/problems/fibd/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 10, 2025
"""

"""Sample Dataset Problem
Recall the definition of the Fibonacci numbers from “Rabbits and Recurrence Relations”, 
which followed the 'recurrence relation' Fₙ = Fₙ₋₁ + Fₙ₋₂ and assumed that each pair of rabbits reaches maturity in one month and produces a single pair of offspring (one male, one female) each subsequent month.

Our aim is to somehow modify this recurrence relation to achieve a 'dynamic programming' solution in the case that all rabbits die out after a fixed number of months. 
See Figure 4 for a depiction of a rabbit tree in which rabbits live for three months (meaning that they reproduce only twice before dying).

    Given: Positive integers n ≤ 100 and m ≤ 20.
    Return: The total number of pairs of rabbits that will remain after the n-th month if all rabbits live for m months.
    
This problem uses Fₙ = Fₙ₋₁ + Fₙ₋₂ - Fₙ₋₍ₘ₊₁₎, where Fₙ₋₍ₘ₊₁₎ are the rabbits that die.
    
    """


def rabbit_population_mortal(months, lifespan):
    """
    Calculates the total number of rabbit pairs after a given number of months, if each rabbit lives for only a fixed number of months ('lifespan').

    Args:
        months: Positive integer (n ≤ 100) representing the total number of months.
        lifespan: Positive integer (m ≤ 20) representing how many months each rabbit lives.
    
    Returns:
       The total number of rabbit pairs that remain alive after n months.   
    """

    # If no months are simulated, there are no rabbits
    if months == 0:
        return 0

     # Track rabbits born in each month (they live for 'lifespan' months)
    history = [0] * (months + 1)
    history[1] = 1          # Start with 1 pair in month 1

    # Calculate population from month 2 onward
    for month in range(2, months + 1):
        # Rabbits that can reproduce: born from max(1, month-lifespan) to month-2
        history[month] = sum(history[max(1, month - lifespan):month - 1])
        """month - lifespan would be the month when the oldest living rabbits were born. 
            Use max(1, month - lifespan) to ensure we never go below index 1 (since month 1 is the start).
            The full slice gives all the rabbits that are:
               -Old enough to reproduce (at least 1 month old, so not including month - 1)
               -Still alive (not older than lifespan months)"""

    # Total alive extracts rabbits born in the last 'lifespan' months and sums them to get the total living population
    total_alive = sum(history[max(1, months - lifespan + 1):months + 1])
    return total_alive

print(rabbit_population_mortal(6,3))
"""Output: 4"""


"""Actual Dataset Problem"""
print(rabbit_population_mortal(80,16))
"""Output: 23105735806966033"""