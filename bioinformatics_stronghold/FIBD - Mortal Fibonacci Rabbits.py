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


def rabbit_population_mortal(months, lifespan, offspringpairs):
    """
    Calculates the total number of rabbit pairs after a given number of months, if each rabbit lives for only a fixed number of months ('lifespan').

    Args:
        months: Positive integer (n ≤ 100) representing the total number of months.
        lifespan: Positive integer (m ≤ 20) representing how many months each rabbit lives.
    
    Returns:
       The total number of rabbit pairs that remain alive after n months.   
    """

    # Base cases: month 1 and 2 each start with 1 pair
    prev_gen = 1        # rabbits from two months ago (mature)
    current_gen = 1     # rabbits from last month (total so far)

    # Keep track of total rabbits per month to handle mortality
    history = [1,1]     # population record for each month

    # Calculate population from month 3 onward
    for month in range(2, months):
        if month < lifespan:
            # Before any rabbits die → normal Fibonacci growth
            next_gen = current_gen + prev_gen
        else:
            # Rabbits that die = those born lifespan months ago
            next_gen = current_gen + prev_gen - history[-(lifespan + 1)]

                # current_gen: total pairs alive last month
                # prev_gen: number of mature pairs that reproduced this month
                # history[-(lifespan + 1)]: rabbits born lifespan months ago, who die this month
           
        # Update for next iteration
        prev_gen, current_gen = current_gen, next_gen
        history.append(next_gen)

    return current_gen

    print(rabbit_population_mortal(6,3))
