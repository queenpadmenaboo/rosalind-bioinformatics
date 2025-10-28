"""
Rosalind Problem: FIB - Rabbits and Recurrence Relations
Problem URL: http://rosalind.info/problems/fib/

Author: Bunsree Patel
Date: October 27, 2025
"""


"""In 1202, Leonardo of Pisa (commonly known as Fibonacci) considered a mathematical exercise regarding the reproduction of a population of rabbits.
He made several assumptions about the population in order to calculate how many pairs of rabbits would remain in one year. After a year, the rabbit populations boasts 144 pairs.

The problem introduces a new computational topic, which involves building up large solutions from smaller ones.

A sequence is an ordered collection of objects (usually numbers), which area allowed to repeat.
Sequences can be finite or infinite. Two examples are the finite sequence (π,–√2,0,π) and the infinite sequence of odd numbers (1,3,5,7,9,...).

Use the notation aₙ to represent the n-th term of a sequence.

A recurrence relation is a way of defining the terms of a sequence with respect to the values of previous terms.
In the case of Fibonacci's rabbits from the introduction, any given month will contain the rabbits that were alive the previous month, plus any new offspring. A key observation is that the number of offspring in any month is equal to the number of rabbits that were alive two months prior.
As a result, if Fn represents the number of rabbit pairs alive after the n-th month, then we obtain the Fibonacci sequence having terms Fn that are defined by the recurrence relation Fₙ = Fₙ₋₁ + Fₙ₋₂ (with F₁ = F₂ = 1 to initiate the sequence). 
Although the sequence bears Fibonacci's name, it was known to Indian mathematicians over two millennia ago."""


"""Sample Dataset Problem
    
When finding the n-th term of a sequence defined by a recurrence relation, use the recurrence relation to generate terms for progressively larger values of n.
The problem introduces dynamic programming, which successively builds up solutions by using the answers to smaller cases.

    Given: Positive integers n ≤ 40 and k ≤ 5.
    Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).
    n = 5 months
    k = 3 rabbit pairs

"""
def Fibonacci_Loop_Pythonic(months, offsprings):
    previous, current = 1, 1
    for _ in range(2, months):
        previous, current = current, current + (previous * offsprings)
    return current

print(Fibonacci_Loop_Pythonic(6,2))
print(Fibonacci_Loop_Pythonic(35, 5))

"""
range (2, months) means:
range(start, stop) generates a sequence of integers starting at 'start' abd stopping before 'stop'.
So range(2, months) produces:
    2, 3, 4, ..., months - 1
    for example range (2, 6) --> [2, 3, 4, 5], the loop will run months - 2 times (since the first two months are already defined by base cases).

**Note, index 2 corresponds to the 3rd month - i.e. the next month after our base cases.**
Python loop index 0 | Corresponding month (not used)
Python loop index 1 | Corresponding month (not used)
Python loop index 2 | Corresponding month - Month 3
Python loop index 3 | Corresponding month - Month 4
Python loop index 4 | Corresponding month - Month 5
Python loop index 5 | Corresponding month - Month 6

Why it starts at 2:
    The Fibonacci-style model that was written starts with Month 1 and Month 2 already known:

    Month 1 --> 1 pair (the first newborns)
    Month 2 --> 1 pair (they've matured, still only 1 pair total)

    From Month 3 onward, we apply the recurrence rule:
    Fₙ = Fₙ₋₁ + (Fₙ₋₂ x k); (with F₁ = F₂ = 1 to initiate the sequence).

c - small (children) rabbits. They have to mature and reproduce in the next cycle only.
C - mature (parents) rabbits. They can reproduce and move to the next cycle.

Month 1: [c]
Month 2: [C]
Month 3: [C c c]
Month 4: [C c c C C]
Month 5: [C c c C C C c c C c c]
Month 6: [C c c C C C c c C c c C c c C C C c c C C]

Output (5,3): 19
Output (6,2): 21

Output (35,5): 875089148811941
"""

"""Actual Dataset Problem"""

# Code rewritten with clearer variable names

def rabbit_population(months, offspringpairs):
    # Base cases: month 1 and 2 each start with 1 pair
    prev_gen = 1 # rabbits from two months ago (mature)
    current_gen = 1 # rabbits from last month (total so far)

    # Calculate population from month 3 onward
    for month in range(2, months):
        prev_gen, current_gen = current_gen, current_gen + (prev_gen *(offspringpairs))
    return current_gen

print(rabbit_population(5,3))
print(rabbit_population(35,5))

"""Output: 875089148811941"""