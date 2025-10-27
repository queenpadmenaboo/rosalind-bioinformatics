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
    Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair)."""

n = 5 
k = 3


