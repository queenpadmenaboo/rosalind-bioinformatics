"""
Rosalind Problem: IPRB - Mendel's First Law
Problem URL: http://rosalind.info/problems/iprb/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: October 31, 2025
"""

"""The contemporary heredity model - blending inheritance, states that an organism must exhibit a blend of its parent's traits.
height
over time, blended traits would simply blend into the average, severely limiting variation.

Mendel --> traits should be divided into discrete building blocks called factors.
He proposed that every factor possesses distinct forms called alleles.

Mendel's first law (law of segregation) --> every organism possess a pair of alleles for a given factor.
For any factor, an organism randomly passes one of its two alleles to each offspring, so that an individual receives one allele from each parent.
If two alleles are same, then it is homozygous for the factor.
If the alleles differ, then the individual is heterozygous.

dominant alleles
recessive alleles
genotype
phenotype

Punnett square represents the different possibilities of an individual's inheritance of two alleles from its parents.

Probability
random variable
outcomes
probability tree diagram
event
"""

"""Sample Dataset Problem

Given: Three positive integers k, m, and n, representing a population containing k + m + n organisms;
k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.

Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele (and thus displaying the dominant phenotype).
Assume that any two organisms can mate.

k = 2
m = 2
n = 2

"""
# Mendel's First Law
# Given k (BB), m (Bb), n (bb)
# Return probability that offspring shows dominant phenotype

k, m, n = 2, 2, 2
total = k + m + n

"""Note that we're sampling ordered pairs of individuals.
Every random draw of two distinct individuals has two possible orders, both equally likely.
In order to capture all possible pairs, we must account for both Bb x bb and bb x Bb

    BB x BB: 100% BB (100% dominant)
    BB x Bb: 50% BB and 50% Bb (100% dominant)
    BB x bb: 100% Bb (100% dominant)
    Bb x BB: 50% BB and 50% Bb (100% dominant)
    bb x BB: 100% Bb (100% dominant)

    Bb x Bb: 25% BB, 50% Bb, 25% bb (75% dominant and 25% recessive)
    Bb x bb: 50% Bb, 50% bb (50% dominant and 50% recessive)
    bb x Bb: 50% Bb, 50% bb (50% dominant and 50% recessive)
    bb x bb: 100% bb (100% recessive)"""

# Option 1: Subtract probability of recessive offspring (bb) because only Bb x Bb, Bb x bb, bb x Bb and bb x bb can yield recessive offspring (bb)
p = 1- (
    (m/total)*((m-1)/(total-1))*0.25 +      # heterozygous x heterozygous (Bb x Bb)
    (m/total)*(n/(total-1))*0.5 +           # heterozygous x recessive (Bb x bb)
    (n/total)*(m/(total-1))*0.5 +           # recessive x heterozygous (bb x Bb)
    (n/total)*((n-1)/(total-1))             # both recessive (bb x bb)
)

print(round(p,5))
"""Output: 0.78333"""

# Option 2: Write as direct sum of dominant probabilities

p = (
    (k/total)*((k-1)/(total-1)) +       # homozygous dominant x homozygous dominant (BB)
    (k/total)*(m/(total-1)) +           # homozygous dominant x heterozygous (BB x Bb)
    (k/total)*(n/(total-1)) +           # homozygous dominant x homozygous recessive (BB x bb)
    (m/total)*(k/(total-1)) +           # heterozygous x recessive (Bb x bb)
    (n/total)*(k/(total-1)) +           # recessive x homozygous dominant (bb x BB)
    (m/total)*((m-1)/(total-1))*0.75 +  # heterozygous x heterozygous (Bb x Bb)
    (m/total)*(n/(total-1))*0.5 +       # heterozygous x recessive (Bb x bb)
    (n/total)*(m/(total-1))*0.5         # recessive x heterozygous (bb x Bb)
)

print(round(p,5))

"""No need to include bb x bb, since it gives a '0%' dominant offspring probability."""

"""Actual Dataset Problem"""
# Using Option 1
k, m, n = 24, 19, 22
total = k + m + n

p = 1- (
    (m/total)*((m-1)/(total-1))*0.25 +      # heterozygous x heterozygous
    (m/total)*(n/(total-1))*0.5 +           # heterozygous x recessive (Bb x bb)
    (n/total)*(m/(total-1))*0.5 +           # recessive x heterozygous (bb x Bb)
    (n/total)*((n-1)/(total-1))             # both recessive (bb x bb)
)

print(round(p,5))
"""Output: 0.76791"""

# Using Option 2
k, m, n = 24, 19, 22
total = k + m + n

p = (
    (k/total)*((k-1)/(total-1)) +       # homozygous dominant x homozygous dominant (BB)
    (k/total)*(m/(total-1)) +           # homozygous dominant x heterozygous (BB x Bb)
    (k/total)*(n/(total-1)) +           # homozygous dominant x homozygous recessive (BB x bb)
    (m/total)*(k/(total-1)) +           # heterozygous x recessive (Bb x bb)
    (n/total)*(k/(total-1)) +           # recessive x homozygous dominant (bb x BB)
    (m/total)*((m-1)/(total-1))*0.75 +  # heterozygous x heterozygous (Bb x Bb)
    (m/total)*(n/(total-1))*0.5 +       # heterozygous x recessive (Bb x bb)
    (n/total)*(m/(total-1))*0.5         # recessive x heterozygous (bb x Bb)
)
print(round(p,5))