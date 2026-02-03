"""
Author: Bunsree Patel
Date: October 27, 2025
"""

def Fibonacci_Loop(number):
    old = 1
    new = 1
    for itr in range(number - 1):
        tmpVal = new
        new = old
        old = old + tmpVal

"""If you have F₁ = 0, Output for Fibonacci_Loop(10) would be 21.
    If you have  F₁ = 1, Output for Fibonacci_Loop(10) would be 55."""