"""
Rosalind Problem: INI2 - Variables and Some Arithmetic
Problem URL: http://rosalind.info/problems/ini2/
Location: Python Village
Author: Bunsree Patel
Date: October 15, 2025
"""

"""A variable is a name that refers to a value stored in the computer's memory. You can create a variable by assigning it a value using the equals sign (=).
In Python, the basic data types are strings and numbers. Strings are sequences of characters enclosed in quotes, while numbers can be integers or floating-point values."
There are two types of numbers: integers (whole numbers, both positive and negative) and floats (fractional numbers with a decimal point)."""


a = 324
b = 24
c = a - b
print('a - b is', c)
"""a, b, and c are all all integers, and 'a-b' is a string"""


"""Common arithmetic operations include addition (+), subtraction (-), multiplication (*), division (/), modulus or division remainder (%), and exponentiation (**)
Here are some examples of these operations:"""
2 + 3 == 5
5 - 2 == 3
3 * 4 == 12
15 / 3 == 5

"""If you divide two integers, Python always rounds down the result to the nearest integer.
To obtain a precise results with decimal points, you can convert one of the integers to a float using the float() function. or express as decimal i.e. 18.0 instead of 18."""

18 % 5 == 3
18.0/5 == 3.6
float(18)/5 == 3.6

"""To multiply a number by an exponent, you can use the exponentiation operator (**). 
For example, 2 raised to the power of 3 is calculated as follows:"""

2 ** 3 == 8

"""A string is an ordered sequence of letters, numbers, and other characters."""
a = "Hello"
b = "World"
"""Notice that the strings are enclosed in double quotes (" "). 
You can also use single quotes (' '), but be consistent."""

a = 'Rosalind'
b = 'Franklin'
c = '!'
print(a + ' ' + b + c*3)


"""Sample Dataset Problem
Given: Two positive integers a and b, each less than 1000.
Return: The integer corresponding to the square of the hypotenuse of the right triangle whose legs have lengths a and b."""

a = 3
b = 5
c = ( a ** 2 + b ** 2 )
print(c)


"""Actual Dataset Problem"""
a = 842
b = 926
c = ( a ** 2 + b ** 2)
print(c)

"""Pythagoras' Theorem, where c is the hypotenuse squared given two sides a and b."""