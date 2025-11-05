"""
Rosalind Problem: INI4 - Conditions and Loops
Problem URL: http://rosalind.info/problems/ini4/
Location: Python Village
Author: Bunsree Patel
Date: October 20, 2025
"""

"""To choose between two actions, then you can use an if/else statement"""
a = 42
if a < 10:
    print('the number is less than 10')
else:
    print('the number is greater or equal to 10')
"""Note the indentation and punctuation especially the colon at the end of the if and else lines."""

a = 5
b = 3
if a + b == 4:
    print('printed when a + b equals four')
print('always printed')
"""The second print statement is not indented, so it is always executed."""

a = 1
b = 3
if a + b == 4:
    print('printed when a + b equals four')

"""A while loop repeats a block of code as long as a certain condition is true. To repeat an action multiple times, you can use a while loop. 
Note the indentation after the while line. The loop will continue as long as the condition greetings <= 3 is true.
Inside the loop, print 'Hello!' multiplied by the current value of greetings, which results in 'Hello!' being printed an increasing number of times.
Then, increment greetings by 1 using greetings = greetings + 1.  This ensures that the loop will eventually terminate when greetings exceeds 3. 
Be careful to avoid infinite loops, which occur when the condition never becomes false."""

greetings = 1
while greetings <= 3:
    print('Hello!' * greetings)
    greetings = greetings + 1


"""Example of an infinite loop (do not run this code):"""
"""greetings = 1
    while greetings <=3:
    print('Hello!' * greetings)
    greetings = greetings + 0 # Bug here"""
"""In this example, greetings is never incremented, so the condition greetings <= 3 remains true indefinitely, causing an infinite loop where Hello! is printed endlessly."""

"""A for loop is used to carry out some action on every element of a list."""
names = ['Alice', 'Bob', 'Charley']
for name in names:
    print('Hello, ' + name)


"""To repeat an action exactly n times, you can use a for loop with range(n). range is a function that creates a list of integers between 0 and n, where n is not included.
In this example, the loop will iterate over the integers from 0 to 9, printing 'i' ten times.""" 
n = 10
for i in range(n):
    print('i')
 
"""You can also specify a starting point for range. In this example, the loop will iterate over the integers from 5 to 11, printing 'i' seven times."""
for i in range(5, 12):
    print('i')


"""Sample Dataset"""
a = 100 
b = 200 
if a < b < 10000:
    total = 0
    for i in range(a, b + 1):
        if i % 2 == 1: # Check if i is odd
            total += i # Add i to total if it is odd
    print(total)
else:
    print('Invalid input')

"""Actual Dataset Problem"""
a = 4004
b = 8427
if a < b < 10000:
    total = 0
    for i in range(a, b + 1):
        if i % 2 == 1: # Check if i is odd
            total += i # Add i to total if it is odd
    print(total)
else:
    print('Invalid input')

"""Simple solution without explicit check for odd numbers"""
a = 4004
b = 8427
print(sum(range(a|1, b+1, 2)))  # a|1 gives the next odd number if a is even

"""Pseudo-code solution"""
sum=0
for i in range(a, b + 1):
    if i % 2 == 1:
        sum += i
print(sum)