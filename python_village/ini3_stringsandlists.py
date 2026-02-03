"""
Rosalind Problem: INI3 - Strings and Lists
Problem URL: http://rosalind.info/problems/ini3/
Location: Python Village
Author: Bunsree Patel
Date: October 16, 2025
"""

"""There are variable types that can store multiple values, the simplest of which is a list.
Assign data to a list using square brackets [] and separate the values with commas i.e. list_name = [item_1, item_2, ..., item_n].
The items of the list can be of any other type: integer, float, string, or even lists of lists.
Any item in a list can be accessed by its index, with the first item having index 0, the second item index 1, and so on. 
Index indicates its place in the list."""

tea_party = ['March Hare', 'Hatter', 'Dormouse', 'Alice']
print(tea_party[2])
"""Output: Dormouse. It was not Hatter, because in Python, indexing begins at 0.
It is easy to change existing list items by reassigning them."""

tea_party[1] = 'Cheshire Cat'
print(tea_party)
"""Output: ['March Hare', 'Cheshire Cat', 'Dormouse', 'Alice']. Hatter was replaced by Cheshire Cat."""

# Use the append(): function to add items to the end of an existing list
tea_party.append('Jabberwocky')
print(tea_party)
"""Output: ['March Hare', 'Cheshire Cat', 'Dormouse', 'Alice', 'Jabberwocky']. Jabberwocky was added to the end of the list.

To obtain only some of a list, use the notation list_name[a:b] (list_name from start index a to end index b-1).
This is called list slicing. The first index is inclusive, while the second index is exclusive."""

tea_party[1:3]
print(tea_party[1:3])
"""Output: ['Cheshire Cat', 'Dormouse']. Only the items at index 1 and 2 were returned.

If the first index of the slice is unspecified, it defaults to 'index 0' (the start of the list).
If the second index is unspecified, it defaults to the length of the list (the end of the list)."""

tea_party[:2]
print(tea_party[:2])
"""Output: ['March Hare', 'Cheshire Cat']. Only the items at index 0 and 1 were returned."""

tea_party[3:]
print(tea_party[3:])
"""Output: ['Alice', 'Jabberwocky']. Only the items at index 3 and beyond were returned.

Negative indices can be used to index a list from the end (backtracking)."""
tea_party[-2:]
print(tea_party[-2:])
"""Output: ['Alice', 'Jabberwocky']. The last two items of the list were returned. 
tea_party[-2:] returns the same output as tea_party[3:].

Strings in Python can be sliced in the same way as lists. A string is a sequence of characters, and each character has an index starting from 0."""

a = 'flimsy'
b = 'miserable'
c = b[0:1] + a[2:]
print(c)
"""Output: 'mimsy'. The first character of b ('m') was concatenated with the substring of a starting from index 2 ('imsy')."""

"""Sample Dataset Problem
Given: A string s of length at most 200 letters and four integers a, b, c, and d.
Return: The slice of this string from indices a through b and c through d (with space in between), inclusively. In other words, we should include elements s[b] and s[d] in our slice."""

user_string = "HumptyDumptysatonawallHumptyDumptyhadagreatfallAlltheKingshorsesandalltheKingsmenCouldntputHumptyDumptyinhisplaceagain"

a = 22
b = 27
c = 97
d = 102

first_part = user_string[a:b+1]
second_part = user_string[c:d+1]

print(first_part, second_part)
""" +1 is added to include the character at the end index since slicing is exclusive of the end index.
Output: Humpty Dumpty"""

"""Actual Dataset Problem"""
user_string = (
    "i3yiUy33wO9Ze0muK7oK9R3lt8NlBzYe5ip73quDEx8JJajhdnygfG8HZQwsxCircaetus4BYhEVxsNMJcm6Sk1mGaqOxoLpG9sU"
    "JKR0bJgi1cdahuricaZWnb4SGlbUV8OKalQy3eeh57tlPwWl9EMoLX37U8oTDPm5Xtpv5zUyH9r3rYrlCydYjCd57KIbX"
)

a = 61
b = 69
c = 110
d = 117

first_part = user_string[a:b+1]
second_part = user_string[c:d+1]

print(first_part, second_part)

"""Output: Circaetus dahurica.
I have split the long string into multiple lines for better readability. 
Also had errors when copying and pasting long strings in Powershell, which added a '/' character in between."""