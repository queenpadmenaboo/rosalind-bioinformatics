"""
Rosalind Problem: MRNA - Inferring mRNA from Protein
Problem URL: http://rosalind.info/problems/mrna/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 6, 2025
"""

"""Sample Dataset Problem
For positive integers a and n, a modulo n (written a  mod  n in shorthand) is the remainder when a  is divided by n.
For example, 29  mod  11 = 7 because 29 = 11 x 2 + 7.

Modular arithmetic is the study of addition, subtraction, multiplication, and division with respect to the modulo operation. We say that a and b are 'congruent' modulo n if a  mod  n = b  mod  n; 
in this case, we use the notation a ≡ b  mod  n.

Two useful facts in modular arithmetic are that if a ≡ b  mod  n and c ≡ d  mod  n, then a + c ≡ b + d  mod  n
and a x c ≡ b x d  mod  n. To check your understanding of these rules, you may wish to verify these relationships for a = 29, b = 73, c = 10, d = 32, and n = 11.
As you will see in this exercise, some Rosalind problems will ask for a (very large) integer solution modulo a smaller number to avoid the computational pitfalls that arise with storing such large numbers.

    Given: A protein string of length at most 1000 aa.
    Return: The total number of different RNA strings from which the protein could have been translated, modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)"""

from ProteinToolkit import count_rna_strings

sample_input = "MA"

protein = sample_input.strip()
result = count_rna_strings(protein)
print(f"Total # of different RNA strings: {result}")
"""Output: 
    Total # of differnt RNA strings: 12"""

sample_input = "MTVKYATRYWDYQACEVCHQDCISQWTVPLTCKMQETSWCDNRWVHLPLCYMISNTFCERGAWLWFEAFQDIDVITWRYFEMGWNNVARQMNNVPSRIFWRNANPESHKDGADEMFHGVCWGGLWETGSYVTGFFMCSMYARYTASVAFIKYHSCCGMLGMVGCFWQSMDIADKYTMMFLSVKQMLCARLNCNYTSEYECADIPYHEKPFWNHWLTVGVPHQMPIEWAAKQFSGEQQCTSFVCKQYVAHTFFERWFDYKTYCKIGQKVHQKHNTWIDQFMETCYVVYHDSSIQTMNMNCTYDEVMPAGRWQAKYETMNKVKHMFWSKLYMPAVTRRIVHVWFRGYVRVLFVSRVHSEPLRKSKTTSFDSYCWYMRTPPYQCDGFAKEPEEQRMCPGMWEYNFHMVFTSWHEWLYKYWPAQMATHEVAGGYWLWIQWDRPLMFITKNCSARHILMPSVMWAQIPHMVPCWTQNVDCRDRPSLCMKWYPWVPQMGITPRKEGKYYANDFCFVIQYRPQTSGKWGVNMCCEYARILLHTACGIYAHQAFPTCIEGFNIFEGEWFEKKIPCLLGSNSVTVTNPYLKADMEWNGDWHPGRRHWDLENMVGIPSKAACEERHGHCVILKLLGKSKWRLWWEKRWDMQRMGKCFLMSRPMLDWEKEYWYEIGRNLGYFIEAVEEDGQLIVMHEIDFGNSIFHIFIFQVVPIMAVTYIPWFFICKEDCMSYACSRQRSPMGKGVLWMETYPLKCHCRVRKCEANSGYFCQKFAKSKHDKQKTFKVWMYPVQWREVRQKWDHKAAYGVPRGVKIYFASAHAYRYYVTNQSGWCEQMQQADSYTIKTHNCFEMRDNGTREFVGHRWLKMSEQDKANTSNYIDEMHVQFDVWEDNNASQLHKHSYPPYFGKMCFYWILTNWWRDDMQRYMSKWLYDQEPQHRGADRNKEWKMYYWYMLQEKYMAAKDWRSHAFNLPVYIKADWQGASGMGKFYVHFNHIKMMALWDDWP"
protein = sample_input.strip()
result = count_rna_strings(protein)
print(f"Total # of different RNA strings: {result}")
"""Output: 
    Total # of differnt RNA strings: 71744"""
