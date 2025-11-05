"""
Rosalind Problem: PRTM - Calculating Protein Mass
Problem URL: http://rosalind.info/problems/prtm/
Location: Bioinformatics Stronghold
Author: Bunsree Patel
Date: November 5, 2025
"""

"""A protein is a 'polypeptide', a chain of amino acids linked by peptide bonds.
Each 'peptide bond' forms when two amino acids link together during protein translation, this reaction releases a molecule of H₂O.
After losing H₂O, amino acids in the chain are called 'residues'.

A polypeptide = n amino acids --> n-1 peptide bonds --> n-1 H₂O molecules removed

The 'monoisotopic mass' of a residue = the mass calculated using the principal (most abundant) 'isotope' of each atom in the amino acid.
The 'average mass' of a residue = the mass calculated by averaging the naturally appearing isotopes of each atom.

In 'mass spectrometry', monoisotopic mass is used more often than average mass;
thus amino acid masses are assumed to be monoisotopic.

Internal peptide mass = sum(monoisotopic residue masses) in Da 'dalton', where 1 Da is defined as 1/12th of the mass of a neutral carbon-12 atom.

For a full protein, total mass = sum(residue masses) + mass of one H₂O (18.01056 Da)."""


"""Sample Dataset Problem

In a 'weighted alphabet', every symbol is assigned a positive real number called a 'weight'. A string formed from a weighted alphabet is called a 'weighted string', and its 'weight' is equal to the sum of the weights of its symbols.

The standard weight assigned to each member of the 20-symbol amino acid alphabet is the monoisotopic mass of the corresponding amino acid.

    Given: A protein string P of length at most 1000 aa.
    Return: The total weight of P. Consult the monoisotopic mass table.
"""
import from ProteinToolkit internal_peptide_mass, full_protein_mass



print(internal_peptide_mass("ACDE"))