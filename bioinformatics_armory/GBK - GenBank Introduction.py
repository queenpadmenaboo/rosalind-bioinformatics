"""
Rosalind Problem: GBK - GenBank Introduction
Problem URL: http://rosalind.info/problems/gbk/
Location: Bioinformatics Armory
Author: Bunsree Patel
Date: November 27, 2025
"""

"""Sample Dataset Problem

GenBank comprises several subdivisions:

Nucleotide: a collection of nucleic acid sequences from several sources.
Genome Survey Sequence (GSS): uncharacterized short genomic sequences.
Expressed Sequence Tags, (EST): uncharacterized short cDNA sequences.
Searching the Nucleotide database with general text queries will produce the most relevant results. You can also use a simple query based on protein name, gene name or gene symbol.

To limit your search to only certain kinds of records, you can search using GenBank's Limits page or alternatively use the Filter your results field to select categories of records after a search.

If you cannot find what you are searching for, check how the database interpreted your query by investigating the Search details field on the right side of the page. This field automatically translates your search into standard keywords.

For example, if you search for Drosophila, the Search details field will contain (Drosophila[All Fields]), and you will obtain all entries that mention Drosophila (including all its endosymbionts). You can restrict your search to only organisms belonging to the Drosophila genus by using a search tag and searching for Drosophila[Organism].

Given: A genus name, followed by two dates in YYYY/M/D format.

Return: The number of Nucleotide GenBank entries for the given genus that were published between the dates specified."""

"""Sample Dataset Problem
Anthoxanthum
2003/7/25
2005/12/27"""

from Bio import Entrez
Entrez.email = "bunsree@gmail.com"
handle = Entrez.esearch(
    db="nucleotide",
    term='"Anthoxanthum"[Organism]',
    mindate="2003/7/25",
    maxdate="2005/12/27",
    datetype="pdat",    # filter by publication date
    retmax=1000       # maximum number of IDS to retrieve
)
# Read the search results
record = Entrez.read(handle)

# Print count and IDS
print("Number of sequences founds:", record["Count"])
print("Sequence IDs:", record["IdList"])
"""Output: 
Number of sequences founds: 7
Sequence IDs: ['33413983', '33413982', '33413981', '33413980', '71056585', '57283843', '57283791']"""

"""Actual Dataset Problem"""

from Bio import Entrez
Entrez.email = "bunsree@gmail.com"
handle = Entrez.esearch(
    db="nucleotide",
    term='"Stibaera"[Organism]',
    mindate="2003/12/06",
    maxdate="2012/06/13",
    datetype="pdat",
    retmax=1000
)
# Read the search results
record = Entrez.read(handle)

# Print count and IDS
print("Number of sequences founds:", record["Count"])
print("Sequence IDs:", record["IdList"])

"""Output: Number of sequences founds: 34
Sequence IDs: ['331143314', '331143308', '374976376', '374973383', '374971581', '374967073', '374959593', '374954509', '374953799', '374949379', '374949145', '374949143', '374935292', '374935290', '374909413', '374908321', '374908317', '356491086', '294470319', '294470317', '294470315', '294470313', '294470311', '294470309', '294470307', '294470305', '294470303', '294470301', '291500179', '291500177', '291500175', '291500173', '291500171', '291500169']"""
