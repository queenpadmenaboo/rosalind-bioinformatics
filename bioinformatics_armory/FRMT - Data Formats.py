"""
Rosalind Problem: FRMT - Data Formats
Problem URL: http://rosalind.info/problems/frmt/
Location: Bioinformatics Armory
Author: Bunsree Patel
Date: November 28, 2025
"""

"""
# Popular Data Formats:
- FASTA (.fas, .fasta)
- NEXUS (.nex, .nexus, .nxs)
- PHYLIP "Phylogeny Interference Package" (.phy)

# FASTA
- string is introduced by a line that begins with '>', followed by string label i.e. Rosalind_2496
- the next line that begins with '>' completes the current string and initiates next string label

# GenBank https://www.ncbi.nlm.nih.gov/genbank/
- world's largest genetic database, hosted by NCBI
- contains an annotated collection of all publicly available nuc strings and their protein translations.
- constantly updated from individual labs and large-scale sequencing centers
- Example for a Protein:
    Accession: 'A2Z669'
- Note: 'Send to:' drop-down at the top right of the page for exporting an entry to a variety of file formats
"""

"""Sample Dataset Problem

GenBank can be accessed at link above.
A detailed description of the GenBank format can be found here.
A tool, from the SMS 2 package, for converting GenBank to FASTA can be found at https://www.bioinformatics.org/sms2/genbank_fasta.html.
    Given: A collection of n (n≤10) GenBank entry IDs.

    Return: The shortest of the strings associated with the IDs in FASTA format."""

from Bio import Entrez
Entrez.email = "bunsree@gmail.com"
handle = Entrez.efetch(db="nucleotide", id=["FJ817486, JX069768, JX469983"], rettype="fasta")
records = handle.read()
print records