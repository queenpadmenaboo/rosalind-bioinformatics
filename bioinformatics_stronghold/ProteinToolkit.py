# Protein Toolkit file
"""
Author: Bunsree Patel
Date: November 3, 2025
"""
import requests
# Allows Python to download data from the internet
import re 
# 're' is Python's built-in module for working with regular expressions

# Helper function to clean FASTA text
def clean_fasta(fasta_text):
    lines = fasta_text.strip().splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines)

# Define function to find N-glycosylation motif positions
def find_n_glycosylation(seq):
    positions = []
    for i in range(len(seq) - 3):
        if seq[i] == 'N' and seq[i + 1] != 'P' and seq[i + 2] in 'ST' and seq[i+3] != 'P':
            positions.append(i +1)
    return positions

"""Input: seq (string) - cleaned amino acid sequence
        Output: positions (list of int) - 1-based starting indices of motif
        N-glycosylation motif:
            -N = Asparagine
            -{P} = any amino acid except Proline
            -[ST] = either Serine or Threonine
            -{P} = any amino acid except Proline
            -Proline has a rigid ring structure that prevents the chain from folding properly, so glycosylation won't happen if a Proline is in either of the {P} positions
            -That is why the first and last {P} are explicitly there --> they are two different positions in the motif."""

# Function to fetch FASTA from UniProt
   
def fetch_uniprot_fasta(uid):
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Could not fetch {uid}")
    return response.text
    
""" Input: uid (string) - UniProt accesion ID
    Output: fasta_text (string) - raw FASTA text"""


