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
    """
    Input: fasta_text (string)
    Output: cleaned sequence (string)
    """
    lines = fasta_text.strip().split("\n")
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence

# Function to fetch FASTA from UniProt
   
def fetch_uniprot_fasta(uid):
    """
    Input: uid (string) - UniProt accesion ID
    Output: fasta_text (string) - raw FASTA text
    """
    url = f"https://uniprot.org/uniprot/{uid}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Could not fetch {uid}: Status {response.status_code}")
    return response.text
  
# Define function to find N-glycosylation motif positions

def find_n_glycosylation(seq):
    """
    Input: seq (string) - cleaned amino acid sequence
    Output: positions (list of int) - 1-based starting indices of motif
    N-glycosylation motif:
    -N = Asparagine
    -{P} = any amino acid except Proline
    -[ST] = either Serine or Threonine
    -{P} = any amino acid except Proline
    -Proline has a rigid ring structure that prevents the chain from folding properly, so glycosylation won't happen if a Proline is in either of the {P} positions
    -That is why the first and last {P} are explicitly there --> they are two different positions in the motif
    """
    motif = r"N[^P][ST][^P]"
    positions = [m.start() + 1 for m in re.finditer(motif, seq)]
    return positions

  