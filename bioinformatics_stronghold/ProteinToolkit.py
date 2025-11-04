# Protein Toolkit file
"""
Author: Bunsree Patel
Date: November 3, 2025
"""
import requests
# Allows Python to download data from the internet

# Helper function to clean FASTA text
def clean_fasta(fasta_text):
    """Remove header lines and whitespace from FASTA text.

    Args:
        fasta_text (str): Raw FASTA format text
    Returns:
        str: Cleaned amino acid sequence
    """
    lines = fasta_text.strip().splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines)

# Define function to find N-glycosylation motif positions in a protein sequence
def find_n_glycosylation(seq):
    """Find N-glycosylation motif positions in a protein sequence.

    The N-glycosylation motif: N-{P}-[ST]-{P}
        -N = Asparagine at position 1
        -{P} = any amino acid except Proline at position 2
        -[ST] = Serine or Threonine at position 3
        -{P} = any amino acid except Proline at position 4
        
    Proline has a rigid ring structure that prevents the chain from folding properly, so glycosylation cannot occur if Proline is at positions 2 or 4.
    
    Args:
        seq (str): Amino acid sequence (single-letter codes)
    Returns:
        list of int: 1-based starting indices of motif occurrences   
    """
    positions = []
    for i in range(len(seq) - 3):
        if (seq[i] == 'N' and 
            seq[i + 1] != 'P' and 
            seq[i + 2] in 'ST' and 
            seq[i + 3] != 'P'):
            positions.append(i + 1)
    return positions

# Function to fetch FASTA from UniProt
def fetch_uniprot_fasta(uid):
    """Fetch FASTA sequence from UniProt database. 
    
    Args:
        uid (str): UniProt accession ID
    Returns:
        str: Raw FASTA text
    Raises:
        ValueError: If the UniProt ID cannot be fetched
        requests.RequestException: If network request fails
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Could not fetch UniProt ID '{uid}' (status code: {response.status_code})")
    return response.text