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
"""Extracts sequence as a single string. Use when raw sequence is needed for motifs, translations, etc."""


# Parse FASTA text into a dictionary of {ID: sequence}
def parse_fasta(fasta_text):
    """Parse FASTA text into a dictionary {id: sequence}.

    Args:
        fasta_text (str): Raw FASTA format text containing one or more sequences.
    Returns:
        dict: Keys are FASTA IDs, values are DNA sequences as strings.
    """
    fasta_dict = {}
    current_id = None

    for line in fasta_text.strip().splitlines():
        if line.startswith(">"):
            current_id = line[1:].strip()
            fasta_dict[current_id] = ""
        else:
            fasta_dict[current_id] += line.strip()
    return fasta_dict
"""Extracts all sequences separately, keeping their IDs. Use when handling multiple sequences individually, e.g., for profile matrices, consensus strings, or batch processing."""  


# Find N-glycosylation motif positions in a protein sequence
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

# Function to fetch FASTA from UniProt by accession ID
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
    """status code 200 means Success, so != i.e. the request failed"""

# Monoisotopic mass table in Da (residue masses after peptide bond formation)
monoisotopic_mass = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406,
    'M': 131.04049,
    'N': 114.04293,
    'P': 97.05276,
    'Q': 128.05858,
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333
}

Water_Mass = 18.01056 # mass of one H₂O molecule, used for full protein mass

"""| Amino Acid     | Mass            | Known Cross-check Source                     |
| -------------- | --------------- | -------------------------------------------- |
| Glycine (G)    | 57.02146        | Standard monoisotopic                        |
| Tryptophan (W) | 186.07931       | Largest natural residue                      |
| Leu / Ile      | both ~113.08406 | Correct: **isobaric pair** (same exact mass) |"""

# Function to compute peptide mass (internal mass, no terminal water)
def internal_peptide_mass(protein_seq):
    total = 0.0
    for aa in protein_seq:
        total += monoisotopic_mass[aa]
    return total

# Function to compute full protein mass (includes terminal water)
def full_protein_mass(protein_seq):
    return internal_peptide_mass(protein_seq) + Water_Mass
