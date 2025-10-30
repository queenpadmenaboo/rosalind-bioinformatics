# DNA Toolkit file
"""
Author: Bunsree Patel
Date: October 24, 2025
"""

# Define nucleotide sets
DNA_Nucleotides = ['A', 'C', 'G', 'T']
RNA_Nucleotides = ['A', 'C', 'G', 'U']


# Check the sequence to make sure it is a valid DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in DNA_Nucleotides:
            return False
    return tmpseq

def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def transcribeDNAtoRNA(dna_seq):
    # Validate input
    for nuc in dna_seq.upper():
        if nuc not in DNA_Nucleotides:
            raise ValueError("Invalid DNA sequence")
    
    # Replace T with U to make RNA
    rna_seq = dna_seq.upper().replace("T", "U")
    return rna_seq

def complementDNA(dna_seq):
    # Validate input
    for nuc in dna_seq.upper():
        if nuc not in DNA_Nucleotides:
            raise ValueError("Invalid DNA sequence")
    
    # Complement mapping
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Build complement strand
    comp_seq = ''.join([complement[nuc] for nuc in dna_seq.upper()])
    return comp_seq

def reversecomplementDNA(dna_seq):
    # Validate input
    dna_seq = validateSeq(dna_seq)

    # Complement mapping
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Build reverse complement strand
    rev_complementDNA = ''.join([complement[nuc] for nuc in dna_seq[::-1]])
    return rev_complementDNA

def gc_content(dna_seq):
    dna_seq = validateSeq(dna_seq)
    # Returns GC content as a percentage
    gc = dna_seq.count('G') + dna_seq.count('C')

    return (gc/ len(dna_seq)) * 100
"""len means length, it's a built-in Python function that returns the number of items in something."""

def parse_fasta(fasta_str):
    # Parses a FASTA string into a dictionary {id: dna_seq}
    fasta_dict = {}     # Initialize empty dictionary to store sequences
    current_id = None   # Tracks the current sequence ID while parsing

    for line in fasta_str.strip().splitlines():     # Loop through each line
        if line.startswith(">"):                    # If line starts with '>', it's a sequence ID
            current_id = line[1:].strip()           # Remove '>' and whitespace to get the ID
            fasta_dict[current_id] = ""             # Initialize empty string for this sequence
        else:
            fasta_dict[current_id] += line.strip()  # Append sequence lines to current ID
            
    return fasta_dict                               # Return dictionary of all sequences

"""In Python, a dictionary is a data structure that stores data as key-value pairs.
    Example: person = {"name": "Z", "age": 24}
    "name" --> key
    "Z" --> value
    "age" --> key
    "24" --> value

    The keys are: "Rosalind_8375", "Rosalind_4074", "Rosalind_0002", etc.

    The values are: the long DNA sequence strings made of A,T,G,C."""

def hamming_distance(s, t):
    """
    Calculate the Hamming distance between two DNA strings of equal length.

    Args:
        s (str): First DNA string.
        t (str): Second DNA string.

    Returns:
        int: Number of positions where s and t differ.
    """
    if len(s) != len(t):
        raise ValueError("Strings must be of equal length to compute Hamming distance.")
    
    distance = 0
    for a, b in zip(s,t):
            if a != b:
                distance += 1
    return distance


# RNA codon translation table
rnacodon_table = {
    'UUU':'F', 'CUU':'L', 'AUU':'I', 'GUU':'V',
    'UUC':'F', 'CUC':'L', 'AUC':'I', 'GUC':'V',
    'UUA':'L', 'CUA':'L', 'AUA':'I', 'GUA':'V',
    'UUG':'L', 'CUG':'L', 'AUG':'M', 'GUG':'V',
    'UCU':'S', 'CCU':'P', 'ACU':'T', 'GCU':'A',
    'UCC':'S', 'CCC':'P', 'ACC':'T', 'GCC':'A',
    'UCA':'S', 'CCA':'P', 'ACA':'T', 'GCA':'A',
    'UCG':'S', 'CCG':'P', 'ACG':'T', 'GCG':'A',
    'UAU':'Y', 'CAU':'H', 'AAU':'N', 'GAU':'D',
    'UAC':'Y', 'CAC':'H', 'AAC':'N', 'GAC':'D',
    'UAA':'Stop', 'CAA':'Q', 'AAA':'K', 'GAA':'E',
    'UAG':'Stop', 'CAG':'Q', 'AAG':'K', 'GAG':'E',
    'UGU':'C', 'CGU':'R', 'AGU':'S', 'GGU':'G',
    'UGC':'C', 'CGC':'R', 'AGC':'S', 'GGC':'G',
    'UGA':'Stop', 'CGA':'R', 'AGA':'R', 'GGA':'G',
    'UGG':'W', 'CGG':'R', 'AGG':'R', 'GGG':'G'
}

# Check the sequence to make sure it is a valid RNA string
def validateRNASeq(rna_seq):
    tmpseq = rna_seq.upper()
    for nuc in tmpseq:
        if nuc not in RNA_Nucleotides:
            return False
    return tmpseq

# Translate an RNA sequence into a protein string
def translate_rna(rna_seq):
    protein = ""
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]
        amino_acid = rnacodon_table.get(codon, '')
        if amino_acid == 'Stop':
            break
        protein += amino_acid
    return protein
