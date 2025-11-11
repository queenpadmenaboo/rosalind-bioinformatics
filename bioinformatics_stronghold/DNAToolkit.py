# DNA Toolkit file
"""
Author: Bunsree Patel
Date: October 24, 2025
"""

# Define nucleotide sets
DNA_Nucleotides = ['A', 'C', 'G', 'T']
RNA_Nucleotides = ['A', 'C', 'G', 'U']


# Check the sequence to make sure it is a valid DNA string
def validateDNASeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in DNA_Nucleotides:
            return False
    return tmpseq

def countDNANucFrequency(dna_seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna_seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def countRNANucFrequency(rna_seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "U": 0}
    for nuc in rna_seq:
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
    for i in range(0, len(rna_seq) - 2, 3):
        codon = rna_seq[i:i+3]
        amino_acid = rnacodon_table.get(codon, '')
        if amino_acid == 'Stop':
            break
        protein += amino_acid
    return protein


def find_substring_locations(s,t):
    """
    Find all locations of substring t in string s.
    Returns 1-indexed positions (biology convention).
    """
    locations = []
    for i in range(len(s) - len(t) + 1):
        if s[i:i+len(t)] == t:
            locations.append(i + 1)
    return locations

# Extract the first k nucleotides from a DNA sequence
def prefix(seq, k):
    return seq[:k]

# Extracts the last k nucleotides from a DNA sequence
def suffix(seq, k):
    return seq[-k:]

# Solves RNA Splicing - removes introns, transribes DNA to RNA, and translates to protein
def solve_splc(fasta_text):
    # Parse all sequences from FASTA format
    sequences = parse_fasta(fasta_text)
    seq_list = list(sequences.values())

    # First sequence is the gene, rest are introns to remove
    gene = seq_list[0]
    introns = seq_list[1:]

    # Remove all introns from the gene
    for intron in introns:
        gene = gene.replace(intron, "")

    # Validate the spliced gene sequence
    gene = validateDNASeq(gene)

    # Transcribe the DNA to RNA
    rna = transcribeDNAtoRNA(gene)

    # Translation RNA to protein string
    protein = translate_rna(rna)

    return protein


def transition_transversion_ratio(fasta_text):
    # Parse sequences and remove whitespace/newlines
    sequences = [seq.replace('\n', '').strip() for seq in parse_fasta(fasta_text).values()]
    s1, s2 = sequences[0], sequences[1]

    # Define transitions
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    # Define transversions
    transversions = {
        ('A', 'C'), ('A', 'T'), 
        ('C', 'A'), ('C', 'G'), 
        ('G', 'C'), ('G', 'T'), 
        ('T', 'A'), ('T', 'G')
    }
    
    transitions_count = 0
    transversions_count = 0

    # Compare psoitions
    for base1, base2 in zip(s1, s2):
        if base1 != base2:
            if (base1, base2) in transitions:
                transitions_count += 1
            elif (base1, base2) in transversions:
                transversions_count += 1 

    if transversions_count == 0:
            return 0.0
        
    return transitions_count / transversions_count
