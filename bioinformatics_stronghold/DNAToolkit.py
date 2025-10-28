# DNA Toolkit file
"""
Author: Bunsree Patel
Date: October 24, 2025
"""

Nucleotides = ['A', 'C', 'G', 'T']


# Check the sequence to make sure it is a valid DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
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
        if nuc not in Nucleotides:
            raise ValueError("Invalid DNA sequence")
    
    # Replace T with U to make RNA
    rna_seq = dna_seq.upper().replace("T", "U")
    return rna_seq

def complementDNA(dna_seq):
    # Validate input
    for nuc in dna_seq.upper():
        if nuc not in Nucleotides:
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
    fasta_dict = {}
    current_id = None

    for line in fasta_str.strip().splitlines():
        if line.startswith(">"):
            current_id = line[1:].strip()
            fasta_dict[current_id] = ""
        else:
            fasta_dict[current_id] += line.strip()
            
    return fasta_dict            