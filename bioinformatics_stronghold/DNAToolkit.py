# DNA Toolkit file

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