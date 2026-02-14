import sys
import os

FILENAME = "felzartamab.fasta"

def read_fasta(filename):
    """Reads a FASTA file and returns the sequence as a single string."""
    if not os.path.exists(filename):
        print(f"Error: {filename} not found.")
        return ""
    
    sequence = ""
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue  # Skip headers and empty lines
            sequence += line
    return sequence

def get_stacked_positions(sequence):

    row_length = 20
    # Iterate through the sequence in blocks of 20
    for start in range(0, len(sequence), row_length):
        end = start + row_length
        block = sequence[start:end]

        # Header Row: Position numbers (01, 02, ... 20)
        positions = [f"{i+1:02}" for i in range(start, min(end, len(sequence)))]
        print(f"Pos: {'  '.join(positions)}")
        print("-" * (len(positions) * 3 + 5))
        
        # Data Row: Amino Acids separated by double spaces
        row_data = "   ".join(block)
        print(f"S01: {row_data}\n")

if __name__ == "__main__":
    # 1. Extract sequence from file
    sequence_data = read_fasta(FILENAME)
    
    # 2. Pass it into the function
    if sequence_data:
        get_stacked_positions(sequence_data)
    else:
        print("Sequence is empty or file missing.")