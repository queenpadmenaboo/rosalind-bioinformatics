from Bio import SeqIO
SeqIO.convert("test.fastq", "fastq", "test.fasta", "fasta")

# Output will be in 'test.fasta' file under same file structure as 'test.fastq' and 'convert.py'