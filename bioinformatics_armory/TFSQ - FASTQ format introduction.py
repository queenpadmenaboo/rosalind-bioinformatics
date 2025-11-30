"""
Rosalind Problem: TFSQ - FASTQ format introduction
Problem URL: http://rosalind.info/problems/tfsq/
Location: Bioinformatics Armory
Author: Bunsree Patel
Date: November 29, 2025
"""

"""The Phred quality score Q is calculated using the following equation:
Q = -10logP
Q = -10 * math.log10(P)

where P is the probability that the corresponding base call is incorrect. 
Thus, for Q=10, the probability that the base is called incorrectly is 1 in 10, and the accuracy is 90%.
For Q=20, the probability of error is 1 in 100 and the accuracy is 99%, and so on.

'FASTQ format' - standard format for storing the output of high throughput sequencing instruments.
    - provides both the sequence and the per-base quality scores for each read
    - an extension of FASTA format, generally with suffix '.fq' or '.fastq'
    - Phred scores are usually 2-digit numbers i.e. 10, 50
    - Nucleotides are represented by a single symbol
    - When listed together, an offset of '33' is added to the score, which is then represented by the corresponding symbol from the ASCII table:      
        e.g. quality of a base with a Phred score of 10 (90% accuracy) is denoted by "+" (ASCII symbol #43)

    -Each record in a FASTQ file typically consists of four lines:
        - A line starting with @ that contains the sequence identifier
        - The actual sequence
        - A line starting with + containing an optional sequence identifier or comment
        - A line with quality scores encoded as ASCII symbols
        ***Note*** that the 2nd and 4th line should always have the same length.         


'ASCII' (American Standard Code for Information Interchange) - a character-encoding table where each character is represented by the number from 0 to 127.
    - used for computer storage
    - first 32 codes are reserved for control characters (non-printable information) i.e. backspace, vertical tab
    - the other 95 are printable characters i.e. letters, punctuation marks, etc."""

"""Sample Dataset Problem
Sometimes it's necessary to convert data from FASTQ format to FASTA format. 
For example, you may want to perform a BLAST search using reads in FASTQ format obtained from your brand new Illumina Genome Analyzer.

Links:

    - A FASTQ to FASTA converter can be accessed from the Sequence conversion website at https://sequenceconversion.bugaco.com/converter/biology/sequences/index.php.
    - A free GUI converter developed by BlastStation is available at https://www.blaststation.com/intl/en/q2a-pro.php for download or as an add-on to Google Chrome.
    - There is a FASTQ to FASTA converter in the Galaxy web platform at https://usegalaxy.org/root?tool_id=cshl_fastq_to_fasta. Note that you should register in the Galaxy and upload your file prior to using this tool.

    Given: FASTQ file

    Return: Corresponding FASTA records"""

"""Sample Dataset"""

# Creat fastq file
fastq_content = """@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!*((((***+))%%%++)(%%%%).1***-+*****))**55CCF>>>>>>CCCCCCC65
"""

# Save file to disk
file_path = "C:/Users/meeko/Desktop/TFSQ sample.fastq"
with open(file_path, "w") as f:
    f.write(fastq_content)

from bioblend.galaxy import GalaxyInstance

# Server info
GALAXY_URL = "https://usegalaxy.org"  # or your Galaxy instance
API_KEY = "7e36e36f43bd157d7489a3ed48406f84"
HISTORY_ID = "bbd44e69cb8906b5297928531850d7f9"

gi = GalaxyInstance(url=GALAXY_URL, key=API_KEY)

# Upload a file
gi.tools.upload_file("path/to/local/file.txt", HISTORY_ID)

print("Upload complete!")
