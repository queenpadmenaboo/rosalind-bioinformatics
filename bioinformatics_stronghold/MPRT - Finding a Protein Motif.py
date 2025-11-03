"""
Rosalind Problem: MPRT - Finding a Protein Motif
Problem URL: http://rosalind.info/problems/mprt/

Author: Bunsree Patel
Date: November 3, 2025
"""

"""
A 'protein domain' is a structural and functional unit the corresponds to a sequence of amino acids that can fold and function independently.

- Each domain corresponds to a single specific function (e.g., DNA binding, catalysis).
- Some proteins have a single domain i.e. myoglobin and the Cytochrome comples; others are multifunctional with several domains.
- It is possible to create 'chimeric proteins' thorugh artificial fusion of domains.

Like species, proteins evolve into 'homologous' groups called 'protein families', which share similar domain sets and functions.

A' motif' is a short, conserved amino acid pattern within a domain, critical to its function - also known as blocks, signatures, or fingerprints.

Proteins discovered in laboratories worldwide are collected in public databases.
'UniProt' is a key public repository for protein data, providing annotations for function, domain structures, post-translational modifications, taxonomy, and literature references."""

"""Sample Dataset Problem

To allow for the presence of its varying forms, a protein motif is represented by a shorthand as follows: [XY] means "either X or Y" and {X} means "any amino acid except X." 
For example, the N-glycosylation motif is written as N{P}[ST]{P}.

You can see the complete description and features of a particular protein by its access ID "uniprot_id" in the UniProt database, by inserting the ID number into:

http://www.uniprot.org/uniprot/uniprot_id

Alternatively, you can obtain a protein sequence in 'FASTA format' by following:

http://www.uniprot.org/uniprot/uniprot_id.fasta

For example, the data for protein B5ZC00 can be found at http://www.uniprot.org/uniprot/B5ZC00.

    Given: At most 15 UniProt Protein Database access IDs.
    Return: For each protein possessing the N-glycosylation motif, output its given access ID followed by a list of locations in the protein string where the motif can be found.
"""


import re
import requests

from ProteinToolkit import clean_fasta, find_n_glycosylation, fetch_uniprot_fasta

# FASTA format always starts with '>' character, followed by an identifier and optional description

protein_A2Z669 = """>sp|A2Z669|CSPLT_ORYSI CASP-like protein 5A2 OS=Oryza sativa subsp. indica OX=39946 GN=OsI_33147 PE=3 SV=1
MRASRPVVHPVEAPPPAALAVAAAAVAVEAGVGAGGGAAAHGGENAQPRGVRMKDPPGAP
GTPGGLGLRLVQAFFAAAALAVMASTDDFPSVSAFCYLVAAAILQCLWSLSLAVVDIYAL
LVKRSLRNPQAVCIFTIGDGITGTLTLGAACASAGITVLIGNDLNICANNHCASFETATA
MAFISWFALAPSCVLNFWSMASR"""

# sp → Swiss-Prot section of UniProt
# A2Z669 → UniProt accession number
# CSPLT_ORYSI → protein short name / identifier
# CASP-like protein 5A2 → protein name
# OS=Oryza sativa subsp. indica → organism name (rice, indica subspecies)
# OX=39946 → taxonomy ID
# GN=OsI_33147 → gene name
# PE=3 → protein existence evidence level (3 = inferred from homology)
# SV=1 → sequence version
# Each line of the sequence contains amino acid residues in single-letter codes.
# Line breaks are allowed anywhere — tools read it as one continuous sequence.


protein_B5ZC00 = """>sp|B5ZC00|SYG_UREU1 Glycine--tRNA ligase OS=Ureaplasma urealyticum serovar 10 (strain ATCC 33699 / Western) OX=565575 GN=glyQS PE=3 SV=1
MKNKFKTQEELVNHLKTVGFVFANSEIYNGLANAWDYGPLGVLLKNNLKNLWWKEFVTKQ
KDVVGLDSAIILNPLVWKASGHLDNFSDPLIDCKNCKARYRADKLIESFDENIHIAENSS
NEEFAKVLNDYEISCPTCKQFNWTEIRHFNLMFKTYQGVIEDAKNVVYLRPETAQGIFVN
FKNVQRSMRLHLPFGIAQIGKSFRNEITPGNFIFRTREFEQMEIEFFLKEESAYDIFDKY
LNQIENWLVSACGLSLNNLRKHEHPKEELSHYSKKTIDFEYNFLHGFSELYGIAYRTNYD
LSVHMNLSKKDLTYFDEQTKEKYVPHVIEPSVGVERLLYAILTEATFIEKLENDDERILM
DLKYDLAPYKIAVMPLVNKLKDKAEEIYGKILDLNISATFDNSGSIGKRYRRQDAIGTIY
CLTIDFDSLDDQQDPSFTIRERNSMAQKRIKLSELPLYLNQKAHEDFQRQCQK"""

protein_P07204 = """>sp|P07204|TRBM_HUMAN Thrombomodulin OS=Homo sapiens OX=9606 GN=THBD PE=1 SV=2
MLGVLVLGALALAGLGFPAPAEPQPGGSQCVEHDCFALYPGPATFLNASQICDGLRGHLM
TVRSSVAADVISLLLNGDGGVGRRRLWIGLQLPPGCGDPKRLGPLRGFQWVTGDNNTSYS
RWARLDLNGAPLCGPLCVAVSAAEATVPSEPIWEEQQCEVKADGFLCEFHFPATCRPLAV
EPGAAAAAVSITYGTPFAARGADFQALPVGSSAAVAPLGLQLMCTAPPGAVQGHWAREAP
GAWDCSVENGGCEHACNAIPGAPRCQCPAGAALQADGRSCTASATQSCNDLCEHFCVPNP
DQPGSYSCMCETGYRLAADQHRCEDVDDCILEPSPCPQRCVNTQGGFECHCYPNYDLVDG
ECVEPVDPCFRANCEYQCQPLNQTSYLCVCAEGFAPIPHEPHRCQMFCNQTACPADCDPN
TQASCECPEGYILDDGFICTDIDECENGGFCSGVCHNLPGTFECICGPDSALARHIGTDC
DSGKVDGGDSGSGEPPPSPTPGSTLTPPAVGLVHSGLLIGISIASLCLVVALLALLCHLR
KKQGAARAKMEYKCAAPSKEVVLQHVRTERTPQRL"""

protein_P20840 = """>sp|P20840|SAG1_YEAST Alpha-agglutinin OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) OX=559292 GN=SAG1 PE=1 SV=2
MFTFLKIILWLFSLALASAININDITFSNLEITPLTANKQPDQGWTATFDFSIADASSIR
EGDEFTLSMPHVYRIKLLNSSQTATISLADGTEAFKCYVSQQAAYLYENTTFTCTAQNDL
SSYNTIDGSITFSLNFSDGGSSYEYELENAKFFKSGPMLVKLGNQMSDVVNFDPAAFTEN
VFHSGRSTGYGSFESYHLGMYCPNGYFLGGTEKIDYDSSNNNVDLDCSSVQVYSSNDFND
WWFPQSYNDTNADVTCFGSNLWITLDEKLYDGEMLWVNALQSLPANVNTIDHALEFQYTC
LDTIANTTYATQFSTTREFIVYQGRNLGTASAKSSFISTTTTDLTSINTSAYSTGSISTV
ETGNRTTSEVISHVVTTSTKLSPTATTSLTIAQTSIYSTDSNITVGTDIHTTSEVISDVE
TISRETASTVVAAPTSTTGWTGAMNTYISQFTSSSFATINSTPIISSSAVFETSDASIVN
VHTENITNTAAVPSEEPTFVNATRNSLNSFCSSKQPSSPSSYTSSPLVSSLSVSKTLLST
SFTPSVPTSNTYIKTKNTGYFEHTALTTSSVGLNSFSETAVSSQGTKIDTFLVSSLIAYP
SSASGSQLSGIQQNFTSTSLMISTYEGKASIFFSAELGSIIFLLLSYLLF"""

# Test cleaning function
cleaned_seq = clean_fasta(protein_A2Z669)
print("Testing clean_fasta:")
print(cleaned_seq)

# Put proteins in a dictionary
proteins = { 
    "A2Z669": protein_A2Z669,
    "B5ZC00": protein_B5ZC00,
    "P07204_TRBM_HUMAN": protein_P07204,
    "P20840_SAG1_YEAST": protein_P20840
}

# Loop through all proteins and find N-glycosylation motif positions

print("Finding N-glycosylation motifs:")
for uid, fasta_text in proteins.items():
    seq = clean_fasta(fasta_text)
    positions = find_n_glycosylation(seq)
    if positions:
        print(uid)
        print(*positions)
    # 'uid' is the protein's unique identifier like 'A2Z669'
    # 'fasta_text' is the raw FASTA string for that protein
    # 'seq' now contains full protein sequence as a single string
    # 'positions' uses regular expressions to search for the motif pattern in 'seq'
    # 're.finditer' returns all occurences of 'motif' in 'seq'
    # 'm.start()' gives the starting index (0-based), +1 converts to 1-based indexing
    # 'if positions' only prints results if the motif is found at least once
    # 'print(uid)' - prints the protein ID
    # 'print(*positions) - prints all starting positions of the motif, separated by spaces


"""Actual Dataset Problem"""
import re
import requests
from ProteinToolkit import clean_fasta, find_n_glycosylation, fetch_uniprot_fasta

uniprot_ids = """
B0RU89
Q5U1Y9
P01044_KNH1_BOVIN
P04921_GLPC_HUMAN
B2G8U6
P04441_HG2A_MOUSE
P20840_SAG1_YEAST
Q66GC7
P07359_GPBA_HUMAN
Q05557
P01880_DTC_HUMAN
P00740_FA9_HUMAN
"""

uniprot_ids = [uid.strip() for uid in uniprot_ids.strip().split("\n") if uid.strip()]

for uid in uniprot_ids:
    accession = uid.split("_")[0]
    fasta_text = fetch_uniprot_fasta(accession)
    seq = clean_fasta(fasta_text)
    positions = find_n_glycosylation(seq)
    if positions:
        print(uid)
        print(*positions)