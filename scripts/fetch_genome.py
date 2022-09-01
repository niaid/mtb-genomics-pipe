
from Bio import Entrez
import argparse

"""
fetches reference genome
"""

Entrez.email = "brendan.jeffrey@nih.gov"

# cl args
parser = argparse.ArgumentParser(description='Pull reference genome using Entrez')
parser.add_argument('-i', dest='accession', required=True, help='reference genome accession id')
parser.add_argument('-o', dest='fasta_out', required=True, help='reference genome fasta')

args = parser.parse_args()
accession = args.accession
fasta_out = args.fasta_out

record_fasta = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")

with open(fasta_out, "w") as fasta:
    fasta.write(record_fasta.read())
    
record_fasta.close()
