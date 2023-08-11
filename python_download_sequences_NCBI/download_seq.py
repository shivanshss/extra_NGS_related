#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: bloodmark
"""
import os
from Bio import Entrez, SeqIO
from urllib.error import HTTPError

# Set email address for NCBI Entrez
Entrez.email = "your_email_address"

# Path to the text file with accession numbers and names
input_file = "/test/accession_list.txt"

# Directory to save the downloaded sequences
output_directory = "downloaded.sequences"

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

def download_sequence(accession, name):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()

        output_file = os.path.join(output_directory, f"{name}_{accession}.fasta")
        with open(output_file, "w") as f:
            f.write(record)
        
        print(f"Downloaded and saved {accession} as {output_file}")
    except HTTPError as e:
        print(f"Error downloading {accession}: {e}")
    except Exception as e:
        print(f"An error occurred for {accession}: {e}")

def main():
    with open(input_file, "r") as f:
        for line in f:
            accession, name = line.strip().split("\t")
            download_sequence(accession, name)

if __name__ == "__main__":
    main()
