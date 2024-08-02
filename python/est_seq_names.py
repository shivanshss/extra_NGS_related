#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:17:59 2023

@author: bloodmark
"""

# README
# This script performs BLAST searches for nucleotide sequences in one file against the protein sequences in another file.
# It finds the best match for each nucleotide sequence and writes the results to an output file.
# The script performs the following steps:
# 1. Loads protein sequences from the first file.
# 2. Loads nucleotide sequences from the second file.
# 3. Performs BLAST searches for each nucleotide sequence against the NCBI 'nt' database.
# 4. Extracts the best match for each nucleotide sequence.
# 5. Writes the results to an output file.

# Author: Bloodmark

# Input details:
# - Input protein sequences file: "/media/bloodmark/HDD6_SS_extra/est_gene_name/insect_proteins.fasta"
#   - The file should contain protein sequences in fasta format.
# - Input nucleotide sequences file: "/media/bloodmark/HDD6_SS_extra/est_gene_name/est_seq.fasta"
#   - The file should contain nucleotide sequences in fasta format.

# Output details:
# - Output file: "output.txt"
#   - The file will contain columns: Name_in_file2, Best_match_in_file1.

# Script functioning:
# 1. Loads protein sequences from the first file.
# 2. Loads nucleotide sequences from the second file.
# 3. Performs BLAST searches for each nucleotide sequence against the NCBI 'nt' database.
# 4. Extracts the best match for each nucleotide sequence.
# 5. Writes the results to an output file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input files are in the specified paths and properly formatted.
# Note: Ensure you have internet access to perform BLAST searches.

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

def find_best_match(query_sequence):
    # Perform BLAST search
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence.seq)

    # Parse the BLAST result
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    # Check if there are any alignments
    if blast_record.alignments:
        # Extract the best match
        best_alignment = blast_record.alignments[0]

        # Extract the name of the best match
        best_match_name = best_alignment.hit_def
    else:
        # No alignments found
        best_match_name = "No Match"

    return best_match_name

# Load protein sequences from file1
file1_sequences = list(SeqIO.parse("/media/bloodmark/HDD6_SS_extra/est_gene_name/insect_proteins.fasta", "fasta"))

# Load nucleotide sequences from file2
file2_sequences = list(SeqIO.parse("/media/bloodmark/HDD6_SS_extra/est_gene_name/est_seq.fasta", "fasta"))

# Open a file for writing the output
with open("output.txt", "w") as output_file:
    # Write the header
    output_file.write("Name_in_file2\tBest_match_in_file1\n")

    # Iterate over nucleotide sequences in file2
    for file2_sequence in file2_sequences:
        # Find the best match in the NCBI 'nt' database
        best_match_name = find_best_match(file2_sequence)

        # Write the results to the file
        output_file.write(f"{file2_sequence.id}\t{best_match_name}\n")

print("Output written to output.txt")

