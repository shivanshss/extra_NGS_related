#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 18:58:45 2023

@author: bloodmark
"""

# README
# This script reads a multi-fasta file, concatenates the sequences, and writes the merged sequence to a new fasta file with a custom sequence name.
# The script performs the following steps:
# 1. Reads the multi-fasta file line by line.
# 2. Concatenates the sequence lines, skipping header lines.
# 3. Writes the merged sequence to a new fasta file with a custom sequence name.

# Author: Bloodmark

# Input details:
# - Input multi-fasta file: "/media/bloodmark/HDD6_SS_extra/new_final_chr/20231025/chr/unplaced_contigs.fasta"
#   - The file should contain multiple sequences in fasta format.

# Output details:
# - Output fasta file: "/media/bloodmark/HDD6_SS_extra/new_final_chr/20231025/chr/unplaced_contigs_scaff.fasta"
#   - The file will contain the merged sequence with the custom name "unplaced_contigs".

# Script functioning:
# 1. Reads the multi-fasta file line by line.
# 2. Concatenates the sequence lines, skipping header lines.
# 3. Writes the merged sequence to a new fasta file with a custom sequence name.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

# Define the input and output file names
input_file = "/media/bloodmark/HDD6_SS_extra/new_final_chr/20231025/chr/unplaced_contigs.fasta"
output_file = "/media/bloodmark/HDD6_SS_extra/new_final_chr/20231025/chr/unplaced_contigs_scaff.fasta"

merged_sequence_name = "unplaced_contigs"

# Initialize the merged sequence
merged_sequence = ""

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith(">"):
            # Skip header lines
            continue
        # Concatenate the sequence lines
        merged_sequence += line.strip()

    # Write the merged sequence with the custom name
    outfile.write(f">{merged_sequence_name}\n{merged_sequence}\n")

