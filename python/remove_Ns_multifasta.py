# README
# This script processes a fasta file to remove 'N' characters from the sequences.
# The script performs the following steps:
# 1. Reads the input fasta file.
# 2. Removes 'N' characters from each sequence.
# 3. Writes the cleaned sequences to a new fasta file.

# Author: Bloodmark

# Input details:
# - Input fasta file: "/home/bloodmark/targets_tcas52.fasta"
#   - The file should contain sequences in fasta format.

# Output details:
# - Output fasta file: "/home/bloodmark/targets_tcas52_No_Ns.fasta"
#   - The file will contain the cleaned sequences without 'N' characters.

# Script functioning:
# 1. Reads the input fasta file.
# 2. Removes 'N' characters from each sequence.
# 3. Writes the cleaned sequences to a new fasta file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

from Bio import SeqIO

input_file = "/home/bloodmark/targets_tcas52.fasta"  # replace with your input file name
output_file = "/home/bloodmark/targets_tcas52_No_Ns.fasta"  # replace with your desired output file name

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        # Remove Ns from the sequence
        clean_seq = str(record.seq).replace("N", "")
        # Write the cleaned sequence to the output file
        record.seq = clean_seq
        SeqIO.write(record, outfile, "fasta")

print("Finished removing Ns from sequences.")

