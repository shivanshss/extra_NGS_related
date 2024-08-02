#!/bin/bash

# README
# This script processes a multifasta file and splits each sequence into a separate fasta file.
# Each sequence will be written to a new file named after the sequence header (without the '>' character and any additional descriptions).
# The script requires one input argument which is the multifasta file to be processed.

# Author: Bloodmark

# Input details:
# - Multifasta file provided as a command line argument.

# Output details:
# - Individual fasta files for each sequence in the multifasta file. The files are named after the sequence headers.

# Script functioning:
# 1. Ensure a multifasta file is provided as an argument.
# 2. Check if the input file exists.
# 3. Read the multifasta file line by line.
# 4. For each sequence header (lines starting with '>'), create a new output file named after the header.
# 5. For each sequence line, append it to the current output file.

# Example of running the script:
# bash split_multifasta.sh <multifasta_file>

# Note: Ensure the input multifasta file is in the correct format.

# Code starts below:

# Ensure a file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <multifasta_file>"
    exit 1
fi

input_file="$1"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found!"
    exit 1
fi

# Read the multifasta file and split each sequence into a separate file
while read line
do
    if [[ $line == '>'* ]]; then
        # Extract the sequence name without '>'
        seq_name=$(echo $line | cut -d '>' -f 2 | cut -d ' ' -f 1)
        # Create or overwrite a new file
        output_file="${seq_name}.fasta"
        echo $line > "$output_file"
    else
        # Append the sequence to the current file
        echo $line >> "$output_file"
    fi
done < "$input_file"

