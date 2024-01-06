#!/bin/bash

# Define the input file containing three columns
input_file="/home/bloodmark/workarea/synteny/synteny_filepath.txt"

# Loop through each line in the input file
while IFS=$'\t' read -r folder_name fasta_file1 fasta_file2; do
    # Create directories
    mkdir -p "/home/bloodmark/workarea/synteny/$folder_name"

    # Change directory
    cd "/home/bloodmark/workarea/synteny/$folder_name"

    # Run sibeliaz
    sibeliaz "$fasta_file1" "$fasta_file2"

    # Run maf2synteny
    maf2synteny sibeliaz_out/alignment.maf

    # Move back to the original directory if needed
    cd -

done < "$input_file"
