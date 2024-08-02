# README
# This script reads a set of sequences from a multifasta file and filters them based on a list of sequence names.
# The script performs the following steps:
# 1. Installs and loads the necessary Bioconductor package (Biostrings).
# 2. Reads sequences from a multifasta file.
# 3. Reads a list of sequence names from a text file.
# 4. Filters the sequences to include only those whose names are in the list.
# 5. Writes the filtered sequences to a new multifasta file.

# Author: Bloodmark

# Input details:
# - Input multifasta file: "/media/bloodmark/HDD6_SS_extra/DiscoverY/contigs/contigs_01B.fasta"
# - Input text file with sequence names: "/media/bloodmark/HDD6_SS_extra/DiscoverY/Final_plots/01B_y_contig_names.txt"

# Output details:
# - Output multifasta file: "01B_y_contigs.fasta" containing the filtered sequences.

# Script functioning:
# 1. Install and load the Bioconductor Biostrings package if not already installed.
# 2. Read sequences from the input multifasta file.
# 3. Read the list of sequence names from the input text file.
# 4. Filter the sequences based on the list of sequence names.
# 5. Write the filtered sequences to a new multifasta file.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install and load the Bioconductor Biostrings package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")

library(Biostrings)

# Read sequences from the multifasta file
sequences <- readDNAStringSet('/media/bloodmark/HDD6_SS_extra/DiscoverY/contigs/contigs_01B.fasta')

# Read sequence names from the text file
sequence_names <- readLines('/media/bloodmark/HDD6_SS_extra/DiscoverY/Final_plots/01B_y_contig_names.txt')

# Subset the sequences to include only those whose names are in the list
selected_sequences <- sequences[names(sequences) %in% sequence_names]

# Write the selected sequences to a new multifasta file
writeXStringSet(selected_sequences, file = '01B_y_contigs.fasta')

