# README
# This script processes multiple fasta files, creates a multiple sequence alignment, generates a consensus sequence, and calculates summary statistics.
# The script performs the following steps:
# 1. Installs and loads the necessary Bioconductor packages (Biostrings and DECIPHER).
# 2. Reads sequences from multiple fasta files.
# 3. Creates a multiple sequence alignment from the sequences.
# 4. Generates a consensus sequence based on a threshold.
# 5. Writes the consensus sequence to a fasta file.
# 6. Calculates and prints summary statistics (number of base pairs, number of sequences, and base percentages).

# Author: Bloodmark

# Input details:
# - List of multifasta files:
#   - "/home/bloodmark/12B_y_contigs.fasta"
#   - "/home/bloodmark/12E_y_contigs.fasta"
#   - "/home/bloodmark/13C_y_contigs.fasta"
#   - "/home/bloodmark/13E_y_contigs.fasta"
#   - "/home/bloodmark/18B_y_contigs.fasta"
#   - "/home/bloodmark/18A_y_contigs.fasta"

# Output details:
# - Consensus sequence written to "chrY_consensus.fasta".
# - Summary statistics printed to the console.

# Script functioning:
# 1. Install and load the Bioconductor Biostrings and DECIPHER packages if not already installed.
# 2. Read sequences from the multifasta files.
# 3. Create a multiple sequence alignment.
# 4. Generate a consensus sequence from the alignment based on a threshold.
# 5. Write the consensus sequence to a fasta file.
# 6. Calculate and print summary statistics (number of base pairs, number of sequences, and base percentages).

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

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

library(DECIPHER)

# List of your multifasta files (replace with your file paths)
fasta_files <- c("/home/bloodmark/12B_y_contigs.fasta", "/home/bloodmark/12E_y_contigs.fasta", "/home/bloodmark/13C_y_contigs.fasta", "/home/bloodmark/13E_y_contigs.fasta", "/home/bloodmark/18B_y_contigs.fasta", "/home/bloodmark/18A_y_contigs.fasta")

# Read sequences from the multifasta files
seq_list <- readDNAStringSet(fasta_files)

# Create a multiple sequence alignment
alignment <- AlignSeqs(seq_list)

# Generate a consensus sequence from the alignment
consensus <- ConsensusSequence(alignment)

# Set the threshold (at least 3 samples must have a consensus at a position)
threshold <- 3

# Generate the consensus sequence based on the threshold
final_consensus <- consensus >= threshold

# Print or save the final consensus sequence
cat(as.character(final_consensus), "\n")

# Write the consensus sequence to a fasta file
writeXStringSet(DNAStringSet(final_consensus), "chrY_consensus.fasta")

# Calculate summary statistics
num_base_pairs <- length(final_consensus)
num_sequences <- length(seq_list)
base_percentages <- table(as.character(final_consensus)) / num_base_pairs * 100

# Print or save the summary statistics
cat("Number of Base Pairs:", num_base_pairs, "\n")
cat("Number of Sequences:", num_sequences, "\n")
cat("Base Percentages:\n")
print(base_percentages)

