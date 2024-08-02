# README
# This script filters sequences from a fasta file based on a list of sequence names provided in a separate file.
# The script reads the list of sequence names, matches them to the sequence identifiers in the fasta file,
# and writes the matching sequences to a new fasta file.

# Author: Bloodmark

# Input details:
# - Input fasta file: "/media/bloodmark/HDD6_SS_extra/transcriptome_assembly/JNS_RNAseq/SOAP_assembly_fasta_files/SOAP_31.fasta"
#   - The file should be in fasta format.
# - Names file: "est_seq_names"
#   - The file should contain a list of sequence names, one per line.

# Output details:
# - Output fasta file: "est_seq.fasta"
#   - The file will contain sequences from the input fasta file that match the names in the names file.

# Script functioning:
# 1. Install and load the Bioconductor package Biostrings if not already installed.
# 2. Define a function to filter the fasta file based on a list of sequence names.
# 3. Read the list of sequence names.
# 4. Read the fasta file and extract sequence identifiers.
# 5. Identify sequence identifiers in the fasta file that match the names in the list.
# 6. Read the full fasta sequences and identify sequences corresponding to the matching headers.
# 7. Write the matching sequences to the output fasta file.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install and load the Bioconductor package Biostrings if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}

library(Biostrings)

# Function to filter fasta file based on a list of sequence names
filter_fasta <- function(input_fasta, output_fasta, names_file) {
  # Read the list of names
  target_names <- tolower(trimws(readLines(names_file)))
  
  # Read the fasta file and extract sequence identifiers
  fasta_headers <- readLines(input_fasta)
  fasta_headers <- fasta_headers[grep("^>", fasta_headers)]
  fasta_headers <- tolower(trimws(sub("^>(\\S+).*", "\\1", fasta_headers)))
  
  # Display some information for debugging
  cat("First 10 names from the list:", head(target_names, 10), "\n")
  cat("First 10 sequence identifiers from the fasta file:", head(fasta_headers, 10), "\n")
  
  # Identify sequence identifiers in the fasta file that match the names in the list
  matching_indices <- which(fasta_headers %in% target_names)
  
  # Display some information for debugging
  cat("Number of matching sequence identifiers:", length(matching_indices), "\n")
  
  # Read the full fasta sequences
  fasta_sequences <- readDNAStringSet(input_fasta)
  
  # Identify sequences corresponding to the matching headers
  matching_sequences <- fasta_sequences[matching_indices]
  
  # Display some information for debugging
  cat("Number of matching sequences:", length(matching_sequences), "\n")
  
  # Write the matching sequences to the output fasta file
  writeXStringSet(matching_sequences, output_fasta, format = "fasta")
}

# Example usage
input_fasta_file <- "/media/bloodmark/HDD6_SS_extra/transcriptome_assembly/JNS_RNAseq/SOAP_assembly_fasta_files/SOAP_31.fasta"  # Replace with your input file name
output_fasta_file <- "est_seq.fasta"  # Replace with your output file name
names_file <- "est_seq_names"  # Replace with your list of names file name

filter_fasta(input_fasta_file, output_fasta_file, names_file)

