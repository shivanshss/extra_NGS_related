# README
# This script processes all FASTA files in a specified directory to count the total length of sequences and the number of 'N' characters in each sequence.
# The script performs the following steps:
# 1. Defines functions to count 'N' characters in sequences and read FASTA files.
# 2. Reads each FASTA file in the specified directory.
# 3. Calculates the total length of sequences and the count of 'N' characters.
# 4. Outputs the results for each file.

# Author: Bloodmark

# Input details:
# - Directory containing FASTA files: "/home/bloodmark/workarea/new_final_chr/chromosomes/"
#   - The directory should contain files with the extension ".fasta".

# Output details:
# - A data frame printed to the console with the following columns:
#   - File: The name of the FASTA file (without path and extension).
#   - TotalLength: The total length of the sequences in the file.
#   - NCount: The number of 'N' characters in the sequences.

# Script functioning:
# 1. Define a function to count the number of 'N' characters in a sequence.
# 2. Define a function to read a FASTA file and count 'N's in each sequence.
# 3. Define a function to process all FASTA files in a directory.
# 4. Process the specified directory and output the results.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified directory and properly formatted.

# Code starts below:

# Function to count the number of 'N' characters in a sequence
count_Ns <- function(sequence) {
  n_count <- sum(strsplit(sequence, NULL)[[1]] == 'N')
  return(n_count)
}

# Function to read a FASTA file and count 'N's in each sequence
count_Ns_in_fasta <- function(file_path) {
  # Read the FASTA file
  fasta_data <- readLines(file_path)
  
  # Initialize variables to store sequence information
  current_sequence <- NULL
  total_length <- 0
  
  # Iterate through the lines in the FASTA file
  for (line in fasta_data) {
    if (startsWith(line, ">")) {
      # If line starts with '>', it's a header line
      if (!is.null(current_sequence)) {
        # If there was a previous sequence, update total length
        total_length <- total_length + nchar(current_sequence)
      }
      
      # Reset current sequence for the new header
      current_sequence <- ""
    } else {
      # If the line is not a header, it's part of the sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Update total length for the last sequence in the file
  total_length <- total_length + nchar(current_sequence)
  
  # Count 'N's in the entire sequence
  n_count <- count_Ns(current_sequence)
  
  return(list(total_length = total_length, n_count = n_count))
}

# Function to process all FASTA files in a directory
process_fasta_directory <- function(directory_path) {
  # Get a list of all files in the directory
  fasta_files <- list.files(directory_path, pattern = "\\.fasta$", full.names = TRUE)
  
  # Initialize a data frame to store results
  results <- data.frame(File = character(), TotalLength = numeric(), NCount = numeric(), stringsAsFactors = FALSE)
  
  # Iterate through each FASTA file
  for (file_path in fasta_files) {
    # Get results for the current file
    file_results <- count_Ns_in_fasta(file_path)
    
    # Extract file name without the path and extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Add results to the data frame
    results <- rbind(results, c(file_name, file_results$total_length, file_results$n_count))
  }
  
  return(results)
}

# Path to directory containing FASTA files
fasta_directory_path <- "/home/bloodmark/workarea/new_final_chr/chromosomes/"
result <- process_fasta_directory(fasta_directory_path)

# Print the final results
print(result)

