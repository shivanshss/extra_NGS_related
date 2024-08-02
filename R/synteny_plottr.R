# README
# This script reads fasta files, displays sequence names and lengths, and creates a synteny plot using the syntenyPlotteR package.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Reads the fasta files and displays sequence information.
# 3. Creates a synteny plot using the syntenyPlotteR package.

# Author: Bloodmark

# Input details:
# - Fasta file 1: "/home/bloodmark/workarea/new_final_chr/NC_007417.3_RagTag.fa"
# - Fasta file 2: "/home/bloodmark/workarea/pubref/pubref_chromosomes/LG2.fa"
# - Data file: "/home/bloodmark/workarea/synteny/LG2/5000/lg2_blocks_coords_syntenyplottr.txt"

# Output details:
# - Displays sequence names and lengths for each fasta file.
# - A synteny plot created using syntenyPlotteR.

# Script functioning:
# 1. Install and load the Bioconductor package Biostrings and the syntenyPlotteR package.
# 2. Define a function to read fasta files and display sequence names and lengths.
# 3. Read and display information for the fasta files.
# 4. Prepare data for the synteny plot and call the draw.eh function.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install and load the required Bioconductor package
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("Biostrings")

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("Farre-lab/syntenyPlotteR")

library(Biostrings)
library(syntenyPlotteR)

# Function to read fasta file and display sequence names and lengths
display_fasta_info <- function(fasta_file) {
  # Read the fasta file
  sequences <- readDNAStringSet(fasta_file)
  
  # Display sequence names and lengths
  cat("Sequence Name\tSequence Length\n")
  for (i in seq_along(names(sequences))) {
    seq_name <- names(sequences)[i]
    seq_length <- width(sequences[i])
    cat(seq_name, "\t", seq_length, "\n")
  }
}

# Fasta file paths (replace with your actual file paths)
fasta_file1 <- "/home/bloodmark/workarea/new_final_chr/NC_007417.3_RagTag.fa"
fasta_file2 <- "/home/bloodmark/workarea/pubref/pubref_chromosomes/LG2.fa"
file_data <- read.table(file = "/home/bloodmark/workarea/synteny/LG2/5000/lg2_blocks_coords_syntenyplottr.txt", header = FALSE)

# Display information for the first fasta file
cat("Information for", fasta_file1, ":\n")
display_fasta_info(fasta_file1)

# Display information for the second fasta file
cat("\nInformation for", fasta_file2, ":\n")
display_fasta_info(fasta_file2)

chr_id <- c("NC_007417.3_RagTag", "NC_007417.3")
chr_len <- c("14807813", "15265516")
sp_id <- c("lg2_new", "lg2_old")

data_file <- data.frame(chr_id, chr_len, sp_id)

# Write the data frame to a temporary file
temp_file <- tempfile(fileext = ".txt")
write.table(data_file, file = temp_file, sep = "\t", row.names = FALSE, col.names = FALSE)

# Call the draw.eh function with the temporary file
draw.eh("output", temp_file, file_data, fileformat = "png", colour = "lightblue", inverted.colour = "lightpink", w = 5.5, h = 10, ps = 10)

# Remove the temporary file after using it
unlink(temp_file)

