# README
# This script generates random DNA sequences, performs a multiple sequence alignment, 
# and visualizes the alignment using various R packages including Biostrings, msa, and ggmsa.
# The script first installs and loads the necessary libraries, generates random sequences,
# performs the alignment, and then plots the alignment.

# Author: Bloodmark

# Input details:
# - The script generates 25 random DNA sequences each of 1000 bases long. No external input file is required.

# Output details:
# - The script produces a plot of the multiple sequence alignment.

# Script functioning:
# 1. Install and load necessary libraries.
# 2. Define a function to generate random DNA sequences.
# 3. Generate 25 random DNA sequences of 1000 bases each.
# 4. Create a DNAStringSet object from the generated sequences.
# 5. Perform multiple sequence alignment using the msa package.
# 6. Convert the alignment to a data frame for visualization.
# 7. Plot the alignment using ggmsa with specified parameters.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure internet connectivity is available for installing the necessary packages.

# Install necessary packages
install.packages("Biostrings")
install.packages("msa")
install.packages("ggplot2")
install.packages("ggmsa")

# Load libraries
library(Biostrings)
library(msa)
library(ggmsa)

# Generate random sequences
generate_random_sequence <- function(length) {
  paste(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
}

sequences <- sapply(1:25, function(x) generate_random_sequence(1000))
names(sequences) <- paste0("Seq", 1:25)

# Create DNAStringSet
dna_set <- DNAStringSet(sequences)

# Perform multiple sequence alignment
alignment <- msa(dna_set)

# Convert alignment to data frame for visualization
alignment_df <- as(alignment, "data.frame")

# Plot the alignment
ggmsa(alignment_df, start = 1, end = 1000, color = "Chemistry_NT") +
  ggtitle("Multiple Sequence Alignment")

