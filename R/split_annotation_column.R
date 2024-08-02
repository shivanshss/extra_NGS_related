# README
# This script processes an Excel file and retrieves gene names from NCBI based on accession numbers.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Reads the data from an Excel file and performs data cleaning.
# 3. Splits the content of a specific column.
# 4. Writes the processed data to a tab-delimited file.
# 5. Filters the accession numbers and writes them to separate files.
# 6. Queries NCBI for gene names based on the accession numbers and writes the results to a file.
# 7. Plays a sound notification upon completion.

# Author: Bloodmark

# Input details:
# - Input Excel file: "/media/bloodmark/HDD6_SS_extra/final_gff/chr3.tab"
# - Input text file with accession numbers: "/home/bloodmark/workarea/20231011_maker/chr4_gene_codes.tab"

# Output details:
# - Processed data written to "chr3_split.tab".
# - Filtered accession numbers written to "chr4_gi_accessions.txt" and "chr4_non_gi_accessions.txt".
# - Gene names written to "chr3_gene_names.txt".

# Script functioning:
# 1. Install and load the necessary packages (readxl, rentrez, beepr, Biostrings).
# 2. Read the data from the Excel file and clean the data.
# 3. Split the content of a specific column and convert the data to a data frame.
# 4. Convert columns to numeric and calculate a new column.
# 5. Write the combined data to a tab-delimited file.
# 6. Filter the accession numbers and write them to separate files.
# 7. Query NCBI for gene names and write the results to a file.
# 8. Play a sound notification upon completion.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("Biostrings")
install.packages("readxl")
install.packages("rentrez")
install.packages("beepr")

# Load the required packages
library(readxl)
library(rentrez)
library(beepr)

# Load the data from the Excel file
data <- read.table(file = "/media/bloodmark/HDD6_SS_extra/final_gff/chr3.tab", header = TRUE, col.names = c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6", "Col7", "Col8", "Col9", "Col10", "Col11", "Col12"), sep = "\t")

# Replace blank spaces with "NA" in the data
null_cells <- is.na(data)
data[null_cells] <- NA

# Split the content of each row at every ';'
split_data <- strsplit(data$Col9, ";")

# Find the maximum number of columns generated
max_columns <- max(sapply(split_data, length))

# Create a matrix to store the split values, filling empty elements with NA
split_data_matrix <- matrix("", nrow = length(split_data), ncol = max_columns)

for (i in 1:length(split_data)) {
  split_data_matrix[i, 1:length(split_data[[i]])] <- split_data[[i]]
}

# Convert the matrix to a data frame
split_data_df <- as.data.frame(split_data_matrix)

# Set column names
colnames(split_data_df) <- paste("Col", 1:max_columns, sep = "")

# Convert Col5 and Col4 to numeric if they are stored as characters
non_numeric_rows <- which(is.na(as.numeric(data$Col5)))
cat("Rows with non-numeric values in Col5: ", non_numeric_rows, "\n")

# Print the values in Col5 for the identified rows
cat("Non-numeric values in Col5: ", data$Col5[non_numeric_rows], "\n")

data$Col5 <- as.numeric(data$Col5)
data$Col4 <- as.numeric(data$Col4)
data$Col6 <- data$Col5 - data$Col4

# Create a data frame with selected columns from data and split_data_df
combined_data <- data.frame(
  Col1 = data$Col1,
  Col2 = data$Col2,
  Col4 = data$Col4,
  Col5 = data$Col5,
  Col6 = data$Col6,
  split_data_df$Col2,
  split_data_df$Col3,
  split_data_df$Col4,
  split_data_df$Col5
)

# Specify the file path
output_file <- "chr3_split.tab"

# Write the combined data frame to a tab-delimited file
write.table(combined_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Process accession numbers
file <- readLines("/home/bloodmark/workarea/20231011_maker/chr4_gene_codes.tab")
accession_list <- gsub("^Name=", "", file)
writeLines(accession_list, "chr4_gene_codes.tab")

# Filter rows starting with "gi"
gi_rows <- grep("^gi", file, value = TRUE)
writeLines(gi_rows, "chr4_gi_accessions.txt")

# Filter rows that do not start with "gi"
non_gi_rows <- grep("^(?!gi)", file, perl = TRUE, value = TRUE)
writeLines(non_gi_rows, "chr4_non_gi_accessions.txt")

# Define your list of NCBI accession codes
accession_list <- split_data_df$Col2

# Remove the "Name=" prefix from the accession codes
accession_list <- gsub("^Name=", "", accession_list)

# Create an empty list to store gene names
gene_names <- list()

# Loop through the accession list
for (accession in accession_list) {
  tryCatch({
    # Define a list of databases to query (you can add more if needed)
    databases <- c("protein", "nucleotide")
    
    # Initialize a variable to store the result
    gene_name <- NULL
    
    # Loop through the databases
    for (db in databases) {
      gene_names_chunk <- entrez_summary(db = db, id = accession, parsed = TRUE)
      
      if (!is.null(gene_names_chunk$gene)) {
        gene_name <- gene_names_chunk$gene
        break  # Exit the loop once a gene name is found
      }
    }
    
    if (!is.null(gene_name)) {
      gene_names[[accession]] <- gene_name
    } else {
      gene_names[[accession]] <- "NA"
    }
  }, warning = function(w) {
    cat("Warning:", conditionMessage(w), "\n")
  }, error = function(e) {
    cat("Error:", conditionMessage(e), "\n")
  })
}

# Open a file for writing (change "output.txt" to your desired file name)
output_file <- file("chr3_gene_names.txt", "w")

for (accession in names(gene_names)) {
  # Create a line of text to write
  line <- paste(accession, ": ", gene_names[[accession]], "\n")
  
  # Write the line to the file
  writeLines(line, output_file)
}

# Close the file
close(output_file)

# Play a sound notification upon completion
beep()

