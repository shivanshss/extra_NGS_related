# README
# This script reads two data files into data frames, merges them based on matching values in specific columns,
# sorts the merged data frame by a specified column, and prints the resulting data frame.
# The script assumes that the input files are properly formatted and located in the specified paths.

# Author: Bloodmark

# Input details:
# - Input file 1: "/media/bloodmark/HDD6_SS_extra/genome_assembly_figure_data/synteny/TCAS/5000/blocks_coords_cleaned1.txt"
#   - The file should be tab-delimited and contain six columns without a header row.
# - Input file 2: "/home/bloodmark/seqid"
#   - The file should be tab-delimited and contain at least one column named 'Seq_id'.

# Output details:
# - A merged data frame based on matching values in Col2 of file1 and Seq_id of file2.
# - The merged data frame sorted by Col1.
# - The script prints the first few rows of the sorted merged data frame.

# Script functioning:
# 1. Read the first file into a data frame.
# 2. Read the second file into another data frame.
# 3. Merge the two data frames based on matching values in Col2 of file1 and Seq_id of file2.
# 4. Sort the merged data frame by Col1.
# 5. Print the first few rows of the sorted merged data frame.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure that the input files are in the specified paths and properly formatted.

# Code starts below:

# Read the first file into a data frame
file1 <- tryCatch({
  read.table("/media/bloodmark/HDD6_SS_extra/genome_assembly_figure_data/synteny/TCAS/5000/blocks_coords_cleaned1.txt", header = FALSE, col.names = c("Col1", "Col2", "Col3", "Col4", "Col5", "Col6"), sep = "\t")
}, error = function(e) {
  message("Error reading file1: ", e)
  NULL
})

# Read the second file into another data frame
file2 <- tryCatch({
  read.table("/home/bloodmark/seqid", header = TRUE)
}, error = function(e) {
  message("Error reading file2: ", e)
  NULL
})

# Proceed only if both files are read successfully
if (!is.null(file1) && !is.null(file2)) {
  # Merge the two data frames based on matching values in Col2 of file1 and Seq_id of file2
  merged_data <- merge(file1, file2, by.x = "Col2", by.y = "Seq_id", all.x = TRUE)
  
  # Sort the merged data frame by Col1
  merged_data_sorted <- merged_data[order(merged_data$Col1), ]
  
  # Print the first few rows of the sorted merged data frame
  print(head(merged_data_sorted))
} else {
  message("File reading failed, merging skipped.")
}

