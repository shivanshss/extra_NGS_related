# README
# This script processes GFF files by removing overlapping annotations, merging attributes, and prioritizing annotations with "Liftoff" in column 2.
# The script performs the following steps:
# 1. Installs and loads necessary packages.
# 2. Reads the GFF file and removes overlapping annotations.
# 3. Reads and processes another GFF file, removes duplicate rows, and concatenates the two files.
# 4. Removes overlapping annotations with a priority for "Liftoff" annotations.
# 5. Writes the processed data to new GFF files.

# Author: Bloodmark

# Input details:
# - Input GFF file 1: "/home/bloodmark/workarea/final_maker/maker_output_2023/round3_gff/round1_all_corrected3.gff"
#   - The file should be in GFF format.
# - Input GFF file 2: "/home/bloodmark/workarea/new_final_chr/annotation/liftoff_maker_newref/maker_output_loc_corr_tsv.gff"
#   - The file should be in GFF format.

# Output details:
# - Output GFF file 1: "maker_overlap_removed.gff"
#   - The file will contain annotations from the first input GFF file with overlapping annotations removed.
# - Output GFF file 2: "maker_overlap_removed_liftoff_correct.gff"
#   - The file will contain concatenated annotations from both input GFF files with overlapping annotations removed and prioritized.
# - Output GFF file 3: "maker_annotation_col_corr.gff"
#   - The file will contain the final annotations with an additional column extracted from the attributes.

# Script functioning:
# 1. Install and load the necessary packages.
# 2. Read the first GFF file, remove overlapping annotations, and write the result to a new GFF file.
# 3. Read the second GFF file, filter rows, remove duplicates, and concatenate with the first GFF file.
# 4. Remove overlapping annotations with a priority for "Liftoff" annotations.
# 5. Write the processed data to new GFF files.
# 6. Extract "Name=" information from the attributes and write to a new GFF file.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install and load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

if (!requireNamespace("progress", quietly = TRUE)) {
  install.packages("progress")
}
library(progress)

# Read the GFF file into a data frame
gff <- read.table("/home/bloodmark/workarea/final_maker/maker_output_2023/round3_gff/round1_all_corrected3.gff", sep="\t", comment.char="#", header=FALSE, stringsAsFactors=FALSE)

# Set column names
colnames(gff) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

# Define a function to merge overlapping annotations based on priority
merge_annotations <- function(df) {
  df <- arrange(df, source %in% c("protein2genome", "est2genome", "repeat"))  # Sort by source in descending order of priority
  merged_annotation <- paste(df$attributes, collapse=";")  # Merge attributes
  df[1, 9] <- merged_annotation  # Update the merged column in the first row
  return(df[1, ])  # Return the first row with merged information
}

# Group by genomic coordinates and apply the merge_annotations function
result <- gff %>%
  group_by(seqid, start, end) %>%
  do(merge_annotations(.))

# Write the result to a new GFF file
write.table(result, file="maker_overlap_removed.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("GenomicRanges")

library(GenomicRanges)

# Read the GFF files
file1 <- read.table("/home/bloodmark/workarea/new_final_chr/annotation/liftoff_maker_newref/maker_output_loc_corr_tsv.gff", header = FALSE, stringsAsFactors = FALSE, fill = TRUE, comment.char = "#")

# Filter out rows based on the values in column 1
correct_rows <- grepl("ragtag|chromosome|scaff|contigs|NC_|unplaced_contigs|y_chromosome", file1$V1, ignore.case = TRUE)
file1 <- file1[correct_rows, ]

# Remove duplicate rows from file1
file1 <- unique(file1)

file1 <- file1[!duplicated(file1[, c(1, 4, 5)]), ]

# Read the second GFF file
file2 <- read.table("maker_overlap_removed.gff", header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

# If file2 has more than 9 columns, merge columns 9 onwards into a single column
if (ncol(file2) > 9) {
  # Merge columns 9 onwards into a single column
  file2$V9 <- apply(file2[, 9:ncol(file2)], 1, function(x) paste(x, collapse = ";"))
  
  # Keep only the first 9 columns
  file2 <- file2[, 1:9, drop = FALSE]
}

# Concatenate the two files
merged_file <- rbind(file1, file2)

# Sort the merged file by chromosome and start position
merged_file <- merged_file[order(merged_file$V1, as.numeric(merged_file$V4)), ]

unique(merged_file$V2)
unique(merged_file$V3)

# Initialize an empty data frame to store the non-overlapping annotations
result <- data.frame(matrix(ncol = ncol(merged_file), nrow = 0))

# Create a progress bar
pb <- progress_bar$new(total = length(unique(merged_file$V1)), format = "[:bar] :percent :eta", clear = FALSE)

# Loop through each unique chromosome in the merged file
for (chromosome in unique(merged_file$V1)) {
  # Subset the merged file for the current chromosome
  subset_df <- merged_file[merged_file$V1 == chromosome, ]
  
  # Initialize a variable to store the last end position
  last_end <- -Inf
  
  # Loop through each row in the subset
  for (i in 1:nrow(subset_df)) {
    # Check for coordinate overlap with the last annotation
    if (subset_df$V4[i] <= last_end) {
      # Overlap found, prioritize annotations with "Liftoff" in column 2
      liftoff_annotations <- subset_df[subset_df$V2 == "Liftoff", ]
      if (nrow(liftoff_annotations) > 0) {
        result <- rbind(result, liftoff_annotations[1, ])
      } else {
        result <- rbind(result, subset_df[i, ])
      }
    } else {
      # No overlap found, retain all annotations
      result <- rbind(result, subset_df[i, ])
    }
    
    # Update the last end position
    last_end <- max(last_end, subset_df$V5[i])
  }
  
  # Update the progress bar
  pb$tick()
}

# Close the progress bar
pb$close()

# Write the result to a new GFF file
write.table(result, file = "maker_overlap_removed_liftoff_correct.gff", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Extract the "Name=" information from the last column
result$NewColumn <- sub('.*Name=([^;]+);.*', '\\1', result$V9)

# Write the new data frame to a new file
write.table(result, "maker_annotation_col_corr.gff", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

data <- read.csv(file="maker_annotation_col_corr.gff", header = FALSE)

