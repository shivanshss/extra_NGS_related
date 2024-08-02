# README
# This script processes and merges two data files: a TSV file and a CSV file. The steps include:
# 1. Reading the data files into data frames.
# 2. Cleaning and transforming the data.
# 3. Merging the data frames based on a common column.
# 4. Handling missing values and removing duplicates.
# 5. Writing the merged data to a new file.

# Author: Bloodmark

# Input details:
# - Input TSV file: "/home/bloodmark/workarea/20231011_maker/chr4_split.tab"
# - Input CSV file: "/home/bloodmark/workarea/20231011_maker/chr4_non_gi_gene_names.txt"

# Output details:
# - Output file: "chr4.gff" containing the merged and cleaned data.

# Script functioning:
# 1. Installs and loads the required packages (dplyr, readr, data.table).
# 2. Reads the input TSV and CSV files into data frames.
# 3. Cleans and transforms the data by replacing blank spaces with "NA" and modifying specific columns.
# 4. Merges the data frames using a left join.
# 5. Fills missing values with "NA" and removes duplicate rows.
# 6. Writes the cleaned and merged data to a new file in the same format.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install required packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# Load the required packages
library(dplyr)
library(readr)
library(data.table)

# Read the TSV file (file1) into a data frame
file1 <- fread("/home/bloodmark/workarea/20231011_maker/chr4_split.tab", sep = "\t", header = TRUE)

# Replace blank spaces with "NA" in file1
file1[file1 == ""] <- "NA"

# Read the CSV file (file2) into a data frame
file2 <- read.csv("/home/bloodmark/workarea/20231011_maker/chr4_non_gi_gene_names.txt", header = TRUE)

# Remove "name=" from column 6 in file1
file1$split_data_df.Col2 <- gsub("Name=", "", file1$split_data_df.Col2)

# Rename the column in file2 to match the column in file1
colnames(file2)[colnames(file2) == "Accession.Number"] <- "split_data_df.Col2"

# Perform a left join to merge the data frames
merged_data <- left_join(file1, file2, by = "split_data_df.Col2")

# Fill missing values with 'NA'
merged_data[is.na(merged_data$Gene.Name), "Gene.Name"] <- "NA"

# Remove duplicate rows
merged_data <- distinct(merged_data)

# Write the result to a new file with the same format
write.table(merged_data, "chr4.gff", sep = "\t", row.names = FALSE, quote = FALSE)

