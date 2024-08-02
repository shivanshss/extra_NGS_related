# README
# This script processes a GFF file by removing duplicate rows, computing new columns, extracting unique values,
# creating separate files for each unique value in a specific column, and generating a summary table of counts.
# The script assumes that the input file is properly formatted and located in the specified path.

# Author: Bloodmark

# Input details:
# - Input file: "maker_annotation_col_corr.gff"
#   - The file should be a tab-separated values (TSV) file without a header row.

# Output details:
# - Separate files for each unique value in column 3.
# - A summary table printed to the console with counts of combinations of unique values in columns 2 and 3.

# Script functioning:
# 1. Install and load the dplyr package if not already installed.
# 2. Read the input TSV file into a data frame.
# 3. Remove duplicate rows.
# 4. Compute a new column V6 as the difference between columns V5 and V4.
# 5. Extract unique values in columns V2 and V3.
# 6. Create separate files for each unique value in column V3.
# 7. Generate a summary table with counts of combinations of unique values in columns V2 and V3.
# 8. Print the summary table.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure that the input file is in the specified path and properly formatted.

# Code starts below:

# Install and load the dplyr package if not already installed
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

library(dplyr)

# Read the TSV file
data <- tryCatch({
  read.csv("maker_annotation_col_corr.gff", header = FALSE, sep = "\t")
}, error = function(e) {
  message("Error reading the TSV file: ", e)
  NULL
})

# Proceed only if data is read successfully
if (!is.null(data)) {
  # Remove duplicate rows
  data <- distinct(data)
  
  # Compute a new column V6 as the difference between columns V5 and V4
  data <- mutate(data, V6 = V5 - V4)
  
  # Get unique values in columns V2 and V3
  unique_values_col2 <- unique(data$V2)
  unique_values_col3 <- unique(data$V3)
  
  # Initialize an empty data frame for the summary
  summary_table <- data.frame(Type = character(), Value = character(), Count = numeric(), stringsAsFactors = FALSE)
  
  # Create separate files for each unique value in column V3
  for (value in unique_values_col3) {
    subset_df <- data[data$V3 == value, ]
    tryCatch({
      write.table(subset_df, file = paste0(value, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }, error = function(e) {
      message("Error writing file for value ", value, ": ", e)
    })
  }
  
  # Loop through unique values in columns V2 and V3
  for (col2_value in unique_values_col2) {
    for (col3_value in unique_values_col3) {
      count <- sum(data$V2 == col2_value & data$V3 == col3_value)
      summary_table <- rbind(summary_table, data.frame(Type = col2_value, Value = col3_value, Count = count))
    }
  }
  
  # Print the summary table
  print(summary_table)
} else {
  message("Data reading failed, processing skipped.")
}

