# README
# This script processes a TSV file to create a new data frame with reference and query chromosome names, positions, and annotation types.
# The script performs the following steps:
# 1. Reads the input TSV file.
# 2. Processes each pair of rows in the input data to extract relevant information.
# 3. Determines the annotation type based on the input data.
# 4. Creates a new data frame with the extracted information.
# 5. Writes the result to a new TSV file.

# Author: Bloodmark

# Input details:
# - Input TSV file: "/home/bloodmark/workarea/synteny/LG2/5000/blocks_coords_cleaned1.txt"
#   - The file should contain the necessary columns for processing.

# Output details:
# - Output TSV file: "lg2_blocks_coords_plotsr.txt"
#   - The file will contain columns: Reference_chromosome_name, Reference_start_position, Reference_end_position, Query_chromosome_name, Query_start_position, Query_end_position, Annotation_type.

# Script functioning:
# 1. Reads the input TSV file into a data frame.
# 2. Initializes vectors to store the new columns.
# 3. Processes each pair of rows in the input data.
# 4. Determines the annotation type based on the input data.
# 5. Creates a new data frame with the new columns.
# 6. Writes the result to a new TSV file.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

# Read the input TSV file
data <- read.table(file="/home/bloodmark/workarea/synteny/LG2/5000/blocks_coords_cleaned1.txt", header = FALSE)
ref_name="lg2_new"
query_name="lg2_old"

# Initialize vectors to store the new columns
ref_chr_name <- character()
ref_start_pos <- integer()
ref_end_pos <- integer()
query_chr_name <- character()
query_start_pos <- integer()
query_end_pos <- integer()
annotation_type <- character()

# Process each pair of rows in the input data
for (i in seq(1, nrow(data), by = 2)) {
  ref_chr_name <- c(ref_chr_name, paste(ref_name, data[i, 2], sep = ""))
  ref_start_pos <- c(ref_start_pos, data[i, 4])
  ref_end_pos <- c(ref_end_pos, data[i, 5])
  
  query_chr_name <- c(query_chr_name, paste(query_name, data[i + 1, 2], sep = ""))
  query_start_pos <- c(query_start_pos, data[i + 1, 4])
  query_end_pos <- c(query_end_pos, data[i + 1, 5])
  
  # Determine annotation type
  if (data[i, 6] > data[i + 1, 6]) {
    annotation_type <- c(annotation_type, "INS")
  } else if (data[i, 6] < data[i + 1, 6]) {
    annotation_type <- c(annotation_type, "DEL")
  } else {
    annotation_type <- c(annotation_type, NA)
  }
}

# Create a data frame with the new columns
result_data <- data.frame(
  Reference_chromosome_name = ref_chr_name,
  Reference_start_position = ref_start_pos,
  Reference_end_position = ref_end_pos,
  Query_chromosome_name = query_chr_name,
  Query_start_position = query_start_pos,
  Query_end_position = query_end_pos,
  Annotation_type = annotation_type
)

# Write the result to a new TSV file
write.table(result_data, file = "lg2_blocks_coords_plotsr.txt", sep = "\t", row.names = FALSE, quote = FALSE)

