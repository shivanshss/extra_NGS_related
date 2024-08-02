# README
# This script processes multiple Fst windowed files, combines them into a single data frame,
# and identifies the top 100 rows based on the WEIGHTED_FST value.
# The script assumes that each input file is in the working directory and has a header row.

# Author: Bloodmark

# Input details:
# - Input files: Multiple windowed Fst files located in the working directory
#   - Each file should be tab-delimited and contain a header row with at least one column named 'WEIGHTED_FST'

# Output details:
# - A combined data frame of all input files
# - A data frame containing the top 100 rows sorted by WEIGHTED_FST in descending order

# Script functioning:
# 1. Set the working directory to the location of your data files.
# 2. Define a list of file names to be processed.
# 3. Initialize an empty list to store data frames.
# 4. Read each file into a data frame and store it in the list.
# 5. Combine all data frames into a single data frame.
# 6. Identify the top 100 rows based on the WEIGHTED_FST value.
# 7. Print or further process the top 100 data frame as needed.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure that the input files are in the specified working directory and properly formatted.

# Code starts below:

# Set the working directory where your data files are located
setwd("/media/bloodmark/HDD6_SS_extra/ch2_vcf/ch2_stats/fst/10kb_window")

# List of file names
file_names <- c(
  "pop01_vs_pop02.windowed.weir.fst", "pop01_vs_pop11.windowed.weir.fst",
  "pop01_vs_pop12.windowed.weir.fst", "pop01_vs_pop13.windowed.weir.fst",
  "pop01_vs_pop18.windowed.weir.fst", "pop01_vs_pop20.windowed.weir.fst",
  "pop01_vs_pop24.windowed.weir.fst", "pop02_vs_pop11.windowed.weir.fst",
  "pop02_vs_pop12.windowed.weir.fst", "pop02_vs_pop13.windowed.weir.fst",
  "pop02_vs_pop18.windowed.weir.fst", "pop02_vs_pop20.windowed.weir.fst",
  "pop02_vs_pop24.windowed.weir.fst", "pop11_vs_pop12.windowed.weir.fst",
  "pop11_vs_pop13.windowed.weir.fst", "pop11_vs_pop18.windowed.weir.fst",
  "pop11_vs_pop20.windowed.weir.fst", "pop11_vs_pop24.windowed.weir.fst",
  "pop12_vs_pop13.windowed.weir.fst", "pop12_vs_pop18.windowed.weir.fst",
  "pop12_vs_pop20.windowed.weir.fst", "pop12_vs_pop24.windowed.weir.fst",
  "pop13_vs_pop18.windowed.weir.fst", "pop13_vs_pop20.windowed.weir.fst",
  "pop13_vs_pop24.windowed.weir.fst", "pop18_vs_pop20.windowed.weir.fst",
  "pop18_vs_pop24.windowed.weir.fst", "pop20_vs_pop24.windowed.weir.fst"
)

# Create an empty list to store all data frames
dfs <- list()

# Read each file into a data frame and store in the list
for (file in file_names) {
  df <- tryCatch({
    read.table(file, header = TRUE)  # Adjust parameters based on your actual data file
  }, error = function(e) {
    message(sprintf("Error reading file %s: %s", file, e))
    NULL
  })
  
  if (!is.null(df)) {
    dfs[[file]] <- df
  }
}

# Combine all data frames into one
combined_df <- do.call(rbind, dfs)

# Step 2: Identify Top 100 WEIGHTED_FST Values
# Sort combined_df by WEIGHTED_FST in descending order and select top 100
top_100 <- combined_df[order(-combined_df$WEIGHTED_FST), ][1:100, ]

# Print or further process top_100 dataframe as needed
print(top_100)

# Optionally, save the top 100 dataframe to a file
write.table(top_100, "top_100_weighted_fst.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

