#!/bin/bash

# README
# This script processes depth files in a specified directory by adding column headers and then merging them using a Python script.
# The script iterates over each file in the directory, adds column headers, and appends the original content.
# Finally, it uses a Python script to merge the files based on the first two columns if the number of rows is consistent across all files.

# Author: Bloodmark

# Input details:
# - Directory containing the depth files: 
# - Each file in the directory should be a regular text file without column headers.

# Output details:
# - Processed files with headers: Original files with "_with_headers.txt" appended to their names.
# - Merged output file: merged_output.txt in the specified directory.

# Script functioning:
# 1. Define the base directory for the operation.
# 2. Iterate over each file in the base directory.
# 3. Add column headers to each file and save as a new file.
# 4. Use a Python script to merge the processed files based on the first two columns.
# 5. If the number of rows is not consistent across all files, print an error message.

# Example of running the script:
# bash your_script_name.sh

# Note: Ensure Python 3 and pandas are installed and accessible.

# Code starts below:

# Define the directory where you want to perform the operation
base_dir="/media/bloodmark/HDD6_SS_extra/chapter2/w_newref_old/3.depth/LG10"

# Iterate over each file in the base directory
for file in "$base_dir"/*; do
    # Skip if not a regular file
    [[ -f "$file" ]] || continue

    filename=$(basename "$file")
    echo "Processing file: $filename"

    # Add column names to the file
    echo "Chr, Pos, ${filename%.*}" > "${file}_with_headers.txt"
    cat "$file" >> "${file}_with_headers.txt"
done

# Merge files using Python script
python3 <<EOF
import pandas as pd
import glob

# Get all files with headers
files_with_headers = glob.glob("$base_dir/*_with_headers.txt")

# Read all files into a list of dataframes
dfs = [pd.read_csv(file, sep=', ') for file in files_with_headers]

# Check if the number of rows is the same for all files
if all(len(df) == len(dfs[0]) for df in dfs):
    # If the number of rows is the same, merge based on column 1 and 2
    merged_df = pd.concat([df.iloc[:, :2] for df in dfs], axis=1)
    # Add depth columns
    for df in dfs:
        depth_column = df.columns[-1]
        merged_df[depth_column] = df[depth_column]
    # Write merged data to file
    merged_df.to_csv("$base_dir/merged_output.txt", sep=',', index=False)
else:
    print("Number of rows is not the same for all files. Cannot merge.")
EOF

echo "Processing completed for files in folder: $base_dir"

