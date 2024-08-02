#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 21:59:28 2023

@author: bloodmark
"""

# README
# This script reads a gap fill detail file, splits its content based on headers, and writes each section to a separate file.
# The script performs the following steps:
# 1. Reads the input file line by line.
# 2. Creates a new output file each time a header line (starting with ">") is encountered.
# 3. Writes the relevant lines to the respective output file.
# 4. Closes all output files at the end.

# Author: Bloodmark

# Input details:
# - Input file: "/media/bloodmark/HDD6_SS_extra/tgs_gap_closer/7_ont/7_ont.gap_fill_detail"
#   - The file should contain sections starting with headers (lines starting with ">").

# Output details:
# - Output files: Multiple files with names based on headers (without ">") followed by ".txt".
#   - Each file contains the section of the input file corresponding to its header.

# Script functioning:
# 1. Reads the input file line by line.
# 2. Creates a new output file each time a header line (starting with ">") is encountered.
# 3. Writes the relevant lines to the respective output file.
# 4. Closes all output files at the end.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

input_file = "/media/bloodmark/HDD6_SS_extra/tgs_gap_closer/7_ont/7_ont.gap_fill_detail"

with open(input_file, 'r') as file:
    lines = file.readlines()

output_files = {}
current_output_file = None

for line in lines:
    if line.startswith(">"):
        # Close the current output file if it exists
        if current_output_file:
            current_output_file.close()

        # Create a new output file with the filename
        filename = line.strip().replace(">", "") + ".txt"
        current_output_file = open(filename, 'w')
        output_files[filename] = current_output_file

    if current_output_file:
        current_output_file.write(line)

# Close all output files
for output_file in output_files.values():
    output_file.close()

print("Files have been split and saved.")

