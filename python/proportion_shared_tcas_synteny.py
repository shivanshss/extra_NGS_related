#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:58:50 2024

@author: bloodmark
"""

# README
# This script reads a delta file and extracts information about chromosome alignments.
# It calculates the proportion of the genome shared between each pair of new and old chromosomes and
# prints a matrix representing these proportions.

# Author: Bloodmark

# Input details:
# - Input delta file: "/media/bloodmark/HDD6_SS_extra/genome_assembly_figure_data/nucmer/tcas6_final_all.delta"
#   - The file should contain alignment information in delta format.

# Output details:
# - Prints a matrix with the proportion of genome shared between new and old chromosomes.
# - The matrix includes the chromosome names in rows and columns.

# Script functioning:
# 1. Reads the delta file and extracts chromosome alignment information.
# 2. Calculates the proportion of the genome shared between each pair of new and old chromosomes.
# 3. Creates a matrix representing these proportions.
# 4. Prints the matrix with row and column names.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

# Initialize variables and data structures
chromosome_names = set()  # To store unique chromosome names
chromosome_data = {}  # To store information for each chromosome pair

# Read the delta file and extract relevant information
with open("/media/bloodmark/HDD6_SS_extra/genome_assembly_figure_data/nucmer/tcas6_final_all.delta", "r") as delta_file:
    lines = delta_file.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            # Extract chromosome names and positions
            header = lines[i].split()
            new_chromosome = header[0][1:]  # Remove '>' from the new chromosome name
            old_chromosome = header[1]
            i += 1
            if i < len(lines):
                new_positions = list(map(int, lines[i].split()))
                i += 1
                if i < len(lines):
                    old_positions = list(map(int, lines[i].split()))
                    
                    # Check if both positions lists contain enough elements
                    if len(new_positions) >= 2 and len(old_positions) >= 2:
                        # Calculate the proportion of genome shared
                        shared_proportion = len(set(range(new_positions[0], new_positions[1] + 1)).intersection(range(old_positions[0], old_positions[1] + 1))) / float(new_positions[1] - new_positions[0] + 1)
                        
                        # Store the information in the dictionary
                        chromosome_data[(new_chromosome, old_chromosome)] = shared_proportion
                        
                        # Add chromosome names to the set
                        chromosome_names.add(new_chromosome)
                        chromosome_names.add(old_chromosome)
        
        i += 1

# Create the matrix for the proportion of genome shared
matrix = [[0.0] * len(chromosome_names) for _ in range(len(chromosome_names))]
row_index = {name: index for index, name in enumerate(chromosome_names)}

for (new_chromosome, old_chromosome), shared_proportion in chromosome_data.items():
    matrix[row_index[new_chromosome]][row_index[old_chromosome]] = shared_proportion

# Print the matrix with row and column names
print("new genome chromosome names in rows:")
print(list(chromosome_names))
print("old genome chromosome names in columns:")
print(list(chromosome_names))
print("Proportion of genome shared matrix:")
for row in matrix:
    print(row)

