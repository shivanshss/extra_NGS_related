#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:48:38 2023

@author: bloodmark
"""

# README
# This script reads a delta file and extracts information about insertions, deletions, and gaps.
# The script performs the following steps:
# 1. Reads the delta file line by line.
# 2. Extracts fields from lines starting with '>'.
# 3. Identifies insertions, deletions, and gaps based on specific markers in the fields.
# 4. Prints the information about insertions, deletions, and gaps.

# Author: Bloodmark

# Input details:
# - Input delta file: '/home/bloodmark/workarea/nucmer/tcas6_final_all.delta'
#   - The file should contain alignment information in delta format.

# Output details:
# - The script prints information about insertions, deletions, and gaps.

# Script functioning:
# 1. Reads the delta file line by line.
# 2. Extracts fields from lines starting with '>'.
# 3. Identifies insertions, deletions, and gaps based on specific markers in the fields.
# 4. Prints the information about insertions, deletions, and gaps.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

delta_file_path = '/home/bloodmark/workarea/nucmer/tcas6_final_all.delta'

with open(delta_file_path, 'r') as delta_file:
    for line in delta_file:
        if line.startswith('>'):
            fields = line.split()
            if fields[0] == '>1':
                if fields[1] == 'I':
                    print(f"Insertion at position {fields[2]} in sequence 2, length {fields[3]}")
                elif fields[1] == 'D':
                    print(f"Deletion at position {fields[2]} in sequence 1, length {fields[3]}")
                elif fields[1] == 'G':
                    print(f"Gap at position {fields[2]} in both sequences, length {fields[3]}")

