#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:24:15 2023

@author: bloodmark
"""

# README
# This script reads a Nucmer coordinates file, identifies various types of structural variations (deletion, insertion, translocation, inversion), and writes the identified regions to a BED file.
# The script performs the following steps:
# 1. Reads the Nucmer coordinates file line by line.
# 2. Skips header lines and empty lines.
# 3. Extracts relevant fields from each line.
# 4. Identifies structural variations based on the extracted fields.
# 5. Writes the identified structural variations to a BED file.

# Author: Bloodmark

# Input details:
# - Input Nucmer coordinates file: "/home/bloodmark/workarea/nucmer/NC_007417.3.coords"
#   - The file should contain alignment coordinates in a tab-separated format.

# Output details:
# - Output BED file: "NC_007417.bed"
#   - The file will contain columns: chromosome, start position, end position, types of structural variations.

# Script functioning:
# 1. Reads the Nucmer coordinates file line by line.
# 2. Skips header lines and empty lines.
# 3. Extracts relevant fields from each line.
# 4. Identifies structural variations based on the extracted fields.
# 5. Writes the identified structural variations to a BED file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

with open("/home/bloodmark/workarea/nucmer/NC_007417.3.coords", "r") as nucmer_file, open("NC_007417.bed", "w") as bed_file:
    header_found = False

    for line in nucmer_file:
        if line.startswith("[S1]") or not line.strip():
            header_found = True if line.startswith("[S1]") else False
            continue  # Skip header and empty lines

        if not header_found:
            continue  # Skip lines until the header is found

        fields = line.strip().split("\t")

        # Check if the line has the expected number of fields
        if len(fields) < 9:
            print(f"Skipping line: {line}")
            continue

        chrom = "NC_007417.3"  # Update with your actual chromosome name
        try:
            chromStart = int(fields[0])
            chromEnd = int(fields[1])
            conditions_fulfilled = []

            # Check conditions for deletion, insertion, translocation, and inversion
            if chromEnd < chromStart:
                conditions_fulfilled.append("DEL")  # Deletion
            if chromEnd > chromStart:
                conditions_fulfilled.append("INS")  # Insertion
            if fields[4] != fields[5] or fields[6] != fields[7]:
                conditions_fulfilled.append("TRA")  # Translocation
            if (chromStart > chromEnd and int(fields[2]) < int(fields[6])) or (chromStart < chromEnd and int(fields[2]) > int(fields[6])):
                conditions_fulfilled.append("INV")  # Inversion

            if conditions_fulfilled:
                bed_file.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{','.join(conditions_fulfilled)}\n")

        except ValueError:
            print(f"Skipping line: {line}")
            continue

