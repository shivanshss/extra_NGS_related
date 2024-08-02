#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:19:39 2023

@author: bloodmark
"""

# README
# This script reads a SNPs file, identifies gaps, and writes the gap regions to a BED file.
# The script performs the following steps:
# 1. Reads the SNPs file line by line.
# 2. Extracts relevant fields from lines starting with "Contig".
# 3. Identifies gap regions where the third field is a dot (".").
# 4. Writes the identified gap regions to a BED file.

# Author: Bloodmark

# Input details:
# - Input SNPs file: "NC_007417.snps.txt"
#   - The file should contain SNP information in a tab-separated format.

# Output details:
# - Output BED file: "gaps.bed"
#   - The file will contain columns: chromosome, start position, end position.

# Script functioning:
# 1. Reads the SNPs file line by line.
# 2. Extracts relevant fields from lines starting with "Contig".
# 3. Identifies gap regions where the third field is a dot (".").
# 4. Writes the identified gap regions to a BED file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

with open("NC_007417.snps.txt", "r") as snps_file, open("gaps.bed", "w") as bed_file:
    for line in snps_file:
        if line.startswith("Contig"):
            fields = line.split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[3])
            if fields[2] == ".":
                bed_file.write(f"{chrom}\t{start}\t{end}\n")

