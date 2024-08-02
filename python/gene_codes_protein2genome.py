#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 19:06:25 2023

@author: bloodmark
"""

# README
# This script extracts accession numbers from a GFF file, determines the database for each accession number,
# fetches the gene names from the NCBI database, and writes the results to an output CSV file.
# The script performs the following steps:
# 1. Reads the accession numbers from the ninth column of a GFF file.
# 2. Determines the database for each accession number based on its format.
# 3. Fetches gene names from the NCBI database with retries for each accession number.
# 4. Writes the accession numbers and their corresponding gene names to an output CSV file.

# Author: Bloodmark

# Input details:
# - Input GFF file: '/home/bloodmark/protein_match.tsv'
#   - The file should contain accession numbers in the ninth column.

# Output details:
# - Output CSV file: '/home/bloodmark/protein_match_names.txt'
#   - The file will contain columns: Accession Number, Gene Name.

# Script functioning:
# 1. Reads the accession numbers from the input GFF file.
# 2. Determines the database for each accession number.
# 3. Fetches gene names from the NCBI database with retries for each accession number.
# 4. Writes the accession numbers and their corresponding gene names to an output CSV file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.
# Note: Ensure you have internet access to fetch gene names from the NCBI database.

import csv
from Bio import Entrez
from Bio import SeqIO
import time
import re

# Set your email for the Entrez service
Entrez.email = "shivansh.bioscience@gmail.com"

# Function to determine the database for a given accession number
def determine_database(accession):
    if accession.startswith("XP_") or accession.startswith("NP_"):
        return "protein"
    elif re.match(r'^gi\|\d+\|\w+\|\w+\.\d+\|$', accession):
        return "nucleotide"
    else:
        return None

# Function to extract gene name from the description
def extract_gene_name(description):
    parts = description.split('[')
    if len(parts) > 1:
        gene_name = parts[0]
    else:
        gene_name = description
    return gene_name.strip()

# Function to fetch gene names with retries for a batch of accession numbers
def fetch_gene_names_for_batch(accession_batch, max_retries=3):
    gene_names = []
    for accession in set(accession_batch):
        database = determine_database(accession)
        if database is None:
            gene_names.append((accession, 'Error: Unsupported accession format'))
            continue
        retries = 0
        while retries < max_retries:
            try:
                handle = Entrez.efetch(db=database, id=accession, rettype="gb", retmode="text")
                record = SeqIO.read(handle, "genbank")
                description = record.description
                gene_name = extract_gene_name(description)
                gene_names.append((accession, gene_name))
                break
            except Exception as e:
                retries += 1
                if retries < max_retries:
                    time.sleep(1)
                else:
                    gene_names.append((accession, f'Error: {str(e)}'))
                    print(f"Error for accession {accession}: {str(e)}")
    return gene_names

# Read the input GFF file
input_file = '/home/bloodmark/protein_match.tsv'
output_file = '/home/bloodmark/protein_match_names.txt'

# Extract accession numbers from the ninth column of the GFF file
accession_numbers = []
with open(input_file, 'r') as gff_file:
    for line_number, line in enumerate(gff_file, start=1):
        if not line.startswith('#'):
            columns = line.strip().split('\t')
            if len(columns) >= 9:
                match = re.search(r'(\bgi\|\d+\|\w+\|\w+\.\d+\|)', columns[8])
                if match:
                    accession_numbers.append(match.group(1))
        if line_number == 10:
            break

# Fetch gene names in batches
batch_size = 10
gene_names = []
for i in range(0, len(accession_numbers), batch_size):
    accession_batch = accession_numbers[i:i + batch_size]
    batch_gene_names = fetch_gene_names_for_batch(accession_batch)
    gene_names.extend(batch_gene_names)

# Write results to an output CSV file
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Accession Number', 'Gene Name'])
    for accession, gene_name in gene_names:
        writer.writerow([accession, gene_name])

print(f"Gene names for {len(accession_numbers)} accession numbers have been fetched and saved to {output_file}.")

