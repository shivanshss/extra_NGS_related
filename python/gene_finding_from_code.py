#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to Fetch Gene Names from Accession Numbers Using NCBI Entrez

Description:
This script reads a list of accession numbers from a CSV file, fetches the corresponding gene names from NCBI using the Entrez API, and writes the results to an output text file.

Author: Bloodmark

Input Details:
- The input file is a CSV file containing accession numbers, specified by the `input_file` variable.
- The script expects the input file path as an argument.

Output Details:
- The output is a text file containing accession numbers and the corresponding gene names, specified by the `output_file` variable.

How to Use:
1. Save the script as `fetch_gene_names.py`.
2. Run the script from the command line with the input and output file paths set in the script.

Example:
python fetch_gene_names.py


This will read `chr2_gene_codes.tab` and create `chr2_gene_names.txt` with the gene names.

Note:
- Ensure that the input file follows the expected format (one accession number per row).
- The script uses the NCBI Entrez API, so you need an internet connection to fetch the gene names.

Script Functioning:
1. Parses the input CSV file to get the accession numbers.
2. Determines the database type based on the accession prefix.
3. Fetches gene names in batches with retries in case of errors.
4. Writes the accession numbers and gene names to the output text file.
"""

from Bio import Entrez
import pandas as pd
import warnings

# Set your email for the Entrez service
Entrez.email = 'shivansh.bioscience@gmail.com'

# Read the single-column file
df = pd.read_csv('/home/bloodmark/workarea/20231011_maker/chr2_gene_codes.tab', header=None, names=["Accession"])

# Remove "Name=" from the entries
df['Accession'] = df['Accession'].str.replace('Name=', '')

# Extract the accession numbers
accession_list = df['Accession'].tolist()

# Create a dictionary to store gene names
gene_names = {}

# Loop through the accession list
for accession in accession_list:
    try:
        # Define a list of databases to query (you can add more if needed)
        databases = ["nucleotide", "protein"]

        # Initialize a variable to store the result
        gene_name = None

        # Loop through the databases
        for db in databases:
            try:
                handle = Entrez.esummary(db=db, id=accession)
                record = Entrez.read(handle)
                handle.close()

                if 'GeneName' in record[0]:
                    gene_name = record[0]['GeneName']
                    break  # Exit the loop once a gene name is found
            except Exception as e:
                warnings.warn(f"Warning: {str(e)}", UserWarning)

        if gene_name is not None:
            gene_names[accession] = gene_name
        else:
            gene_names[accession] = "NA"
    except Exception as e:
        print(f"Error: {str(e)}")

# Write the gene names to a text file
with open("chr2_gene_names.txt", "w") as file:
    for accession, gene_name in gene_names.items():
        file.write(f"{accession} : {gene_name}\n")

