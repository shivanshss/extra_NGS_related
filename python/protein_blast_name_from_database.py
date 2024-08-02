#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:24:10 2024

@author: bloodmark
"""

# README
# This script performs BLAST searches for protein sequences in a multifasta file and writes the top 3 hits to a TSV file.
# The script performs the following steps:
# 1. Reads the protein sequences from a multifasta file.
# 2. Performs BLAST searches for each protein sequence against the NCBI 'nr' database.
# 3. Extracts information about the top 3 hits for each protein sequence.
# 4. Writes the extracted information to a TSV file.

# Author: Bloodmark

# Input details:
# - Input multifasta file: "/media/bloodmark/HDD6_SS_extra/final_maker/round3/tcas6_final_all.maker.output/tcas6_final_all.all.maker.non_overlapping_ab_initio.proteins.fasta"
#   - The file should contain protein sequences in fasta format.

# Output details:
# - Output TSV file: "tcas6_final_all.all.maker.non_overlapping_ab_initio_blast_results.tsv"
#   - The file will contain columns: Protein Name, Database Protein Name, Size, Percent Match, Percent Identity, Organism Name.

# Script functioning:
# 1. Reads the protein sequences from the input multifasta file.
# 2. Performs BLAST searches for each protein sequence against the NCBI 'nr' database.
# 3. Extracts information about the top 3 hits for each protein sequence.
# 4. Writes the extracted information to a TSV file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.
# Note: Ensure you have internet access to perform BLAST searches.

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Define the path to your multifasta file
fasta_file = "/media/bloodmark/HDD6_SS_extra/final_maker/round3/tcas6_final_all.maker.output/tcas6_final_all.all.maker.non_overlapping_ab_initio.proteins.fasta"
output_file = "tcas6_final_all.all.maker.non_overlapping_ab_initio_blast_results.tsv"

# Read the multifasta file
records = list(SeqIO.parse(fasta_file, "fasta"))

# Open the output TSV file for writing
with open(output_file, "w") as tsv_file:
    # Write the header line to the TSV file
    tsv_file.write("Protein Name\tDatabase Protein Name\tSize\tPercent Match\tPercent Identity\tOrganism Name\n")
    
    # Perform BLAST for each protein in the multifasta
    for record in records:
        result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)
        blast_records = NCBIXML.parse(result_handle)
        
        # Get the top 3 hits
        top_hits = []
        for blast_record in blast_records:
            for alignment in blast_record.alignments[:3]:
                title = alignment.title
                length = alignment.length
                percent_match = (alignment.length / len(record.seq)) * 100
                percent_identity = alignment.hsps[0].identities / alignment.hsps[0].align_length * 100
                
                # Extract the organism name from the title
                organism_name = title.split("[")[1].split("]")[0]
                
                top_hits.append({
                    "Protein Name": record.id,
                    "Database Protein Name": title,
                    "Size": length,
                    "Percent Match": percent_match,
                    "Percent Identity": percent_identity,
                    "Organism Name": organism_name
                })
        
        # Write the top hits to the TSV file
        for hit in top_hits:
            tsv_file.write(f"{hit['Protein Name']}\t{hit['Database Protein Name']}\t{hit['Size']}\t{hit['Percent Match']:.2f}\t{hit['Percent Identity']:.2f}\t{hit['Organism Name']}\n")

print(f"Results have been written to {output_file}")

