#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:17:59 2023

@author: bloodmark
"""
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

def find_best_match(query_sequence, target_sequences):
    # Perform BLAST search
    result_handle = NCBIWWW.qblast("blastn", "nt", query_sequence.seq)

    # Parse the BLAST result
    blast_records = NCBIXML.read(result_handle)
    result_handle.close()

    # Check if there are any alignments
    if blast_records.alignments:
        # Extract the best match
        best_alignment = blast_records.alignments[0]

        # Extract the name of the best match
        best_match_name = best_alignment.title.split()[1]
    else:
        # No alignments found
        best_match_name = "No Match"

    return best_match_name

# Load protein sequences from file1
file1_sequences = list(SeqIO.parse("/media/bloodmark/HDD6_SS_extra/est_gene_name/insect_proteins.fasta", "fasta"))

# Load nucleotide sequences from file2
file2_sequences = list(SeqIO.parse("/media/bloodmark/HDD6_SS_extra/est_gene_name/est_seq.fasta", "fasta"))

# Open a file for writing the output
with open("output.txt", "w") as output_file:
    # Write the header
    output_file.write("Name_in_file2\tBest_match_in_file1\n")

    # Iterate over nucleotide sequences in file2
    for file2_sequence in file2_sequences:
        # Find the best match in file1
        best_match_name = find_best_match(file2_sequence, file1_sequences)

        # Write the results to the file
        output_file.write(f"{file2_sequence.id}\t{best_match_name}\n")

print("Output written to output.txt")
