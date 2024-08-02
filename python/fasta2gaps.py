#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:23:32 2023

@author: bloodmark
"""

# README
# This script performs sequence alignment between two sequences from FASTA files, identifies unaligned regions, and writes those regions to a BED file.
# The script performs the following steps:
# 1. Reads sequences from two FASTA files.
# 2. Performs global sequence alignment.
# 3. Identifies unaligned regions in the alignment.
# 4. Writes the unaligned regions to a BED file.

# Author: Bloodmark

# Input details:
# - Input FASTA files:
#   - '/home/bloodmark/workarea/new_final_chr/chrY/y_chr_scaff.fasta'
#   - '/home/bloodmark/workarea/pubref/pubref_chromosomes/LGY.fa'
#   - These files should contain sequences in FASTA format.

# Output details:
# - Output BED file: 'LGY.bed'
#   - The file will contain columns: chromosome, start position, end position.

# Script functioning:
# 1. Reads sequences from two FASTA files.
# 2. Performs global sequence alignment.
# 3. Identifies unaligned regions in the alignment.
# 4. Writes the unaligned regions to a BED file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input files are in the specified paths and properly formatted.

from Bio import SeqIO
from Bio import pairwise2

# Function to perform sequence alignment
def align_sequences(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
    return alignments[0]

# Function to identify unaligned regions
def identify_unaligned_regions(alignment, seq_length):
    unaligned_regions = []
    start = None
    for i in range(len(alignment.seqA)):
        if alignment.seqA[i] == '-' and alignment.seqB[i] == '-':
            if start is None:
                start = i
        else:
            if start is not None:
                end = i - 1
                unaligned_regions.append((start, end))
                start = None

    # Check if the last region extends to the end of the sequence
    if start is not None:
        end = seq_length - 1
        unaligned_regions.append((start, end))

    return unaligned_regions

# Function to write BED file
def write_bed_file(regions, output_file):
    with open(output_file, 'w') as bed_file:
        for start, end in regions:
            bed_file.write(f"chr1\t{start}\t{end}\n")

# Read fasta files
fasta1 = SeqIO.read("/home/bloodmark/workarea/new_final_chr/chrY/y_chr_scaff.fasta", "fasta")
fasta2_iterator = SeqIO.parse("/home/bloodmark/workarea/pubref/pubref_chromosomes/LGY.fa", "fasta")
fasta2 = next(fasta2_iterator, None)

if fasta2 is None:
    print("Error: No records found in the second FASTA file.")
else:
    # Perform sequence alignment
    alignment = align_sequences(fasta1.seq, fasta2.seq)

    # Identify unaligned regions
    seq_length = len(fasta1.seq)
    unaligned_regions = identify_unaligned_regions(alignment, seq_length)

    # Write BED file
    write_bed_file(unaligned_regions, "LGY.bed")

