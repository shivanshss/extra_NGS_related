#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 20:25:05 2024

@author: bloodmark
"""

# README
# This script generates random DNA sequences, aligns them using Clustal Omega, and plots the alignment.
# The script performs the following steps:
# 1. Generates random DNA sequences.
# 2. Writes the sequences to a fasta file.
# 3. Performs alignment using Clustal Omega.
# 4. Reads the alignment and converts it to a matrix.
# 5. Plots the alignment.

# Author: Bloodmark

# Input details:
# - Number of sequences: 25
# - Sequence length: 1000

# Output details:
# - Output fasta file: "sequences.fasta"
# - Aligned fasta file: "aligned.fasta"
# - Multiple sequence alignment plot.

# Script functioning:
# 1. Generates random DNA sequences of specified length.
# 2. Creates SeqRecord objects for each sequence.
# 3. Writes the sequences to a fasta file.
# 4. Calls Clustal Omega to perform the alignment.
# 5. Reads the aligned sequences.
# 6. Converts the alignment to a matrix and plots it.

# Example of running the script:
# python your_script_name.py

# Note: Ensure Clustal Omega is installed and accessible in your system's PATH.

import random
import os
import matplotlib.pyplot as plt
import numpy as np
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

# Generate random sequences
sequences = [generate_random_sequence(1000) for _ in range(25)]

# Create SeqRecord objects
seq_records = [SeqRecord(Seq(seq), id=f"Seq{i+1}") for i, seq in enumerate(sequences)]

# Write sequences to a fasta file
SeqIO.write(seq_records, "sequences.fasta", "fasta")

# Perform alignment using Clustal Omega with force overwrite
clustalomega_cline = ClustalOmegaCommandline(infile="sequences.fasta", outfile="aligned.fasta", verbose=True, auto=True, force=True)
os.system(str(clustalomega_cline))

# Read the alignment
alignment = AlignIO.read("aligned.fasta", "fasta")

# Convert alignment to a matrix
alignment_array = np.array([list(rec) for rec in alignment], np.dtype('S1'))

# Plot the alignment
plt.figure(figsize=(20, 10))
plt.imshow(alignment_array == b'-', cmap='Greys', interpolation='none', aspect='auto')
plt.xlabel("Position")
plt.ylabel("Sequence")
plt.title("Multiple Sequence Alignment")
plt.show()

