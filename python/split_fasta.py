#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Script to Split FASTA File by Gene Name

#Description:
#This script reads a multi-FASTA file and splits it into individual files based on gene names. Each output file contains sequences corresponding to a single gene. The gene name is extracted from the header line of each sequence in the input file.

#Author: Bloodmark

#Input Details:
#- The input should be a multi-FASTA file with the following format in the header line: `>sequence_id|gene_name|other_info`.
#- The script expects the input file path as a command-line argument.

#Output Details:
#- The output will be a set of files, each named after the gene name extracted from the headers.
#- Each output file will contain the sequences corresponding to that gene.

#How to Use:
#1. Save the script as `split_fasta.py`.
#2. Run the script from the command line with the input file path as an argument:

import sys

def split_fasta_by_genename(input_file):
    with open(input_file, 'r') as infile:
        outfile = None

        for line in infile:
            if line.startswith(">"):
                if outfile:
                    outfile.close()
                genename = line.strip().split('|')[1]
                filename = genename + ".txt"
                outfile = open(filename, 'w')
                outfile.write(line)
            else:
                if outfile:
                    outfile.write(line)

        if outfile:
            outfile.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_fasta_file>")
    else:
        split_fasta_by_genename(sys.argv[1])

