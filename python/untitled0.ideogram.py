#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:40:55 2023

@author: bloodmark
"""

from Bio import SeqIO
import subprocess

def find_mutations(file1, file2, output_bed):
    # Step 1: Create temporary files for bedtools
    temp_file1 = "temp_file1.bed"
    temp_file2 = "temp_file2.bed"

    # Step 2: Convert FASTA files to BED format using bedtools
    convert_to_bed(file1, temp_file1)
    convert_to_bed(file2, temp_file2)

    # Step 3: Use bedtools to find differences
    subprocess.run(["bedtools", "intersect", "-a", temp_file1, "-b", temp_file2, "-v", ">", "differences.bed"], shell=True)

    # Step 4: Parse differences and write to output BED file
    parse_differences("differences.bed", output_bed)

    # Step 5: Remove temporary files
    subprocess.run(["rm", temp_file1, temp_file2, "differences.bed"], shell=True)

def convert_to_bed(input_fasta, output_bed):
    with open(input_fasta, "r") as fasta_file, open(output_bed, "w") as bed_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            bed_file.write(f"{record.id}\t0\t{len(record)}\n")

def parse_differences(input_bed, output_bed):
    with open(input_bed, "r") as input_file, open(output_bed, "w") as output_file:
        for line in input_file:
            fields = line.strip().split("\t")
            start = int(fields[1])
            end = int(fields[2])
            mutation_type = determine_mutation_type(start, end)
            output_file.write(f"{fields[0]}\t{start}\t{end}\t{mutation_type}\n")

def determine_mutation_type(start, end):
    # Implement your logic to determine mutation type based on start and end positions
    # This could involve comparing with the original sequence in file1
    # Return 'ins', 'del', 'inv', or 'dup' accordingly
    return "mutation_type"

if __name__ == "__main__":
    file1 = "/home/bloodmark/workarea/new_final_chr/NC_007416.3_RagTag.fa"
    file2 = "/home/bloodmark/workarea/pubref/pubref_chromosomes/LGX.fa"
    output_bed = "LG2.bed"
    find_mutations(file1, file2, output_bed)
