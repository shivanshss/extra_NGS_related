#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:22:35 2023

@author: bloodmark
"""
def parse_delta_file(delta_file, output_bed):
    with open(delta_file, "r") as delta:
        lines = delta.readlines()

    with open(output_bed, "w") as bed_file:
        in_diff_section = False
        ref_sequence_id = None
        skipped_lines = 0  # Initialize the counter for skipped lines
        for line in lines:
            if line.startswith("NUCMER"):
                in_diff_section = True
                continue
            if in_diff_section:
                fields = line.strip().split()
                if not fields or fields[0].startswith('>'):
                    # Skip empty lines and lines starting with '>'
                    skipped_lines += 1  # Increment the counter for skipped lines
                    continue

                # Handle lines with only -1 values
                if all(value == "-1" for value in fields):
                    skipped_lines += 1  # Increment the counter for skipped lines
                    continue

                # Ensure there are enough fields in the line
                if len(fields) >= 2:
                    mutation_type = "SNP"  # Assume SNP by default
                    start, end = int(fields[0]), int(fields[1])
                    
                    if len(fields) >= 3:
                        length = int(fields[2])
                        if length > 0:
                            mutation_type = "INS"
                            end = start + length
                        elif length < 0:
                            mutation_type = "DEL"
                            start = end

                    bed_file.write(f"{ref_sequence_id}\t{start}\t{end}\t{mutation_type}\n")
                else:
                    print(f"Skipping line: {line.strip()} (Insufficient fields)")

        print(f"Skipped {skipped_lines} lines with insufficient fields.")

if __name__ == "__main__":
    delta_file = "/home/bloodmark/workarea/nucmer/tcas6_final_all.delta"
    output_bed = "tcas6_allchr.bed"
    parse_delta_file(delta_file, output_bed)
