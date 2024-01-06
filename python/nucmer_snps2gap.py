#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:19:39 2023

@author: bloodmark
"""

with open("NC_007417.snps.txt", "r") as snps_file, open("gaps.bed", "w") as bed_file:
    for line in snps_file:
        if line.startswith("Contig"):
            fields = line.split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[3])
            if fields[2] == ".":
                bed_file.write(f"{chrom}\t{start}\t{end}\n")