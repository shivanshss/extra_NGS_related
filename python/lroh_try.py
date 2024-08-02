#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to Read and Plot Long Runs of Homozygosity (LROH)

Description:
This script reads LROH data from a specified file, extracts the start and end positions of LROH segments, and plots them using matplotlib. The LROH segments are displayed as shaded regions along the x-axis representing genomic positions.

Author: Bloodmark

Input Details:
- The input file is a .LROH file containing LROH segment data, specified by the `lroh_file` variable.

Output Details:
- The script generates a plot showing the LROH segments along the genomic positions.

How to Use:
1. Save the script as `plot_lroh.py`.
2. Update the `lroh_file` variable with the path to your .LROH file.
3. Run the script from the command line.

Example:
python plot_lroh.py


This will read the specified .LROH file and generate a plot of the LROH segments.

Note:
- Ensure that the .LROH file follows the expected format with start and end positions in the second and third columns respectively.
- The script uses matplotlib for plotting, so ensure that the package is installed.

Script Functioning:
1. Parses the .LROH file to extract start and end positions of LROH segments.
2. Plots the LROH segments using matplotlib.

Script:
"""

import matplotlib.pyplot as plt

def read_lroh_file(lroh_file):
    """
    Reads the LROH file and extracts the start and end positions of LROH segments.
    
    Parameters:
    lroh_file (str): Path to the .LROH file.
    
    Returns:
    list: List of tuples containing start and end positions of LROH segments.
    """
    lroh_data = []
    with open(lroh_file, 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                continue
            # Split the line into columns
            columns = line.strip().split()
            # Check if the line has the required number of columns
            if len(columns) < 3:
                continue
            # Extract start and end positions of the LROH segment
            try:
                start = int(columns[1])
                end = int(columns[2])
                lroh_data.append((start, end))
            except ValueError:
                # Skip lines with invalid data
                continue
    return lroh_data

def plot_lroh(lroh_data):
    """
    Plots the LROH segments.
    
    Parameters:
    lroh_data (list): List of tuples containing start and end positions of LROH segments.
    """
    fig, ax = plt.subplots(figsize=(10, 5))
    for start, end in lroh_data:
        ax.axvspan(start, end, alpha=0.3, color='blue')
    ax.set_xlabel("Position")
    ax.set_ylabel("LROH")
    ax.set_title("Long Runs of Homozygosity")
    plt.show()

if __name__ == "__main__":
    lroh_file = "/media/bloodmark/HDD6_SS_extra/w_newref/5.vcfstats/vcfstats_filter4/out.LROH"  # Replace with the path to your .LROH file
    lroh_data = read_lroh_file(lroh_file)
    plot_lroh(lroh_data)

