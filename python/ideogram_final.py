#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 10:22:14 2023

@author: bloodmark
"""

from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def create_ideogram(file_path):
    ideogram_data = []

    # Read data from the file, skipping the header
    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            values = line.strip().split('\t')
            ideogram_data.append({
                'ref_chrom': values[0],
                'ref_start': int(values[1]),
                'ref_end': int(values[2]),
                'query_chrom': values[3],
                'query_start': int(values[4]),
                'query_end': int(values[5]),
                'annotation_type': values[6]
            })

    # Calculate the maximum length for the ideogram
    max_length = max(entry['ref_end'] for entry in ideogram_data)

    # Set up the plot
    fig, ax = plt.subplots(figsize=(8, 2))
    ax.set_xlim(0, max_length)
    ax.set_ylim(0, 1)
    ax.set_yticks([])

    # Plot rectangles for each annotation
    for entry in ideogram_data:
        rect = patches.Rectangle((entry['ref_start'], 0),
                                 entry['ref_end'] - entry['ref_start'],
                                 1, linewidth=1, edgecolor='none',
                                 facecolor=get_color(entry['annotation_type']))
        ax.add_patch(rect)

    plt.show()

def get_color(annotation_type):
    # Assign a unique color for each annotation type
    color_map = {
        'INS': 'blue',
        'DEL': 'red',
        'NA': 'gray',
    }

    return color_map.get(annotation_type, 'gray')


if __name__ == "__main__":
    file_path = "/home/bloodmark/workarea/synteny/LG2/5000/lg2_blocks_coords_plotsr.txt"  # Replace with your file path
    create_ideogram(file_path)
