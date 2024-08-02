#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 18:54:01 2023

@author: bloodmark
"""

# README
# This script processes contig data to identify and plot contigs potentially associated with chromosome Y.
# The script performs the following steps:
# 1. Reads and processes contig data from a CSV file.
# 2. Filters contigs based on length and coverage.
# 3. Plots the proportion shared with females vs. coverage.
# 4. Identifies contigs that meet specific criteria for chromosome Y.
# 5. Writes the names of potential chromosome Y contigs to a file.

# Author: Bloodmark

# Input details:
# - Input CSV file: "/media/bloodmark/HDD6_SS_extra/DiscoverY/annotated_male_contigs/01B/01B_male_contigs_names.txt"
#   - The file should contain contig data with columns: Chr_name, Length, Proportion, Coverage.

# Output details:
# - Output plot: "01B_discoverY.png"
#   - A scatter plot of proportion shared with females vs. coverage.
# - Output text file: "01B_y_contig_names.txt"
#   - A file containing the names of contigs potentially associated with chromosome Y.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set the plot size
plt.rcParams['figure.figsize'] = (12, 8)

# Load the contig data
ctgs_all = pd.read_csv('/media/bloodmark/HDD6_SS_extra/DiscoverY/annotated_male_contigs/01B/01B_male_contigs_names.txt', names=["Chr_name", "Length", "Proportion", "Coverage"], delimiter=' ', header=None)
ctgs_mapped = ctgs_all[ctgs_all['Chr_name'] != "unmapped"]
print(ctgs_mapped.tail())
print(ctgs_mapped.shape)

# Filter contigs based on length
ctgs = ctgs_mapped[ctgs_mapped['Length'] > 100]
print(ctgs.shape)

# Flip the proportions (Proportions from now on will be proportion shared with female)
ctgs.loc[:, "Proportion"] = ctgs["Proportion"].apply(lambda x: 1 - x)
print("Checking dataframe...")
print(ctgs.head())

# Add a Label column which is 1 if not chrY, 0 if chrY
ctgs.loc[:, "Label"] = ctgs["Chr_name"].apply(lambda x: 0 if x == "chrY" else 1)
ctgs.head()

# Remove outlier contigs with super high proportion
ctgs = ctgs[ctgs['Proportion'] <= 1.0]
ctgs.shape
ctgs['Length'].sum()

# Remove outlier contigs with super high coverage
# Drops about 2.5% of all contigs
ctgs_less_500_cvg = ctgs[ctgs['Coverage'] < 1000]
ctgs_less_500_cvg.shape
print("Total length of all contigs after removing outliers: ", ctgs_less_500_cvg['Length'].sum())

# Print total length of Y contigs after removing outliers
non_outlier_Y_ctgs = ctgs_less_500_cvg[ctgs_less_500_cvg['Chr_name'] == "chrY"]
print("Total length of Y contigs after removing outliers: ", non_outlier_Y_ctgs['Length'].sum())
print("# of Y ctgs: ", non_outlier_Y_ctgs.shape)

# Extract 'Proportion' and 'Coverage' columns
proportion = ctgs_less_500_cvg['Proportion']
coverage = ctgs_less_500_cvg['Coverage']

# Define conditions for labeling points green
green_condition = (proportion < 0.2) & (coverage >= 10) & (coverage <= 200)

# Create a scatter plot for non-Y contigs (blue points)
plt.scatter(proportion[~green_condition], coverage[~green_condition], c='blue', marker='o', label='Non-Y Contigs')

# Create a scatter plot for Y contigs (green points)
plt.scatter(proportion[green_condition], coverage[green_condition], c='green', marker='o', label='Y Contigs')

# Set labels and title
plt.xlabel('Proportion')
plt.ylabel('Coverage')
plt.title('Proportion vs Coverage for 01B')

# Set axis limits
plt.xlim(0, 1)  # Adjust these values based on your data
plt.ylim(0, 1000)  # Adjust these values based on your data

# Show the legend
plt.legend()

# Save the plot with all labels into a file
plt.savefig("01B_discoverY.png", bbox_inches='tight')

# Show the plot
plt.show()

# Filter the DataFrame based on the conditions
filtered_ctgs = ctgs_less_500_cvg[(ctgs_less_500_cvg['Proportion'] < 0.2) & (ctgs_less_500_cvg['Coverage'] >= 10) & (ctgs_less_500_cvg['Coverage'] <= 200)]

# Extract the 'chr_name' values from the filtered DataFrame
chr_names = filtered_ctgs['Chr_name']

# Write the 'chr_name' values to a file named "01B_y_contig_names.txt"
with open("01B_y_contig_names.txt", "w") as file:
    for name in chr_names:
        file.write(name + '\n')

print("Files have been saved.")

