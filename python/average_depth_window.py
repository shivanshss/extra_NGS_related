# README
# This script calculates the average sequencing depth in 1kb windows from a depth information file.
# The script performs the following steps:
# 1. Reads the depth information from a file.
# 2. Calculates the average depth in 1kb windows.
# 3. Saves the calculated average depths to a new file.

# Author: Bloodmark

# Input details:
# - Input file: 'depth_info.ldepth'
#   - The file should contain columns: CHROM, POS, N_SITES, DEPTH.

# Output details:
# - Output file: 'average_depth_1kb_windows.csv'
#   - The file will contain columns: CHROM, START, END, AVG_DEPTH.

# Script functioning:
# 1. Reads the depth information from the input file.
# 2. Calculates the average depth in 1kb windows.
# 3. Saves the calculated average depths to a new file.

# Example of running the script:
# python your_script_name.py

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

import pandas as pd
import numpy as np

# Load depth information
depth_data = pd.read_csv('depth_info.ldepth', sep='\t', header=0)
depth_data.columns = ['CHROM', 'POS', 'N_SITES', 'DEPTH']

# Initialize window size
window_size = 1000

# Function to calculate average depth in windows
def calculate_average_depth(depth_data, window_size):
    # Initialize variables
    window_starts = np.arange(1, depth_data['POS'].max(), window_size)
    avg_depths = []

    for start in window_starts:
        end = start + window_size - 1
        window_data = depth_data[(depth_data['POS'] >= start) & (depth_data['POS'] <= end)]
        if not window_data.empty:
            avg_depth = window_data['DEPTH'].mean()
            avg_depths.append([window_data['CHROM'].iloc[0], start, end, avg_depth])
        else:
            avg_depths.append([depth_data['CHROM'].iloc[0], start, end, 0])

    return pd.DataFrame(avg_depths, columns=['CHROM', 'START', 'END', 'AVG_DEPTH'])

# Calculate average depth in 1kb windows
avg_depth_df = calculate_average_depth(depth_data, window_size)

# Save the result
avg_depth_df.to_csv('average_depth_1kb_windows.csv', index=False)

print("Average depth calculation completed and saved to 'average_depth_1kb_windows.csv'.")

