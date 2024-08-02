# README
# This script visualizes pairwise synteny alignments using the syntenyPlotteR package. It performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Sets the working directory to the location of the input files.
# 3. Draws pairwise synteny plots using the draw.pairwise function from the syntenyPlotteR package.

# Author: Bloodmark

# Input details:
# - Input files should be located in the specified directory: "/home/bloodmark/Downloads/syntenyPlotteR/data"
# - Required input files:
#   - "lengths.txt": A file containing lengths of the chromosomes or contigs.
#   - "alignment_1.txt": A file containing the first alignment data.
#   - "alignment_2.txt": (Optional) A file containing the second alignment data for a second pairwise plot.

# Output details:
# - Pairwise synteny plots saved as PNG files:
#   - "trial_w_3.png": Synteny plot with two alignment files.
#   - "trial_w_2.png": Synteny plot with one alignment file.

# Script functioning:
# 1. Installs and loads the necessary packages (syntenyPlotteR, ggplot2, tidyr, tidyverse, stringr).
# 2. Sets the working directory to the location of the input files.
# 3. Draws pairwise synteny plots using the draw.pairwise function.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input files are in the specified paths and properly formatted.

# Code starts below:

# Install required packages if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("marta-fb/syntenyPlotteR")
install.packages("tidyverse")
install.packages("stringr")

# Load required libraries
library(syntenyPlotteR)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)

# Update all packages without asking for confirmation
update.packages(ask = FALSE)

# Get the current working directory
cat("Current working directory:", getwd(), "\n")

# Set the working directory to the location of the input files
setwd("/home/bloodmark/Downloads/syntenyPlotteR/data")

# Draw pairwise synteny plots
draw.pairwise("trial_w_3", "lengths.txt", "alignment_1.txt", "alignment_2.txt", fileformat = "png", w = 13, h = 5)
draw.pairwise("trial_w_2", "lengths.txt", "alignment_1.txt", fileformat = "png", w = 13, h = 5)

