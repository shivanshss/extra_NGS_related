# README
# This script uses the syntenyPlotteR package to reformat synteny data and generate visualizations.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Reformats synteny data from a MAF file.
# 3. Generates various types of visualizations using the reformatted data.

# Author: Bloodmark

# Input details:
# - Input MAF file: "/home/bloodmark/workarea/synteny/LG2/sibeliaz_out/alignment.maf"
#   - The file should be a Multiple Alignment Format (MAF) file.

# Output details:
# - Reformatted synteny data file: "reformat_lg2_maf" in the specified directory.
# - Various visualizations in PNG format.

# Script functioning:
# 1. Install and load the devtools and syntenyPlotteR packages.
# 2. Reformat the synteny data from the input MAF file.
# 3. Generate visualizations using the reformat.syntenyData and draw functions from syntenyPlotteR.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure internet connectivity is available for installing the necessary packages.

# Code starts below:

# Install and load the devtools package if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

library(devtools)

# Install and load the syntenyPlotteR package from GitHub
devtools::install_github("Farre-lab/syntenyPlotteR")
library(syntenyPlotteR)

# Reformat synteny data from the MAF file
reformat.syntenyData("/home/bloodmark/workarea/synteny/LG2/sibeliaz_out/alignment.maf", "reformat_lg2_maf", directory = "/home/bloodmark/workarea/synteny/LG2/")

# Variables for visualization functions
output <- "output"
chrRange <- "chrRange"  # Placeholder: replace with actual chromosome range
data_file <- "data_file"  # Placeholder: replace with actual data file path
sizefile <- "sizefile"  # Placeholder: replace with actual size file path
file_data <- "file_data"  # Placeholder: replace with actual file data path
colours.default <- c("lightblue", "lightpink")  # Example default colours, modify as needed

# Evolution Highway visualization
draw.eh(output, chrRange, data_file, fileformat = "png", colour = "lightblue", inverted.colour = "lightpink", w = 5.5, h = 10, ps = 10)

# Ideogram visualization
draw.ideogram(file_data, sizefile, output, fileformat = "png", colours = colours.default, w = 8.5, h = 10, ps = 5)

# Linear visualization
draw.linear(output, sizefile, fileformat = "png", colours = colours.default, w = 13, h = 5, opacity = 0.5)

