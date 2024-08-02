# README
# This script creates a linear chromosome plot from structural variant data using ggplot2.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Reads the structural variant data from a file.
# 3. Cleans and processes the data.
# 4. Creates a linear chromosome plot using ggplot2.

# Author: Bloodmark

# Input details:
# - Structural variant data file: "/home/bloodmark/workarea/20231025/LG5/align_out_indels.txt"
#   - The file should be tab-delimited with columns: Scaffold, Type, start, end, and Size.

# Output details:
# - A linear chromosome plot created using ggplot2.

# Script functioning:
# 1. Install and load the necessary packages (if not already installed).
# 2. Read the structural variant data from the input file.
# 3. Clean and process the data by removing rows with NA values and ensuring Size is greater than or equal to 1.
# 4. Create a linear chromosome plot using ggplot2 with custom color palette and size scale.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
BiocManager::install("Biostrings")

install.packages("ggplot2")
install.packages("colorspace")
install.packages("circlize")
install.packages("RColorBrewer")

# Load required libraries
library(karyoploteR)
library(Biostrings)
library(ggplot2)
library(colorspace)
library(circlize)
library(RColorBrewer)

# Read structural variant data from file
sv_data <- read.table("/home/bloodmark/workarea/20231025/LG5/align_out_indels.txt", header = FALSE, sep = "\t", fill = TRUE)

# Display the first few rows of the data
head(sv_data)

# Set column names
colnames(sv_data) <- c("Scaffold", "Type", "start", "end", "Size")

# Remove rows where Size is less than 1 and remove rows with NA values
sv_data <- sv_data[sv_data$Size >= 1, ]
sv_data <- na.omit(sv_data)

# Determine the number of unique values in sv_data$Type
num_unique_types <- length(unique(sv_data$Type))

# Choose a custom color palette with dark, easily visible colors
custom_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")

# Create a linear chromosome plot
chromosome_plot <- ggplot(sv_data, aes(x = start, xend = end, y = Type, color = Type)) +
  geom_segment(aes(yend = Type, size = Size)) +
  scale_color_manual(values = custom_palette) +
  scale_size_continuous(range = c(1, 10)) +
  theme_minimal() +
  labs(x = "Position on Chromosome", y = "Feature Type") +
  scale_y_discrete(limits = unique(sv_data$Type))  # Set a common y-axis scale

# Display the plot
print(chromosome_plot)

