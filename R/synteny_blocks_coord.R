# README
# This script creates synteny plots using ggplot2 and circlize packages.
# The script performs the following steps:
# 1. Installs and loads the necessary packages.
# 2. Reads the synteny data from a TSV file.
# 3. Creates a linear synteny plot using ggplot2.
# 4. Creates a circular synteny plot using circlize.

# Author: Bloodmark

# Input details:
# - Input TSV file: "/home/bloodmark/workarea/synteny/LG2/5000/lg2_blocks_coords_plotsr.txt"
#   - The file should contain columns: Reference_start_position, Reference_end_position, Query_start_position, Query_end_position, and Annotation_type.

# Output details:
# - A linear synteny plot created using ggplot2.
# - A circular synteny plot created using circlize.

# Script functioning:
# 1. Install and load the ggplot2 and circlize packages.
# 2. Read the input TSV file into a data frame.
# 3. Create a linear synteny plot using ggplot2.
# 4. Define a color palette for different annotation types.
# 5. Create a circular synteny plot using circlize.
# 6. Add transparency to the grid lines and add a legend.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

# Install and load required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
library(ggplot2)

if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}
library(circlize)

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

# Read data from the TSV file
data <- read.delim("/home/bloodmark/workarea/synteny/LG2/5000/lg2_blocks_coords_plotsr.txt", header = TRUE, sep = "\t")

# Create a synteny plot using ggplot2
ggplot(data, aes(x = Reference_start_position, xend = Reference_end_position, y = Query_start_position, yend = Query_end_position)) + 
  geom_segment(color = "blue", size = 2, alpha = 0.7) + 
  labs(x = "Reference Position", y = "Query Position", title = "Synteny Plot") + 
  theme_minimal()

# Define color palette for different annotation types
colors <- setNames(c(brewer.pal(8, "Set1")[1], brewer.pal(8, "Set1")[2], "red"), c("DEL", "INS", "NA"))

# Create a circular plot
circos.initialize(factors = data$Reference_start_position, xlim = cbind(data$Reference_start_position, data$Reference_end_position))

circos.trackPlotRegion(factors = data$Reference_start_position, y = data$Query_start_position, panel.fun = function(x, y) {
  circos.rect(xleft = data$Reference_start_position, xright = data$Reference_end_position, ybottom = data$Query_start_position, ytop = data$Query_end_position, col = colors[data$Annotation_type], border = NA)
}, bg.border = NA)

# Add legend
legend("topright", legend = unique(data$Annotation_type), fill = colors[unique(data$Annotation_type)], title = "Annotation Type")

