# README
# This script visualizes the relationship between geographic distance and genetic distance using a scatter plot,
# fits a linear model to the data, and performs a Mantel test to statistically assess the correlation between
# geographic and genetic distances. The script uses the ggplot2 and vegan libraries for plotting and statistical analysis.

# Author: Bloodmark

# Input details:
# - Input file: "/home/bloodmark/genetic_distance_vs_geo_distance.txt"
#   - The file should contain at least two columns: Distance (geographic distance) and Absolute.divergence (genetic distance)

# Output details:
# - Scatter plot with a fitted linear model line
# - Mantel test result printed to the console

# Script functioning:
# 1. Load necessary libraries.
# 2. Load the data from a specified file.
# 3. Compute the logarithm (base 2) of the geographic distance.
# 4. Create a scatter plot of log2 geographic distance against genetic distance with a fitted linear model line.
# 5. Prepare distance matrices for the Mantel test.
# 6. Perform the Mantel test and print the result.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the ggplot2 and vegan libraries are installed and accessible.

# Load necessary libraries
library(ggplot2)
library(vegan)

# Load data
data <- read.table("/home/bloodmark/genetic_distance_vs_geo_distance.txt", header = TRUE)

# Compute log2 of the geographic distance
data$log_distance <- log2(data$Distance)

# Plot using ggplot2
plot <- ggplot(data, aes(x = log_distance, y = Absolute.divergence)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "Log2 Geographic Distance (km)", y = "Genetic Distance (dxy)",
       title = "Isolation by Distance Plot")
print(plot)

# Prepare distance matrices for Mantel test
geo_dist_matrix <- as.dist(matrix(data$Distance, nrow = sqrt(length(data$Distance))))
gen_dist_matrix <- as.dist(matrix(data$Absolute.divergence, nrow = sqrt(length(data$Absolute.divergence))))

# Perform Mantel test
mantel_result <- mantel(geo_dist_matrix, gen_dist_matrix, method="spearman")
print(mantel_result)

