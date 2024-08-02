# README
# This script reads sequence data from a text file, applies log transformation to the sequence lengths,
# creates histograms of the log-transformed lengths, and filters the data based on sequence length.
# The script performs the following steps:
# 1. Reads the sequence data from a text file.
# 2. Applies log transformation (base 10) to the sequence lengths.
# 3. Calculates basic statistics (mean, median, range) for the sequence lengths.
# 4. Creates a histogram of log-transformed sequence lengths.
# 5. Saves the histogram as a PNG image.
# 6. Filters the data for sequences with lengths above a specified threshold.
# 7. Creates a histogram of log-transformed sequence lengths for the filtered data.

# Author: Bloodmark

# Input details:
# - Input text file: "/home/bloodmark/Chr0_RagTag_contigs_corrected.hist"
#   - The file should be tab-delimited with columns: SequenceName and Length.

# Output details:
# - A histogram of log-transformed sequence lengths.
# - A filtered histogram of log-transformed sequence lengths for sequences above a specified threshold.
# - Both histograms are displayed and saved as PNG images.

# Script functioning:
# 1. Reads the sequence data from the input text file.
# 2. Applies log transformation (base 10) to the sequence lengths.
# 3. Calculates basic statistics (mean, median, range) for the sequence lengths.
# 4. Creates and displays a histogram of log-transformed sequence lengths.
# 5. Saves the histogram as a PNG image.
# 6. Filters the data for sequences with lengths above a specified threshold.
# 7. Creates and displays a histogram of log-transformed sequence lengths for the filtered data.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the input file is in the specified path and properly formatted.

# Code starts below:

# Read the data from the text file
data <- read.table("/home/bloodmark/Chr0_RagTag_contigs_corrected.hist", header = FALSE, sep = "\t", col.names = c("SequenceName", "Length"))

# Apply the log transformation (base 10)
data$LogLength <- log10(data$Length)

# Calculate basic statistics
mean_length <- mean(data$Length)
median_length <- median(data$Length)
range_length <- range(data$Length)

cat("Mean Length:", mean_length, "\n")
cat("Median Length:", median_length, "\n")
cat("Range Length:", range_length, "\n")

# Create a histogram of log-transformed sequence lengths
hist(data$LogLength, breaks = 20, main = "Log Sequence Length Histogram", xlab = "Log Sequence Length (base 10)", ylab = "Frequency", col = "blue")

# Save the plot as an image file (e.g., PNG)
png("log_sequence_length_histogram.png", width = 800, height = 600)
hist(data$LogLength, breaks = 20, main = "Log Sequence Length Histogram", xlab = "Log Sequence Length (base 10)", ylab = "Frequency", col = "blue")
dev.off()

# Filter the data for sequences with lengths above 1,000
filtered_data <- data[data$Length > 1000, ]

# Calculate the sum of lengths for the filtered data and the original data
sum_filtered_length <- sum(filtered_data$Length)
sum_original_length <- sum(data$Length)

cat("Sum of Filtered Lengths:", sum_filtered_length, "\n")
cat("Sum of Original Lengths:", sum_original_length, "\n")

# Create a histogram of log-transformed sequence lengths for filtered data
hist(filtered_data$LogLength, breaks = 20, main = "Log Sequence Length Histogram (Length > 1,000)", xlab = "Log Sequence Length (base 10)", ylab = "Frequency", col = "blue")

# Save the plot as an image file (e.g., PNG)
png("log_sequence_length_histogram_filtered.png", width = 800, height = 600)
hist(filtered_data$LogLength, breaks = 20, main = "Log Sequence Length Histogram (Length > 1,000)", xlab = "Log Sequence Length (base 10)", ylab = "Frequency", col = "blue")
dev.off()

