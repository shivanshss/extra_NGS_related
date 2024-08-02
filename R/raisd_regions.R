# README
# This script identifies outlier windows based on Mu values from a RAiSD report, merges nearby outlier windows,
# and saves the merged outlier regions and the top 20 outlier regions based on Mu value to separate files.
# The script uses the dplyr and tidyr libraries for data manipulation.

# Author: Bloodmark

# Input details:
# - Input file: '/media/bloodmark/HDD6_SS_extra/chapter2/w_newref_old/13.RAiSD/Report/RAiSD_Report.my25.NC_007425.3'
#   - The file should be a tab-delimited text file with columns: serial_number, start, end, sweep, SFS, LD, Mu

# Output details:
# - Merged outlier regions file: 'LG10_merged_outlier_regions.txt'
# - Top 20 merged outlier regions file: 'LG10_top_20_merged_outlier_regions.txt'

# Script functioning:
# 1. Load necessary libraries.
# 2. Read the input data file.
# 3. Define a function to identify outliers based on a specified threshold (in standard deviations above the mean Mu value).
# 4. Identify outliers using the defined threshold.
# 5. Sort the identified outliers by Mu value in descending order.
# 6. Define a function to merge nearby windows based on a specified distance threshold.
# 7. Merge nearby outlier windows using the defined distance threshold.
# 8. Print the merged outlier regions.
# 9. Extract and print the top 20 merged outlier regions based on Mu value.
# 10. Save the merged outlier regions and the top 20 merged outlier regions to separate files.

# Example of running the script:
# Rscript your_script_name.R

# Note: Ensure the dplyr and tidyr libraries are installed and accessible.

# Code starts below:

# Load necessary libraries
library(dplyr)
library(tidyr)

# Read the data (assuming the data is in a tab-delimited file called 'data.txt')
data <- read.table('/media/bloodmark/HDD6_SS_extra/chapter2/w_newref_old/13.RAiSD/Report/RAiSD_Report.my25.NC_007425.3', header = FALSE, col.names = c('serial_number', 'start', 'end', 'sweep', 'SFS', 'LD', 'Mu'))

# Define a function to identify outliers based on a threshold
identify_outliers <- function(data, threshold) {
  mean_mu <- mean(data$Mu)
  sd_mu <- sd(data$Mu)
  
  # Identify outliers
  outliers <- data %>% filter(Mu > (mean_mu + threshold * sd_mu))
  return(outliers)
}

# Set a threshold for identifying outliers (e.g., 2 standard deviations above the mean)
threshold <- 2

# Identify outlier windows
outliers <- identify_outliers(data, threshold)

# Sort outliers by Mu value in descending order
outliers <- outliers %>% arrange(desc(Mu))

# Function to merge nearby windows
merge_nearby_windows <- function(data, distance_threshold) {
  # Initialize an empty list to store merged regions
  merged_regions <- list()
  current_region <- data[1,]
  
  for (i in 2:nrow(data)) {
    if (data$start[i] - current_region$end[nrow(current_region)] <= distance_threshold) {
      # If the current window is close enough to the previous region, merge it
      current_region <- rbind(current_region, data[i,])
    } else {
      # If the current window is not close enough, save the previous region and start a new one
      merged_regions <- append(merged_regions, list(current_region))
      current_region <- data[i,]
    }
  }
  # Append the last region
  merged_regions <- append(merged_regions, list(current_region))
  
  # Combine all merged regions into a single data frame
  merged_regions_df <- do.call(rbind, merged_regions)
  return(merged_regions_df)
}

# Set a distance threshold for merging nearby windows (you can adjust this value)
distance_threshold <- 10000  # for example, 10,000 base pairs

# Merge nearby outlier windows
merged_outliers <- merge_nearby_windows(outliers, distance_threshold)

# Print the merged outlier regions
print(merged_outliers)

# Get the top 20 rows based on Mu value
top_20_outliers <- head(merged_outliers, 20)

# Print the top 20 merged outlier regions
print(top_20_outliers)

# Save the merged outlier regions to a file
write.table(merged_outliers, 'LG10_merged_outlier_regions.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(top_20_outliers, 'LG10_top_20_merged_outlier_regions.txt', row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)

