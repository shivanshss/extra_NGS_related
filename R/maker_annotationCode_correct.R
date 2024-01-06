# Install and load the dplyr package if not already installed
#if (!requireNamespace("dplyr", quietly = TRUE)) {
#  install.packages("dplyr")
#}

library(dplyr)

# Read the TSV file
data <- read.csv("maker_annotation_col_corr.gff", header = FALSE)

# Remove duplicate rows
data <- distinct(data)

data <- mutate(data, V6 = V5 - V4)


# Get unique values in column 2 and column 3
unique_values_col2 <- unique(data$V2)
unique_values_col3 <- unique(data$V3)

# Initialize an empty data frame for the summary
summary_table <- data.frame(Type = character(), Count = numeric(), stringsAsFactors = FALSE)

# Get unique values in column 3
unique_values <- unique(data$V3)

# Create separate files for each unique value
for (value in unique_values) {
  subset_df <- data[data$V3 == value, ]
  write.table(subset_df, file = paste0(value, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}


# Loop through unique values in column 2 and column 3
for (col2_value in unique_values_col2) {
  for (col3_value in unique_values_col3) {
    count <- sum(data$V2 == col2_value & data$V3 == col3_value)
    summary_table <- rbind(summary_table, data.frame(Type = col2_value, Value = col3_value, Count = count))
  }
}


# Print the summary table

print(summary_table)

