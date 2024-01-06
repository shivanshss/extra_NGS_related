# Read the input TSV file
data <- read.table(file="/home/bloodmark/workarea/synteny/LG2/5000/blocks_coords_cleaned1.txt", header = FALSE)
ref_name="lg2_new"
query_name="lg2_old"

# Read input file path, reference name, and query name from command line arguments
#args <- commandArgs(trailingOnly = TRUE)
#input_file <- args[1]
#ref_name <- args[2]
#query_name <- args[3]

# Read the input TSV file
#data <- read.table(input_file, header = FALSE)

# Initialize vectors to store the new columns
ref_chr_name <- character()
ref_start_pos <- integer()
ref_end_pos <- integer()
query_chr_name <- character()
query_start_pos <- integer()
query_end_pos <- integer()
annotation_type <- character()

# Process each pair of rows in the input data
for (i in seq(1, nrow(data), by = 2)) {
  ref_chr_name <- c(ref_chr_name, paste(ref_name, data[i, 2], sep = ""))
  ref_start_pos <- c(ref_start_pos, data[i, 4])
  ref_end_pos <- c(ref_end_pos, data[i, 5])
  
  query_chr_name <- c(query_chr_name, paste(query_name, data[i + 1, 2], sep = ""))
  query_start_pos <- c(query_start_pos, data[i + 1, 4])
  query_end_pos <- c(query_end_pos, data[i + 1, 5])
  
  # Determine annotation type
  ref_index <- length(ref_chr_name)
  query_index <- length(query_chr_name)
  if (data[i, 6] > data[i + 1, 6]) {
    annotation_type <- c(annotation_type, "INS")
  } else if (data[i, 6] < data[i + 1, 6]) {
    annotation_type <- c(annotation_type, "DEL")
  } else {
    annotation_type <- c(annotation_type, NA)
  }
}


# Create a data frame with the new columns
result_data <- data.frame(
  Reference_chromosome_name = ref_chr_name,
  Reference_start_position = ref_start_pos,
  Reference_end_position = ref_end_pos,
  Query_chromosome_name = query_chr_name,
  Query_start_position = query_start_pos,
  Query_end_position = query_end_pos,
  Annotation_type = annotation_type
)

# Write the result to a new TSV file
write.table(result_data, file = "lg2_blocks_coords_plotsr.txt", sep = "\t", row.names = FALSE, quote = FALSE)


