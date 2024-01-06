# Install required packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
install.packages("circlize")

# Load required libraries
library(GenomicRanges)
library(circlize)


# Read GFF file
gff_file <- "/home/bloodmark/workarea/sibeliaz_out/blocks_coords.gff"
gff_data <- read.table(gff_file, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

# Read MAF file
maf_file <- "/home/bloodmark/workarea/sibeliaz_out/alignment.maf"
maf_data <- readLines(maf_file)

# Extract information from MAF file
alignment_data <- lapply(strsplit(maf_data, " "), function(x) {
  if (length(x) >= 2 && x[1] == "s") {
    return(c(seqname = x[2], start = as.numeric(x[3]), end = as.numeric(x[3]) + as.numeric(x[4])))
  } else {
    print(paste("Skipping line:", paste(x, collapse = " ")))
    return(NULL)
  }
})


# Combine information from GFF and MAF files
combined_data <- lapply(seq_along(gff_data$V1), function(i) {
  gff_line <- gff_data[i, ]
  maf_line <- alignment_data[[i]]
  
  if (!is.null(maf_line)) {
    return(data.frame(seqname = gff_line[1],
                      start = gff_line[4],
                      end = gff_line[5],
                      id = gff_line[9],
                      aligned_start = maf_line[2],
                      aligned_end = maf_line[3],
                      stringsAsFactors = FALSE))
  } else {
    return(NULL)
  }
})

# Filter out NULL entries
combined_data <- combined_data[!sapply(combined_data, is.null)]

# Print the head of the GFF file
cat("Head of GFF file:\n")
cat(readLines(gff_file, n = 10), sep = "\n")
cat("\n")

# Print the head of the MAF file
cat("Head of MAF file:\n")
cat(readLines(maf_file, n = 20), sep = "\n")
cat("\n")

# Print the first few entries of combined_data
cat("First few entries of combined_data:\n")
for (i in seq_along(combined_data)) {
  cat("Index:", i, "\n")
  print(combined_data[[i]])
  cat("\n")
  if (i == 5) break  # Print only the first 5 entries for brevity
}


# Filter out NULL entries
combined_data <- combined_data[!sapply(combined_data, is.null)]

# Extract sequence names from the combined data
seqnames <- unique(unlist(lapply(combined_data, function(x) x$seqname)))

# Initialize empty vectors
start_vec <- numeric(0)
end_vec <- numeric(0)
id_vec <- character(0)
aligned_start_vec <- numeric(0)
aligned_end_vec <- numeric(0)

# Iterate over combined_data and populate vectors
combined_data <- combined_data[!sapply(combined_data, is.null)]
for (i in seq_along(combined_data)) {
  x <- combined_data[[i]]
  start_vec <- c(start_vec, as.numeric(x$start))
  end_vec <- c(end_vec, as.numeric(x$end))
  id_vec <- c(id_vec, as.character(x$id))
  aligned_start_vec <- c(aligned_start_vec, as.numeric(x$aligned_start))
  aligned_end_vec <- c(aligned_end_vec, as.numeric(x$aligned_end))
}

# Create a GRanges object with explicit sequence names
gr <- GRanges(seqnames = rep(seqnames, each = length(start_vec)),
              ranges = IRanges(start = start_vec, end = end_vec),
              id = id_vec,
              aligned_start = aligned_start_vec,
              aligned_end = aligned_end_vec)





# Print lengths for debugging
cat("Lengths:\n")
cat("seqnames:", length(seqnames), "\n")
cat("start_vec:", length(start_vec), "\n")
cat("end_vec:", length(end_vec), "\n")
cat("id_vec:", length(id_vec), "\n")
cat("aligned_start_vec:", length(aligned_start_vec), "\n")
cat("aligned_end_vec:", length(aligned_end_vec), "\n")


# Print the lengths of vectors
cat("Lengths:\n")
cat("start_vec:", length(start_vec), "\n")
cat("end_vec:", length(end_vec), "\n")
cat("id_vec:", length(id_vec), "\n")
cat("aligned_start_vec:", length(aligned_start_vec), "\n")
cat("aligned_end_vec:", length(aligned_end_vec), "\n")

# Print the head of vectors
cat("Head of vectors:\n")
cat("start_vec:", head(start_vec), "\n")
cat("end_vec:", head(end_vec), "\n")
cat("id_vec:", head(id_vec), "\n")
cat("aligned_start_vec:", head(aligned_start_vec), "\n")
cat("aligned_end_vec:", head(aligned_end_vec), "\n")




