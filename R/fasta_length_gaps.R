# Function to count the number of 'N' characters in a sequence
count_Ns <- function(sequence) {
  n_count <- sum(strsplit(sequence, NULL)[[1]] == 'N')
  return(n_count)
}

# Function to read a FASTA file and count 'N's in each sequence
count_Ns_in_fasta <- function(file_path) {
  # Read the FASTA file
  fasta_data <- readLines(file_path)
  
  # Initialize variables to store sequence information
  current_sequence <- NULL
  total_length <- 0
  
  # Iterate through the lines in the FASTA file
  for (line in fasta_data) {
    if (startsWith(line, ">")) {
      # If line starts with '>', it's a header line
      if (!is.null(current_sequence)) {
        # If there was a previous sequence, update total length
        total_length <- total_length + nchar(current_sequence)
      }
      
      # Reset current sequence for the new header
      current_sequence <- ""
    } else {
      # If the line is not a header, it's part of the sequence
      current_sequence <- paste0(current_sequence, line)
    }
  }
  
  # Update total length for the last sequence in the file
  total_length <- total_length + nchar(current_sequence)
  
  # Count 'N's in the entire sequence
  n_count <- count_Ns(current_sequence)
  
  return(list(total_length = total_length, n_count = n_count))
}

# Function to process all FASTA files in a directory
process_fasta_directory <- function(directory_path) {
  # Get a list of all files in the directory
  fasta_files <- list.files(directory_path, pattern = "\\.fasta$", full.names = TRUE)
  
  # Initialize a data frame to store results
  results <- data.frame(File = character(), TotalLength = numeric(), NCount = numeric(), stringsAsFactors = FALSE)
  
  # Iterate through each FASTA file
  for (file_path in fasta_files) {
    # Get results for the current file
    file_results <- count_Ns_in_fasta(file_path)
    
    # Extract file name without the path and extension
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Add results to the data frame
    results <- rbind(results, c(file_name, file_results$total_length, file_results$n_count))
  }
  
  return(results)
}

# Path to directory containing FASTA files
fasta_directory_path <- "/home/bloodmark/workarea/new_final_chr/chromosomes/"
result <- process_fasta_directory(fasta_directory_path)

# Print the final results
print(result)
