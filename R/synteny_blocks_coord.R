# Install and load required packages
#install.packages("ggplot2")
library(ggplot2)

# Read data from the TSV file
data <- read.delim("/media/bloodmark/HDD6_SS_extra/genome_assembly_figure_data/synteny/LG2/5000/lg2_blocks_coords_plotsr.txt", header = TRUE, sep = "\t")

# Create a synteny plot using ggplot2
ggplot(data, aes(x = Reference_start_position, xend = Reference_end_position, y = Query_start_position, yend = Query_end_position)) + geom_segment(color = "blue", size = 2, alpha = 0.7) + labs(x = "Reference Position", y = "Query Position", title = "Synteny Plot") + theme_minimal()

library(circlize)
library(RColorBrewer)

# Define color palette for different annotation types
colors <- setNames(c(brewer.pal(8, "Set1")[1], brewer.pal(8, "Set1")[2], "red"), c("DEL", "INS", "NA"))

# Find the number of sectors (genomic regions)
num_sectors <- nrow(data)

# Create a circular plot
chordDiagram(
  x = data[, c("Reference_start_position", "Reference_end_position")],
  annotationTrack = "grid",
  preAllocateTracks = 1,
  preAllocateSize = 1,
  grid.col = colors[data$Annotation_type],
  grid.border = colors[data$Annotation_type],
  track.height = 0.1,
  gap.degree = 0.01  # Adjust this value as needed
)


# Add transparency to the grid lines
for (i in seq_len(num_sectors)) {
  annotationTrack(
    x = 1,
    y = i,
    track.height = 0.2,
    name = "grid",
    col = colors[data$Annotation_type][i],
    border = colors[data$Annotation_type][i],
    transparency = 0.5  # Adjust transparency as needed
  )
}

# Add legend
legend(x = "topright", legend = unique(data$Annotation_type), fill = colors[unique(data$Annotation_type)], title = "Annotation Type")
