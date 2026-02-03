###########################################################################

# First, parse and re-format the BirdNET embedding files to a format read by the Tensorflow visualizer.

###########################################################################

# load required packages
#install.packages("readr")  # if not already installed
library(readr)

# load in the BirdNET embeddings and take a look at them
embeddings_file <- "../Embeddings/embeddings.csv"
embeddings <- read_csv(embeddings_file, show_col_types = FALSE)

head(embeddings)

# extract the numerical part into a matrix, and save to a txt file
numbers <- embeddings$embedding

# split comma-separated strings and convert to numeric
numbers_mat <- do.call(
  rbind,
  lapply(numbers, function(num) as.numeric(strsplit(num, ",")[[1]]))
)

# save to tab-separated file
write.table(
  numbers_mat,
  file = "embedding_values_R.tsv",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# parse and reformat the metadata
parseMetadata <- function(filename) {
  split_path <- strsplit(filename, "/")[[1]]
  split_filename <- strsplit(split_path[10], "_")[[1]]
  
  category <- substr(split_filename[length(split_filename)],
                     1,
                     nchar(split_filename[length(split_filename)]) - 4)
  day <- split_filename[1]
  channel <- split_filename[3]
  time <- substr(split_filename[4],
                  1,
                  nchar(split_filename[4]) - 1)
  
  return(c(category, day, channel, time))
}

# apply parsing function
metadata_mat <- t(sapply(embeddings$file_path, parseMetadata))
metadata <- as.data.frame(metadata_mat, stringsAsFactors = FALSE)

# # set column names
colnames(metadata) <- c("Category", "Day", "Channel", "Time")

# save to tab-separated file
write.table(
  metadata,
  file = "metadata_R.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


###########################################################################

# Next, perform dimensionality reduction and visualization.

###########################################################################

library(Rtsne)
library(umap)

# read in embeddings and metadata
numbers_mat <- read.table("embedding_values_R.tsv", header = FALSE)
metadata <- read.table("metadata_R.tsv", header = TRUE)

# PCA
projection_pca <- prcomp(numbers_mat, center = TRUE, scale. = TRUE)$x[, 1:2]

# t-SNE
set.seed(123)
projection_tsne <- Rtsne(numbers_mat, dims = 2)$Y

# UMAP
projection_umap <- umap(numbers_mat)$layout

methods <- c("PCA", "tSNE", "UMAP")
projections <- list(
  PCA = projection_pca,
  tSNE = projection_tsne,
  UMAP = projection_umap
)

png(filename = "projections_R.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3))

for (m in methods) {
  proj <- projections[[m]]
  plot(
    proj[, 1], proj[, 2],
    col = "black", pch = 16,
    xlab = "Dimension 1",
    ylab = "Dimension 2",
    main = paste(m, "Projection")
  )
}

dev.off()

colors <- c(
  noise = "navy",
  wct   = "orange",
  sar   = "springgreen"
)

png(filename = "projections_by_class_R.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1, 3))

for (m in methods) {
  proj <- projections[[m]]
  
  plot(
    proj[, 1], proj[, 2],
    type = "n",
    xlab = "Dimension 1",
    ylab = "Dimension 2",
    main = paste(m, "Projection")
  )
  
  for (cat in names(colors)) {
    idx <- metadata$Category == cat
    points(
      proj[idx, 1], proj[idx, 2],
      col = colors[cat],
      pch = 8   # star marker, similar to matplotlib '*'
    )
  }
  
  legend(
    "topright",
    legend = names(colors),
    col = colors,
    pch = 8,
    bty = "n"
  )
}

dev.off()