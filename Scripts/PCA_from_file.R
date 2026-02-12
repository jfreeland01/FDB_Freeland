#' PCA from file
#'
#' Reads file with samples in columns and variables in rows, and does PCA. 
#' Writes to file scores, loadings, eigenvalues, and summary.
#'
#' @param file Filepath/filename of data matrix
#' @param center Logical, default = TRUE
#' @param scale Logical, default = FALSE
#' @param fread Logical, default = FALSE — use data.table::fread() for large input files
#'
#' @importFrom stats prcomp screeplot
#' @importFrom utils read.delim write.table
#' @importFrom data.table fread
#' @importFrom tibble column_to_rownames
#'
#' @export

library(data.table)
library(tibble)

PCA_from_file <- function(file, center = TRUE, scale = FALSE, fread = FALSE) {
  
  if (fread) {
    data <- data.table::fread(file) %>%
      tibble::column_to_rownames(var = "V1")
  } else {
    data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  }
  
  # Remove rows where all values are 0
  data <- data[rowSums(data == 0) < ncol(data), ]
  
  # Transpose so that samples are rows, genes are columns
  t.data <- t(data)
  
  # Run PCA
  pca <- prcomp(t.data, scale. = scale, center = center)
  
  # Extract outputs
  pca_scores <- cbind(Score = gsub("-", ".", rownames(pca$x)), pca$x)
  pca_loadings <- cbind(Loading = rownames(pca$rotation), pca$rotation)
  pca_evalues <- pca$sdev
  pca_summary <- data.frame(summary(pca)$importance)
  
  # Limit to first 100 PCs max
  pca_scores <- pca_scores[, 1:min(101, ncol(pca_scores)), drop = FALSE]
  pca_loadings <- pca_loadings[, 1:min(101, ncol(pca_loadings)), drop = FALSE]
  
  # Construct file names and save
  name_base <- sub(".txt", "", file)
  
  write.table(pca_scores, paste0(name_base, "_prcomp_scores.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(pca_loadings, paste0(name_base, "_prcomp_loadings.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(pca_evalues, paste0(name_base, "_prcomp_sdev.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(pca_summary, paste0(name_base, "_prcomp_summary.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Scree plot
  screeplot(pca)
}
