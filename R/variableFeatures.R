#Extracted from ArchR
#' Find variable features by cluster
#' 
#' @param mat       a sparse matrix (dgCMatrix) of features x cells
#' @param clusters  a vector of cluster assignments for each cell
#' @param nFeatures the number of variable features to return
#' @param scaleTo   the scale factor for log-normalization (default is 10000)
#' @param method    the method to use for calculating variability ("var" or "vmr")
#' @return          a sparse matrix (dgCMatrix) of the top variable features x cells
#' 
#' @import sparseMatrixStats
findVariableFeaturesByCluster <- function(mat , clusters, nFeatures = 25000, scaleTo = 10000, method = "var") {
  stopifnot(inherits(mat, "Matrix"))
  
  
  # 1. Fast Pseudobulking via Sparse Matrix Multiplication
  cluster_factor <- as.factor(clusters)
  indicatorMat <- Matrix::sparseMatrix(
    i = seq_len(ncol(mat)), 
    j = as.integer(cluster_factor), 
    x = 1,
    dims = c(ncol(mat), nlevels(cluster_factor))
  )
  
  groupMat <- mat %*% indicatorMat
  
  # 2. Log-Normalize (CP10K)
  groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
  groupMat <- log2(groupMat + 1) |> as("dgCMatrix")
  
  # 3. Calculate Variance
  if (method == "var") {
    # Log-variability
    vars <- sparseMatrixStats::rowVars(groupMat)
  } else if (method == "vmr") {
    # Variance-to-mean ratio
    vars <- sparseMatrixStats::rowVars(groupMat) / sparseMatrixStats::rowMeans2(groupMat)
    vars[is.na(vars) | is.nan(vars)] <- 0
  } else {
    stop("Method must be 'var' or 'vmr'")
  }
  
  # 4. Select Top Features
  ranked_idx <- order(vars, decreasing = TRUE)
  top_idx <- head(ranked_idx, nFeatures)
  
  # Return the subsetted sparse matrix with only the highly variable features
  return(rownames(mat)[sort(top_idx)])
}