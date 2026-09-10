#Extracted from ArchR

findVariableFeaturesByCluster <- function(mat , clusters, nFeatures = 25000, scaleTo = 10000, method = "var") {
  stopifnot(inherits(mat, "Matrix"))
  stopifnot(requireNamespace("matrixStats", quietly = TRUE))
  
  # 1. Fast Pseudobulking via Sparse Matrix Multiplication
  cluster_factor <- as.factor(clusters)
  indicatorMat <- Matrix::sparseMatrix(
    i = 1:ncol(sparseMat), 
    j = as.integer(cluster_factor), 
    x = 1,
    dims = c(ncol(sparseMat), nlevels(cluster_factor))
  )
  
  groupMat <- mat %*% indicatorMat
  
  # 2. Log-Normalize (CP10K)
  groupMat <- t(t(groupMat) / colSums(groupMat)) * scaleTo
  groupMat <- log2(groupMat + 1)
  
  # 3. Calculate Variance
  if (method == "var") {
    # Log-variability
    vars <- matrixStats::rowVars(groupMat)
  } else if (method == "vmr") {
    # Variance-to-mean ratio
    vars <- matrixStats::rowVars(groupMat) / rowMeans(groupMat)
  } else {
    stop("Method must be 'var' or 'vmr'")
  }
  
  # 4. Select Top Features
  ranked_idx <- order(vars, decreasing = TRUE)
  top_idx <- head(ranked_idx, nFeatures)
  
  # Return the subsetted sparse matrix with only the highly variable features
  return(mat[sort(top_idx), , drop = FALSE])
}