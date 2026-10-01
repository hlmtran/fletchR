#' Core computational engine: Compute LSI on a matrix
#'
#' @param mat Feature-by-cell matrix (e.g. TF-IDF).
#' @param nDimensions Number of singular values to calculate.
#' @param outliers Logical vector indicating outlier cells to ignore during SVD fitting.
#' @param scaleDims Logical; whether to Z-score LSI dimensions across cells.
#' @return A cell-by-dimension projected matrix.
calcLSI <- function(mat, 
                    nDimensions = 30, 
                    outliers = NULL, 
                    scaleDims = FALSE) {
  if (is.null(outliers)) {
    outliers <- rep(FALSE, ncol(mat))
  }
  message("Running SVD...")
  mat = mat[rowSums(mat[,!outliers])>0,]
  svd <- irlba::irlba(mat[,!outliers], nDimensions, nDimensions)
    
  matSVD = projectSVD(mat,svd$u,svd$d,nDimensions)
  if (scaleDims) {
    message("Scaling matSVD...")
    matSVD <- rowZscores(matSVD)
  } 
  
  return(matSVD)
}

