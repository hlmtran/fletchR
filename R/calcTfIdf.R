#' Core computation: Calculate log(TF-IDF), IDF vector/matrix, and outlier cell IDs
#'
#' @param mat A sparse matrix (features x cells).
#' @param idxOutliers Indices of outlier cells
#' @param scaleTo Numeric scalar; the scale factor for the logTFIDF transformation.
#' @param ... Additional arguments passed to logTFIDF().
#'
#' @return A list containing:
#'   \item{tfidf}{The logTFIDF sparse matrix}
#'   \item{idf}{The computed IDF values}
#'   \item{outliers}{Character vector of outlier column names}
calcTfIdf <- function(mat,idxOutliers = NULL, scaleTo=10000) {
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices")
  
  if (is.null(idxOutliers)){
    idxOutliers <- logical(ncol(mat))
  }
  
  outliers <- colnames(mat)[idxOutliers]
  
  idfMat <- mat[, !idxOutliers, drop = FALSE]
  
  mat <- getTF(mat)
  idfMat <- getIDF(idfMat)
  
  tfidfMat <- logTFIDF(tf = mat, idf = idfMat, scaleTo = scaleTo)
  
  return(list(
    tfidf = tfidfMat,
    idf = idfMat,
    outliers = outliers
  ))
}