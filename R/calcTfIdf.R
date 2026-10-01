#' Core computation: Calculate log(TF-IDF), IDF vector/matrix, and outlier cell IDs
#'
#' @param mat A sparse matrix (features x cells).
#' @param outlierQuantiles Numeric vector of length two or NULL.
#' @param excludeZeros Logical indicating whether to exclude cells with zero total counts when identifying outliers.
#' @param scaleTo Numeric scalar; the scale factor for the logTFIDF transformation.
#' @param ... Additional arguments passed to logTFIDF().
#'
#' @return A list containing:
#'   \item{tfidf}{The logTFIDF sparse matrix}
#'   \item{idf}{The computed IDF values}
#'   \item{outliers}{Character vector of outlier column names}
calcTfIdf <- function(mat, outlierQuantiles = c(0.02, 0.98), excludeZeros = TRUE, scaleTo=10000) {
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices")

  if (!is.null(outlierQuantiles)) {
    idxOutliers <- outlierByQuantile(mat, outlierQuantiles, excludeZeros = excludeZeros) # contains both outliers and 0 columns
    idx0ColSum <- colSums(mat) == 0
    idxOutliers <- idxOutliers & !idx0ColSum # keep only outliers
  } else {
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