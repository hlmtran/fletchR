#' Identify outlier cells by total accessibility
#'
#' Identifies cells whose total counts fall outside user-specified quantiles.
#'
#' @param mat A sparse count matrix.
#' @param quantiles Numeric vector, length two, lower and upper quantiles
#' @param excludeZeros remove colsum = 0
#' 
#' @return logical vector of outlier columns.
#'
#' @export
#'
outlierByQuantile <- function(mat, quantiles=c(0.02, 0.98),excludeZeros = TRUE) {

  #from ArchR
  stopifnot(inherits(mat, "Matrix"))
  if (length(quantiles) != 2 ||
      any(quantiles < 0) ||
      any(quantiles > 1)) {
    stop("quantiles must contain two values between 0 and 1.")
  }

  colSm <- colSums(mat)
  
  if(excludeZeros){
    toTest = colSm[colSm > 0]
  } else {
    toTest = colSm 
  }
  thresholds <- stats::quantile(toTest, probs=sort(quantiles), names=FALSE)
  idxOutlier <- (colSm <= thresholds[1] | colSm >= thresholds[2])
  return(idxOutlier)

}
