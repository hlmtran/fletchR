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
  
  idxOutlier = detectQuantileOutliers(v = colSm,quantiles = quantiles,excludeZeros=excludeZeros)
  return(idxOutlier)
  
}

detectQuantileOutliers = function(v,quantiles,excludeZeros = TRUE,names=FALSE){
  if(excludeZeros){
    toTest = v[v > 0]
  } else {
    toTest = v 
  }
  thresholds <- stats::quantile(toTest, probs=sort(quantiles), names=names)
  idxOutlier <- (v <= thresholds[1] | v >= thresholds[2])
  
  if(excludeZeros){idxOutlier[v == 0] = FALSE }
  return(idxOutlier)
}