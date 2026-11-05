#' Remove LSI dimensions that are correlated with sequencing depth
#'
#' @param projectedMat The cell-by-dimension matrix obtained from LSI.
#' @param depth A numeric vector representing sequencing depth for each cell.
#' @param corCutOff A numeric threshold for correlation; dimensions with correlation above this value will be removed.
#' @return A filtered cell-by-dimension matrix with depth-correlated dimensions removed.
#' @export
removeDepthCorrelatedDims <- function(projectedMat, depth, corCutOff=0.75) {
  
  cors = correlateDimensions(projectedMat = projectedMat, depth = depth)
  message("Checking for depth-correlated columns...")
  toKeep <- which(cors[, 1] <= corCutOff)
  message("Kept ", (length(toKeep)/ncol(projectedMat)*100), "% of columns.")
  projectedMat <- projectedMat[, toKeep, drop = FALSE]
  # colnames(projectedMat) <- paste0("LSI",seq_len(ncol(projectedMat)))
  return(projectedMat)
}


#' Correlate projected dimensions with sequencing depth
#'
#' Computes absolute Pearson correlation for univariate depth or distance
#' correlation for multivariate depth.
#'
#' @param projectedMat Numeric matrix with cells in rows and dimensions in columns.
#' @param depth Numeric vector or matrix of sequencing depths.
#'
#' @return Numeric matrix with one correlation per projected dimension.
#'
#' @importFrom stats cor
#' @importFrom energy dcor
correlateDimensions = function(projectedMat, depth) {
  
  depth = as.matrix(depth)
  
  if (ncol(depth) > 1) {
    cors = apply(projectedMat, 2, function(x)
      energy::dcor(x, depth))
    cors = as.matrix(cors)
    
  } else {
    cors = abs(cor(projectedMat, depth))
                   
  }
  
  return(cors)
}