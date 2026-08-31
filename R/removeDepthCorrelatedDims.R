#' Remove LSI dimensions that are correlated with sequencing depth
#'
#' @param projectedMat The cell-by-dimension matrix obtained from LSI.
#' @param depth A numeric vector representing sequencing depth for each cell.
#' @param corCutOff A numeric threshold for correlation; dimensions with correlation above this value will be removed.
#' @return A filtered cell-by-dimension matrix with depth-correlated dimensions removed.
#' @export
removeDepthCorrelatedDims <- function(projectedMat, depth, corCutOff=0.75) {
  message("Checking for depth-correlated columns...")
  toKeep <- which(abs(cor(projectedMat, depth)[, 1]) <= corCutOff)
  message("Kept ", (length(toKeep)/ncol(projectedMat)*100), "% of columns.")
  projectedMat <- projectedMat[, toKeep, drop = FALSE]
  # colnames(projectedMat) <- paste0("LSI",seq_len(ncol(projectedMat)))
  return(projectedMat)
}
