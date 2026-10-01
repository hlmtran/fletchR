#' scale factors (e.g. PCA) to Z-scores, copied from ArchR
#'
#'
#' @export
#'
ArchRRowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE){
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m),`/`)
  if(limit){
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}