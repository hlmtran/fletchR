#' Identify low accessibility cells
#'
#' Determines which cells observe fewer than a specified number of features
#'
#' @param mat A sparse count matrix.
#' @param prune Minimum number of features observed within a cell.
#'
#' @return A logical vector indicating which columns should be removed.
#'
#' @export
#'
pruneCols = function(mat, prune=1){
  
  cSum <- colSums(mat > 0)  # doesn't binarization already do this?!
  
  toPrune <- (cSum < prune)
  return(toPrune)
  
}