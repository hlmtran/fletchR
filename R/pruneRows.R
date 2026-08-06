#' Identify low-frequency features
#'
#' Determines which features occur in fewer than a specified number of cells.
#'
#' @param mat A sparse count matrix.
#' @param prune Minimum number of cells in which a feature must be observed.
#'
#' @return A logical vector indicating which rows should be removed.
#'
#' @export
#'
pruneRows = function(mat, prune=1){

  rSum <- rowSums(mat > 0)  # doesn't binarization already do this?!
 
  toPrune <- (rSum < prune)
  # if (sum(toPrune) > 0) {
  #   message("Minimum observed feature frequency: ", min(rSum[rSum > 0]))
  #   message("Minimum allowable feature frequency: ", prune)
  #   message("Removing ",sum(toPrune)," of ",length(toPrune)," features.")
  #   # keep <- as(Matrix::Diagonal(x=!toPrune), "sparseMatrix")
  #   # rownames(keep) <- rownames(mat) 
  #   # mat <- keep %*% mat
  #   mat = mat[!toPrune,]
  #   # origRows = origRows[!toPrune]
  #   rSum <- rowSums(mat > 0)
  #   message("Minimum feature frequency is ", min(rSum[rSum > 0]),
  #           " (", sum(!toPrune), "/", length(toPrune), " features retained).")
  # }
  return(toPrune)

}
