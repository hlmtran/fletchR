#' Compute inverse document frequencies
#'
#' Computes the inverse document frequency (IDF) for each feature using the
#' definition:
#'
#' \deqn{
#' \mathrm{IDF} = \frac{N}{n_i}
#' }
#'
#' where \eqn{N} is the number of cells and \eqn{n_i} is the number of cells in
#' which feature \eqn{i} is observed.
#'
#' @param mat A sparse count matrix.
#'
#' @return A sparse vector containing one IDF value per feature.
#'
#' @export
#'
getIDF <- function(mat){

  message("Computing IDF (inverse document frequency) table... ", appendLF=0)
  # idf <- as(((ncol(mat)+1)/(rowSm2+1)) + 1, "sparseVector")
  idf = as(ncol(mat)/rowSums(binarizeMat(mat)), "sparseVector")
  attr(idf, 'names') <- rownames(mat)
  message("done.")
  return(idf)

}
