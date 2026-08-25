#' Compute term-frequency normalization
#'
#' Computes term frequencies by dividing each nonzero entry by the total count
#' of its corresponding column.
#'
#' @param mat A sparse count matrix.
#'
#' @return
#' A sparse matrix containing term-frequency normalized values.
#'
#' @export
#'
getTF <- function(mat){

  # directly borrowed from ArchR
  colSm <- Matrix::colSums(mat)
  colSm[colSm==0] = 1
  # rowSm <- Matrix::rowSums(mat)
  message("Term frequency normalization... ", appendLF=FALSE)
  mat@x <- (mat@x / rep.int(colSm, Matrix::diff(mat@p)))
  message("done.")
  return(mat)

}
