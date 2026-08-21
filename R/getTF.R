#' Compute term-frequency normalization
#'
#' Computes term frequencies by dividing each nonzero entry by the total count
#' of its corresponding column.
#'
#' @param mat A sparse count matrix.
#' @param keep0col Logical indicating to keep columns with colSum = 0
#'
#' @return
#' A sparse matrix containing term-frequency normalized values.
#'
#' @export
#'
getTF <- function(mat,keep0col){

  # directly borrowed from ArchR
  colSm <- Matrix::colSums(mat)
  if(keep0col){
    colSm[colSm==0] = 1
    }
  # rowSm <- Matrix::rowSums(mat)
  message("Term frequency normalization... ", appendLF=FALSE)
  mat@x <- (mat@x / rep.int(colSm, Matrix::diff(mat@p)))
  message("done.")
  return(mat)

}
