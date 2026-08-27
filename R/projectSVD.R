#' Right project matrix
#'
#' Projects observations from a matrix into a lower-dimensional singular value
#' decomposition (SVD) space using supplied left singular vectors and singular
#' values. If the right singular vectors are not supplied, they are calculated
#' from the input matrix.
#'
#' @param mat A numeric or sparse matrix with features as rows and observations
#'   (e.g., cells) as columns.
#' @param u A numeric matrix containing the left singular vectors of the SVD.
#' @param d A numeric vector containing the singular values corresponding to
#'   the columns of `u`.
#' @param nDim Integer specifying the number of SVD dimensions to return.
#' @param v Optional numeric matrix containing the right singular vectors.
#'   If `NULL`, `v` is calculated from `mat`, `u`, and `d`. Default is `NULL`.
#'
#' @return A numeric matrix with observations as rows and SVD dimensions as
#'   columns. Row names are taken from `colnames(mat)`, and columns are named
#'   `"LSI1"`, `"LSI2"`, ..., `"LSIn"`.
#'
#' @details
#' When `v` is not supplied, it is calculated as
#' \deqn{V = X^T U D^{-1},}
#' where \eqn{X} is `mat`, \eqn{U} is `u`, and \eqn{D} is the diagonal matrix
#' containing the singular values `d`.
#'
#' The projected coordinates are then obtained from \eqn{VD}.
#'
#' @export
projectSVD <- function(mat, u, d, nDim, v = NULL) {
  if (is.null(v)) {
    v <- Matrix::t(mat) %*% u %*% diag(1 / d)
  }
  
  svdDiag <- matrix(0, nrow = nDim, ncol = nDim)
  diag(svdDiag) <- d
  
  matSVD <- as.matrix(Matrix::t(svdDiag %*% Matrix::t(v)))
  
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("LSI", seq_len(ncol(matSVD)))
  
  return(matSVD)
}
