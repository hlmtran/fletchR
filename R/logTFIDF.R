#' Compute log-transformed TF-IDF
#'
#' Multiplies a term-frequency matrix by an inverse document frequency vector
#' and applies the ArchR log transformation
#' \eqn{\log(\mathrm{TFIDF} \times scaleTo + 1)}.
#'
#' @param tf A term-frequency normalized sparse matrix.
#' @param idf A sparse vector of inverse document frequencies.
#' @param scaleTo Scaling factor applied prior to the log transformation (1e4)
#'
#' @return
#' A sparse log-transformed TF-IDF matrix.
#'
#' @import Matrix 
#' 
#' @export
logTFIDF <- function(tf, idf, scaleTo=10000){

  origRows <- rownames(tf)
  log_tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  log_tfidf@x <- log(log_tfidf@x * scaleTo + 1) # iffy but consistent with ArchR
  rownames(log_tfidf) <- origRows
  return(log_tfidf)

}
