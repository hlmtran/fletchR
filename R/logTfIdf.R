#' ArchR-style log1p(TF-IDF) transformation but without mandatory binarization
#' 
#' @param mat       a (sparse) matrix of counts (cells as columns) or an SE
#' @param prune     prune (cells >= prune)? (1, by default)
#' @param idf       a pre-trained idf or logTfIdf result (NULL; compute idf)
#' @param binarize  binarize the values for TF-IDF? (TRUE)
#' 
#' @return          a log1p(TF-IDF) version of the (counts) matrix (see Details)
#'
#' @details         TF == (term / sum(terms_in_document));
#'                  IDF == (1 / sum(documents_with_term_not_zero));
#'                  currently only assays(mat)[1] gets transformed for SE/SCEs;
#'                  provided or generated model idf placed in attr(mat, 'idf')
#'
#' @seealso         text2vec::TfIdf
#'
#' @import          Matrix
#'
#' @export
#'
logTfIdf <- function(mat, prune=1, idf=NULL, binarize=TRUE) { 

  if (is(mat, "SummarizedExperiment")) {
    return(logTfIdf(assay(mat), prune=prune, idf=idf, binarize=binarize))
  }
  if (is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  if (binarize) mat <- binarizeMat(mat)
  origRows <- rownames(mat)

  # directly borrowed from ArchR
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  message("Term frequency normalization... ", appendLF=FALSE)
  mat@x <- (mat@x / rep.int(colSm, Matrix::diff(mat@p)))
  message("done.")

  if (!is.null(idf)) {
    message("Using precomputed IDF (inverse document frequency) table...")
    idfnames <- attr(idf, 'names')
    if (!all(idfnames == rownames(mat))) {
      stop("rownames(mat) != attr(idf, 'names'). Subset your matrix first.")
    }
  } else { 
    rowSm2 <- rowSums(mat > 0)  # doesn't binarization already do this?!
    toPrune <- (rowSm2 < prune)
    if (sum(toPrune) > 0) {
      message("Minimum observed feature frequency: ", min(rowSm2[rowSm2 > 0]))
      message("Minimum allowable feature frequency: ", prune)
      message("Removing ",sum(toPrune)," of ",length(toPrune)," features.")
      # keep <- as(Matrix::Diagonal(x=!toPrune), "sparseMatrix")
      # rownames(keep) <- rownames(mat) 
      # mat <- keep %*% mat
      mat = mat[!toPrune,]
      origRows = origRows[!toPrune]
      rowSm2 <- rowSums(mat > 0)
      message("Minimum feature frequency is ", min(rowSm2[rowSm2 > 0]),
              " (", sum(!toPrune), "/", length(toPrune), " features retained).")
    }
    message("Computing IDF (inverse document frequency) table... ", appendLF=0)
    # idf <- as(((ncol(mat)+1)/(rowSm2+1)) + 1, "sparseVector")
    idf = as(ncol(mat)/(rowSm2), "sparseVector")
    attr(idf, 'names') <- rownames(mat)
    message("done.")
  }
  mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
  message("Computing log1p(TF-IDF)... ", appendLF=FALSE)
  mat@x <- log(mat@x * 10001)
  rownames(mat) <- origRows
  attr(mat, 'idf') <- idf
  message("done.")
  return(mat)

}
