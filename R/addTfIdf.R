#' Add an ArchR-style log1p(TF-IDF) matrix to metadata
#'
#' Computes an ArchR-style log-transformed term frequency-inverse document
#' frequency (TF-IDF) matrix from a sparse count matrix stored in a
#' \code{SummarizedExperiment}. The transformed matrix is stored in
#' \code{metadata(x)} together with attributes containing the computed inverse
#' document frequencies (IDF) and the names of outlier cells.
#'
#' @param x A \code{SummarizedExperiment}.
#' @param useMatrix name of assay to use ("counts") 
#' @param excludeChr chroms to exclude by default (chrM, chrX, chrY) 
#' @param subsetLSI Logical indicating whether only features selected for latent
#'   semantic indexing (LSI) should be used.
#' @param binarize Binarize matrix
#' @param outlierQuantiles Numeric vector of length two specifying the lower and
#'   upper quantiles used to identify outlier cells by total accessibility. Set
#'   to \code{NULL} to disable outlier detection.
#' @param assayName name of assay
#' @param ... Additional arguments passed to \code{logTFIDF()}.
#'
#' @return TfIdf matrix in assay
#'
#' @details Outlier cells are excluded when estimating inverse document
#'   frequencies but remain in the transformed matrix. Features observed in
#'   fewer than \code{prune} non-outlier cells are omitted from the stored
#'   TF-IDF matrix. The input \code{SummarizedExperiment} is not subsetted or
#'   otherwise modified apart from the addition of the metadata entry.
#'
#' @seealso \code{\link{logTFIDF}}, \code{\link{TF}}, \code{\link{IDF}}
#'
#' @import SummarizedExperiment
#' @import GenomeInfoDb
#'
#' @export
addTfIdf <- function(x, useMatrix=c("counts","TileMatrix"), excludeChr=c("chrM","chrX","chrY"), subsetLSI=FALSE, binarize=TRUE, outlierQuantiles=c(0.02, 0.98), assayName="TfIdf",...) { 

  if (is(x, "SummarizedExperiment")) {
    # useMatrix <- match.arg(useMatrix)
    mat <- filterAndGetMat(x=x, useMatrix=useMatrix, excludeChr=excludeChr, 
                           subsetLSI=subsetLSI)
  }

  # if (is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  
  if (binarize) mat <- binarizeMat(mat)

  if (!is.null(outlierQuantiles)){
    idxOutliers <- outlierByQuantile(mat,outlierQuantiles,excludeZeros=TRUE) #contains both outliers and 0 columns
    idx0ColSum = colSums(mat) == 0
    idxOutliers = idxOutliers & !idx0ColSum # keep only outliers
  } else {
    idxOutliers <- logical(ncol(mat))
  }

  
  outliers <- colnames(mat)[idxOutliers]
  
  idfMat <- mat[, !idxOutliers, drop = FALSE]
  
  mat <- getTF(mat)
  idfMat <- getIDF(idfMat)
  
  message("Adding ", assayName, " to assays...")
  assay(x,assayName) <- logTFIDF(tf=mat,idf=idfMat , ...)
  metadata(x)[[assayName]] <- list(
    idf = idfMat,
   outliers = outliers 
  )
  # x@metadata[[assayName]][['idf']] = idfMat
  # x@metadata[[assayName]][['outliers']] = outliers
  # attr(assay(x,assayName), 'idf') <- idfMat
  # attr(assay(x,assayName), 'outliers') <- outliers
  return(x)
}
