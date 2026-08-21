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
#' @param prune Minimum number of non-outlier cells in which a feature must be
#'   observed to be retained.
#' @param metadataSlot Character scalar giving the name under which the
#'   transformed matrix is stored in \code{metadata(x)}.
#' @param keep0col Logical indicating to keep columns with colSum = 0
#' @param ... Additional arguments passed to \code{logTFIDF()}.
#'
#' @return The input \code{SummarizedExperiment} with a log-transformed TF-IDF
#'   matrix stored in \code{metadata(x)[[metadataName]]}. The stored matrix has
#'   attributes \code{"idf"} containing the computed inverse document
#'   frequencies and \code{"outliers"} containing the names of cells excluded
#'   from IDF estimation.
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
addTfIdf <- function(x, useMatrix=c("counts","TileMatrix"), excludeChr=c("chrM","chrX","chrY"), subsetLSI=FALSE, binarize=TRUE, outlierQuantiles=c(0.02, 0.98), prune=1, metadataSlot="TfIdf", keep0col=FALSE,  ...) { 

  if (is(x, "SummarizedExperiment")) {
    # useMatrix <- match.arg(useMatrix)
    mat <- filterAndGetMat(x=x, useMatrix=useMatrix, excludeChr=excludeChr, 
                           subsetLSI=subsetLSI)
  }

  # if (is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  
  if (binarize) mat <- binarizeMat(mat)
  
  #remove 0 accessibility cells
  if(!keep0col){
      mat = mat[,!pruneCols(mat,prune=1)]
    }
  
  if (!is.null(outlierQuantiles)){
    idxOutliers <- outlierByQuantile(mat,outlierQuantiles)
  } else {
    idxOutliers <- logical(ncol(mat))
  }
  
  outliers <- colnames(mat)[idxOutliers]
  
  idfMat <- mat[, !idxOutliers, drop = FALSE]
  keep <- !pruneRows(idfMat, prune)
  
  mat <- getTF(mat[keep,],keep0col)
  idfMat <- getIDF(idfMat[keep,])
  message("Adding ", metadataSlot, " to metadata(x)$", metadataSlot, "...")
  metadata(x)[[metadataSlot]] <- logTFIDF(tf=mat,idf=idfMat , ...)
  attr(metadata(x)[[metadataSlot]], 'idf') <- idfMat
  attr(metadata(x)[[metadataSlot]], 'outliers') <- outliers
  return(x)
}
