#' Add an ArchR-style log1p(TF-IDF) matrix to metadata
#'
#' Computes an ArchR-style log-transformed term frequency-inverse document
#' frequency (TF-IDF) matrix from a sparse count matrix stored in a
#' \code{SummarizedExperiment}. The transformed matrix is stored in
#' \code{metadata(x)} together with attributes containing the computed inverse
#' document frequencies (IDF) and the names of outlier cells.
#'
#' @param x A \code{SummarizedExperiment}.
#' @param useMatrix name of assay to use ("TFIDF") 
#' @param excludeChr chroms to exclude by default (chrM, chrX, chrY) 
#' @param subsetLSI Logical indicating whether only features selected for latent
#'   semantic indexing (LSI) should be used.
#' @param binarize Binarize matrix
#' @param outlierQuantiles Numeric vector of length two specifying the lower and
#'   upper quantiles used to identify outlier cells by total accessibility. Set
#'   to \code{NULL} to disable outlier detection.
#' @param prune Minimum number of non-outlier cells in which a feature must be
#'   observed to be retained.
#' @param metadataName Character scalar giving the name under which the
#'   transformed matrix is stored in \code{metadata(x)}.
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
#'
#' @export
addTfIdf <- function(x,useMatrix=c("counts","TileMatrix"),excludeChr=c("chrM","chrX","chrY"),subsetLSI = TRUE,binarize=TRUE,outlierQuantiles=c(0.02,0.98),prune=1,metadataSlot="TFIDF",...) { 
  if (is(x, "SummarizedExperiment")) {
    mat <- filterAndGetMat(x=x,useMatrix=useMatrix,excludeChr=excludeChr,subsetLSI)
  }
  # if (is(idf, "sparseMatrix")) idf <- attr(idf, 'idf') 
  if (!is(mat, "sparseMatrix")) stop("logTfIdf only works on sparse matrices") 
  if (binarize) mat <- binarizeMat(mat)
  # if (subsetLSI){
  #   stopifnot("usedForLSI" %in% names(mcols(x)))
  #   idx <- which(mcols(x)$usedForLSI)
  # }
  if(!is.null(outlierQuantiles)){
  idxOutliers = outlierByQuantile(mat,outlierQuantiles)
  }else{
    idxOutliers = logical(ncol(mat))
  }
  
  outliers <- colnames(mat)[idxOutliers]
  
  idfMat <- mat[, !idxOutliers, drop = FALSE]
  keep <- !pruneRows(idfMat, prune)
  
  mat    = TF(mat[keep,])
  idfMat = IDF(idfMat[keep,])
  
  metadata(x)[[metadataSlot]] <- logTFIDF(tf=mat,idf=idfMat , ...)
  
  attr(metadata(x)[[metadataSlot]], 'idf') <- idfMat
  attr(metadata(x)[[metadataSlot]], 'outliers') <- outliers
  message("Adding idf to metadata(x)$idf...")
  # metadata(x)$idf <- attr(assay(x, metadataName), 'idf')
  return(x)
}

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
#' @keywords internal
TF = function(mat){
  # directly borrowed from ArchR
  colSm <- Matrix::colSums(mat)
  # rowSm <- Matrix::rowSums(mat)
  message("Term frequency normalization... ", appendLF=FALSE)
  mat@x <- (mat@x / rep.int(colSm, Matrix::diff(mat@p)))
  message("done.")
  return(mat)
}

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
#' @return
#' A sparse vector containing one IDF value per feature.
#'
#' @keywords internal
IDF = function(mat){
  message("Computing IDF (inverse document frequency) table... ", appendLF=0)
  # idf <- as(((ncol(mat)+1)/(rowSm2+1)) + 1, "sparseVector")
  idf = as(ncol(mat)/rowSums(binarizeMat(mat)), "sparseVector")
  attr(idf, 'names') <- rownames(mat)
  message("done.")
  return(idf)
}

#' Compute log-transformed TF-IDF
#'
#' Multiplies a term-frequency matrix by an inverse document frequency vector
#' and applies the ArchR log transformation
#' \eqn{\log(\mathrm{TFIDF} \times scaleTo + 1)}.
#'
#' @param tf A term-frequency normalized sparse matrix.
#' @param idf A sparse vector of inverse document frequencies.
#' @param scaleTo Scaling factor applied prior to the log transformation. Default same as ArchR.
#'
#' @return
#' A sparse log-transformed TF-IDF matrix.
#'
#' @keywords internal
logTFIDF = function(tf,idf,scaleTo=10000){
  origRows = rownames(tf)
  log_tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% tf
  log_tfidf@x = log(log_tfidf@x * scaleTo + 1) #seems iffy but consistent with ArchR
  rownames(log_tfidf) <- origRows
  return(log_tfidf)
}

#' Identify low-frequency features
#'
#' Determines which features occur in fewer than a specified number of cells.
#'
#' @param mat A sparse count matrix.
#' @param prune Minimum number of cells in which a feature must be observed.
#'
#' @return
#' A logical vector indicating which rows should be removed.
#'
#' @keywords internal
pruneRows = function(mat,prune=1){
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

#' Identify outlier cells by total accessibility
#'
#' Identifies cells whose total counts fall outside user-specified quantiles.
#'
#' @param mat A sparse count matrix.
#' @param quantiles Numeric vector of length two giving lower and upper
#'   quantiles used to define outliers.
#'
#' @return
#' logical vector of outlier columns.
#'
#' @keywords internal
outlierByQuantile = function(mat,quantiles=c(0.02,0.98)){
  #from ArchR
  stopifnot(inherits(mat, "Matrix"))
  if (
    length(quantiles) != 2 ||
    any(quantiles < 0) ||
    any(quantiles > 1)
  ) {
    stop("quantiles must contain two values between 0 and 1.")
  }
  colSm = colSums(mat)
  thresholds = stats::quantile(
    colSm,
    probs = sort(quantiles),
    names = FALSE
  )
  
  idxOutlier <- colSm <= thresholds[1] | colSm >= thresholds[2]

  return(idxOutlier)
}