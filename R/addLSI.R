#' add an LSI result to an SCE (NOTE: this assumes prefiltering was done!) 
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF") 
#' @param name          name of reducedDim to store the LSI embedding ("LSI") 
#' @param scaleDims     Z-scale the dimensions? (TRUE)
#' @param corCutOff     drop dimensions where cor(LSIdim, nFrag) > corCutOff
#' @param excludeChr    Chromosomes to exclude
#' @param nDimensions   number of dimensions for the SVD (30) 
#' @param depth         colData column holding depth for correlations ('nFrags')
#' @param subsetLSI     subset to mcols(x)$usedForLSI? (TRUE) 
#'
#' @return              SingleCellExperiment with reducedDim(x, name)
#'
#' @details             This function assumes that ArchR selected the features!
#'                      Here we do a simple, noniterative one-pass SVD on TFIDF
#'
#' @seealso             addTfIdf
#' @seealso             ArchR::addIterativeLSI
#'
#' @import              SingleCellExperiment
#' @import              GenomeInfoDb
#' @import              irlba
#'
#' @export
#'
addLSI <- function(x,
  useMatrix="TfIdf",
  name="LSI",
  scaleDims=FALSE,
  corCutOff=0.75,
  excludeChr=c("chrM","chrX","chrY"),
  nDimensions=30,
  depth="nFrags",
  seed = 1,
  subsetLSI = FALSE,...) {
  set.seed(seed)
  stopifnot(depth %in% names(colData(x)))
 
  mat <- filterAndGetMat(x=x,useMatrix=useMatrix,excludeChr=excludeChr,subsetLSI=subsetLSI,prune=c(1,1),replaceZeros = FALSE)
  # message("Running SVD...")

  outliers = getOutliersIdx(x,useMatrix=useMatrix)
  matSVD <- calcLSI(mat=mat, outliers=outliers, nDimensions=nDimensions, scaleDims=scaleDims)
  # mat = mat[rowSums(mat[,!outliers])>0,]
  # svd <- irlba::irlba(mat[,!outliers], nDimensions, nDimensions)
  
  # matSVD = projectSVD(mat,svd$u,svd$d,nDimensions)
  # if (scaleDims) {
  #   # check and see if this is doing it right!
  #   message("Scaling matSVD...")
  #   matSVD <- rowZscores(matSVD)
  # }

  matSVD <- removedepthCorrelatedDims(projectedMat=matSVD, depth=colData(x)[[depth]], corCutOff=corCutOff)

  reducedDim(x, name) <- matSVD
  message("Done.")
  return(x) 
}

#' Get the outlier cells from a SingleCellExperiment object
#'
#' @param x A SingleCellExperiment object.
#' @param useMatrix The name of the assay to use for outlier detection (default is "TfIdf").
#' @return A logical vector indicating which cells are outliers.
getOutliersIdx <- function(x, useMatrix="TfIdf") {
  outliers = metadata(x)[[useMatrix]][["outliers"]] # returns NULL if doesn't exist
  if (is.null(outliers)) {
    message("Either No outliers were detected or metadata(x)[[useMatrix]][['outliers'] does not exist. Confirm prior to exceeding")
  }
  return(colnames(x) %in% outliers)
}

#' Core computational engine: Compute LSI on a matrix
#'
#' @param mat Feature-by-cell matrix (e.g. TF-IDF).
#' @param nDimensions Number of singular values to calculate.
#' @param outliers Logical vector indicating outlier cells to ignore during SVD fitting.
#' @param scaleDims Logical; whether to Z-score LSI dimensions across cells.
#' @return A cell-by-dimension projected matrix.
calcLSI <- function(mat, 
                    nDimensions = 30, 
                    outliers = NULL, 
                    scaleDims = FALSE) {
  if (is.null(outliers)) {
    outliers <- rep(FALSE, ncol(mat))
  }
  message("Running SVD...")
  mat = mat[rowSums(mat[,!outliers])>0,]
  svd <- irlba::irlba(mat[,!outliers], nDimensions, nDimensions)
    
  matSVD = projectSVD(mat,svd$u,svd$d,nDimensions)
  if (scaleDims) {
    message("Scaling matSVD...")
    matSVD <- rowZscores(matSVD)
  } 
  
  return(matSVD)
}

#' Remove LSI dimensions that are correlated with sequencing depth
#'
#' @param projectedMat The cell-by-dimension matrix obtained from LSI.
#' @param depth A numeric vector representing sequencing depth for each cell.
#' @param corCutOff A numeric threshold for correlation; dimensions with correlation above this value will be removed.
#' @return A filtered cell-by-dimension matrix with depth-correlated dimensions removed.
removedepthCorrelatedDims <- function(projectedMat, depth, corCutOff=0.75) {
  message("Checking for depth-correlated columns...")
  toKeep <- which(cor(projectedMat, depth)[, 1] < corCutOff)
  message("Kept ", (length(toKeep)/ncol(projectedMat)*100), "% of columns.")
  projectedMat <- projectedMat[, toKeep, drop = FALSE]
  colnames(projectedMat) <- paste0("LSI",seq_len(ncol(projectedMat)))
  return(projectedMat)
}
