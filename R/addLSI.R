#' add an LSI result to an SCE (NOTE: this assumes prefiltering was done!) 
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF") 
#' @param name          name of reducedDim to store the LSI embedding ("LSI") 
#' @param scaleDims     Z-scale the dimensions? (TRUE)
#' @param corCutOff     drop dimensions where cor(LSIdim, nFrag) > corCutOff
#' @param excludeChr    chroms to exclude by default (chrM, chrX, chrY) 
#' @param features      Features to keep for LSI
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
addLSI <- function(x, metadataSlot="TFIDF", name="LSI", scaleDims=FALSE, corCutOff=0.75, excludeChr=c("chrM","chrX","chrY"), features=NULL ,nDimensions=30, depth="nFrags", seed = 1,...) { 
  set.seed(seed)
  stopifnot(depth %in% names(colData(x)))
  # useMatrix <- match.arg(useMatrix) 
  # keep <- setdiff(seqlevels(x), excludeChr)
  # if (subsetLSI) {
  #   stopifnot("usedForLSI" %in% names(mcols(x)))
  #   idx <- which(mcols(x)$usedForLSI)
  # } else {
  #   idx <- seq_len(nrow(x))
  # }
  # message("Subsetting TF-IDF matrix...") 
  # mat <- assay(keepSeqlevels(x[idx,], keep, pruning.mode="coarse"), useMatrix)
  # mat <- filterAndGetMat(x=x,useMatrix=useMatrix,excludeChr=excludeChr,subsetLSI=subsetLSI)
  mat = metadata(x)[[metadataSlot]]
  if(!is.null(features)){
    mat = mat[rownames(mat) %in% features,]
    }
  message("Running SVD...")
  outliers = colnames(mat) %in% attr(mat, 'outliers')
  svd <- irlba::irlba(mat[,!outliers], nDimensions, nDimensions)
  # svdDiag <- matrix(0, nrow=nDimensions + 5, ncol=nDimensions + 5)
  # diag(svdDiag) <- svd$d
  # matSVD <- t(svdDiag %*% t(svd$v))
  # rownames(matSVD) <- colnames(mat)
  matSVD = projectSVD(mat,svd$u,svd$d,nDimensions)
  if (scaleDims) {
    # check and see if this is doing it right!
    message("Scaling matSVD...")
    matSVD <- rowZscores(matSVD)
  }
  message("Checking for depth-correlated columns...")
  toKeep <- which(cor(matSVD, colData(x)[[depth]])[, 1] < corCutOff)
  message("Kept ", (length(toKeep)/ncol(matSVD)*100), "% of columns.")
  matSVD <- matSVD[, toKeep][, seq_len(nDimensions)]
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  reducedDim(x, name) <- matSVD
  message("Done.")
  return(x) 
}

projectSVD = function(mat,u,d,nDim){
  V <- Matrix::t(mat) %*% u %*% diag(1/d)
  svdDiag <- matrix(0,nrow=nDim,ncol=nDim)# D in UDt(V)
  diag(svdDiag) <- d                                    # as above
  matSVD <- as.matrix(Matrix::t(svdDiag %*% Matrix::t(V)))      # project mat
  rownames(matSVD) <- colnames(mat)                             # cell names
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))       # LSIdim names
  return(matSVD)     
}
