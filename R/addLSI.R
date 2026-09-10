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
  message("Running SVD...")
  outliers = colnames(mat) %in% metadata(x)[[useMatrix]][["outliers"]]#returns NULL if doesn't exist in case of stacking
  if(is.null(outliers)){
    message("Either No outliers were detected or metadata(x)[[useMatrix]][['outliers'] does not exist. Confirm prior to exceeding")
    }
  mat = mat[rowSums(mat[,!outliers])>0,]
  svd <- irlba::irlba(mat[,!outliers], nDimensions, nDimensions)
  
  matSVD = projectSVD(mat,svd$u,svd$d,nDimensions)
  if (scaleDims) {
    # check and see if this is doing it right!
    message("Scaling matSVD...")
    matSVD <- rowZscores(matSVD)
  }
  message("Checking for depth-correlated columns...")
  toKeep <- which(cor(matSVD, colData(x)[[depth]])[, 1] < corCutOff)
  message("Kept ", (length(toKeep)/ncol(matSVD)*100), "% of columns.")
  matSVD <- matSVD[, toKeep,drop=FALSE]
  colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
  reducedDim(x, name) <- matSVD
  message("Done.")
  return(x) 
}