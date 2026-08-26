#' Filter features and extract an assay matrix
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF")
#' @param excludeChr    chroms to exclude by default (chrM, chrX, chrY)
#' @param subsetLSI     subset to mcols(x)$usedForLSI? (FALSE)
#' @param replaceZeros  logical indicating whether filtered features should be replaced with 0 instead
#' @param prune         Minimum number of features in a cell must have and minumum of cells containing a feature
#' c(features,cells)
#' @return              filtered assay matrix
#'
#' @details             Subsets to features marked by `mcols(x)$usedForLSI`
#'                      (if requested), removes excluded chromosomes, and
#'                      returns the requested assay.
#'
#' @import              GenomeInfoDb
#' @import              SingleCellExperiment
#'
#' @keywords internal
filterAndGetMat <- function(
    x,
    useMatrix = c("TFIDF","counts","TileMatrix"),
    excludeChr = c("chrM", "chrX", "chrY"),
    subsetLSI = FALSE,
    prune = c(1,1),
    replaceZeros = TRUE
) {
  
  # keep <- setdiff(seqlevels(x), excludeChr)
  chrKeep = !getChrInd(x=x,chr=excludeChr)
  if (subsetLSI) {
    stopifnot("usedForLSI" %in% names(mcols(x)))
    idx <- mcols(x)$usedForLSI
  } else {
    idx <- rep(TRUE,nrow(x))
  }
  # useMatrix <- match.arg(useMatrix)
  message("Subsetting assay matrix ", useMatrix, "...")
  # mat <- assay(keepSeqlevels(x[idx, ], keep, pruning.mode = "coarse"),
  #              useMatrix)
  mat <- assay(x,useMatrix)
  
  rowToPrune = pruneRows(mat,prune=prune[1])
  colToPrune = pruneCols(mat,prune=prune[2])
  
  keepRows = idx&chrKeep&(!rowToPrune)
  keepCols = !colToPrune
  
  if(replaceZeros){
    origRownames = rownames(mat)
    origColnames = colnames(mat)
    
    mat = as(Matrix::Diagonal(x=keepRows),"sparseMatrix") %*% mat
    mat = mat %*% as(Matrix::Diagonal(x=keepCols),"sparseMatrix")
    
    mat = Matrix::drop0(mat)
    colnames(mat) = origColnames
    rownames(mat) = origRownames
  }else{
    mat <- mat[keepRows, keepCols]
  }
  return(mat)
}

#' get chromosome indices
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param chr           chroms to find
#' @return              return indices of target chromosomes
#'
#' @import              GenomeInfoDb
#' @import              SingleCellExperiment
#'
#' @keywords internal
getChrInd = function(x,chr = c("chrM", "chrX", "chrY")){
  ind = as.character(seqnames(x)) %in% chr
  return(ind) 
}