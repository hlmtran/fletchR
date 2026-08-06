#' Filter features and extract an assay matrix
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param useMatrix     name of assay to use ("TFIDF")
#' @param excludeChr    chroms to exclude by default (chrM, chrX, chrY)
#' @param subsetLSI     subset to mcols(x)$usedForLSI? (TRUE)
#'
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
    subsetLSI = TRUE
) {
  keep <- setdiff(seqlevels(x), excludeChr)
  if (subsetLSI) {
    stopifnot("usedForLSI" %in% names(mcols(x)))
    idx <- which(mcols(x)$usedForLSI)
  } else {
    idx <- seq_len(nrow(x))
  }
  message("Subsetting assay matrix...")
  mat = assay(keepSeqlevels(x[idx, ],keep,pruning.mode = "coarse"),useMatrix)
  return(mat)
}
