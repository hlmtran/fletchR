#' Get ArchR LSI features
#'
#' Extracts the feature names used for latent semantic indexing (LSI) from a
#' SingleCellExperiment created from an ArchR object.
#'
#' @param x a SingleCellExperiment, usually from archRtoSCE
#'
#' @return A character vector containing the names of features used for LSI.
#'
#' @details Features are selected based on the logical `usedForLSI` column
#'   in `rowData(x)`. An error is raised if `usedForLSI` is not present in
#'   the row metadata.
#'
#' @import SingleCellExperiment
#'
#' @export
#'
getArchRLSIFeatures <- function(x) {
  stopifnot("usedForLSI" %in% names(mcols(x)))
  features <- rownames(rowData(x)[rowData(x)$usedForLSI, ])
  return(features)
}

