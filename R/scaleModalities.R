#' Scale modalities in a sparse matrix
#'
#' Scales the stored values of each modality independently while preserving
#' the sparse matrix structure.
#'
#' @param mat A dgCMatrix with features in rows and cells in columns.
#' @param modalities Character vector assigning a modality to each row.
#'   Modalities must occur in consecutive blocks.
#' @param method Scaling method: "Frobenius" or "minMax".
#'
#' @return A scaled dgCMatrix.
scaleModalities = function(mat,modalities,method = "Frobenius"){
  stopifnot(inherits(mat,"dgCMatrix"))
  rowCoordinates = getModalityRowStartEnd(modalities=modalities)
  
  for(modality in unique(modalities)){
    
    rowInd = getSparseIndices(
      rowRanges=rowCoordinates,
      sparseRowIndices=mat@i,
      modality = modality
      )
    mat@x[rowInd] = scaleValues(v=mat@x[rowInd],method = method)
    
  }
  return(mat)
}

#' Extract modality labels from row names
#'
#' Extracts modality labels following "_m" at the end of each row name.
#'
#' @param v Character vector of row names containing "_m" followed by a modality.
#'
#' @return Character vector of modality labels.
getModality = function(v) {
  stopifnot(all(grepl("_m[^_]+$", v)))
  sub(".*_m([^_]+)$", "\\1", v)
}

#' Get modality row boundaries
#'
#' Computes zero-based, row boundaries for consecutive modality blocks.
#'
#' @param modalities Character vector of modality labels.
#'
#' @return Data frame with start and end columns, indexed by modality.
getModalityRowStartEnd = function(modalities){
  
  counts = rle(modalities)
  stopifnot(
    "Modalities must occur in consecutive blocks." =
      anyDuplicated(counts$values) == 0L
  )
  
  counts = setNames(counts$lengths,counts$values)
  ends = cumsum(counts)
  starts = ends - counts
  ends = ends - 1
  df = data.frame(start = starts,end = ends,row.names = names(counts))
  
  return(df)
}

#' Identify sparse indices belonging to a modality
#'
#' @param sparseRowIndices Zero-based row indices from `mat@i` of a `dgCMatrix`.
#' @param rowRanges Data frame returned by `getModalityRowStartEnd()`.
#' @param modality Name of the modality to select.
#'
#' @return Logical vector for subsetting `mat@x`, indicating which stored
#'   values belong to the requested modality.
getSparseIndices = function(rowRanges,sparseRowIndices,modality){
  
  stopifnot(modality %in% rownames(rowRanges))
  sparseIdx = sparseRowIndices <= rowRanges[modality,"end"] & sparseRowIndices >= rowRanges[modality,"start"]
  return(sparseIdx)
}

#' Scale numeric values
#'
#' Applies Frobenius normalization or min-max scaling.
#'
#' @param v Numeric vector of values.
#' @param method Scaling method: "Frobenius" or "minMax".
#'
#' @return Numeric vector of scaled values.
scaleValues = function(v, method = "Frobenius") {
  stopifnot(!all(v==0))
  if (method == "Frobenius") {
    v = v / sqrt(sum(v^2))
    
  } else if (method == "minMax") {
    v = (v - min(v)) / (max(v) - min(v))
    
  } else {
    stop("Invalid method")
  }
  
  return(v)
}