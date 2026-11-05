#' Select features based on across cell accessibility
#' 
#' @param mat           a sparse matrix (dgCMatrix) of features x cells
#' @param filterQuantile numeric value between 0 and 1; the quantile
#'                       of features to keep
#' @param varFeatures    integer; the number of features to return after filtering
#' @param totalFeatures  integer; the total number of features in the matrix
#' @return               character vector of feature names
#' 
findTopFeatures = function(mat,filterQuantile=0.995,varFeatures=25000,totalFeatures=500000,balance=F){
  stopifnot(inherits(mat, "Matrix"))
  if (filterQuantile < 0 || filterQuantile > 1) {
    stop("Quantile must be between 0 and 1.")
  }
  if (!balance && varFeatures < 1000) {
    stop("varFeatures must be at least 1000.")
  }
  # Total counts per feature
  featureTotals <- Matrix::rowSums(mat)
  
  if(balance){
    selected = selectBalancedFeatures(featureTotals = featureTotals,
                                      filterQuantile = filterQuantile,
                                      varFeatures = varFeatures,
                                      totalFeatures = totalFeatures)
  }else {
    selected = selectTopFeatures(featureTotals = featureTotals,
                                 filterQuantile = filterQuantile,
                                 varFeatures = varFeatures,
                                 totalFeatures = totalFeatures)
    
  }
  
  return(rownames(mat)[rownames(mat) %in% selected])
}

#' Select features based on across cell accessibility
#' 
#' @param mat           a sparse matrix (dgCMatrix) of features x cells
#' @param filterQuantile numeric value between 0 and 1; the quantile
#'                       of features to keep
#' @param varFeatures    integer; the number of features to return after filtering
#' @param totalFeatures  integer; the total number of features in the matrix
#' @return               character vector of feature names
#' 
selectBalancedFeatures = function(featureTotals,filterQuantile=0.995,varFeatures=25000,totalFeatures=500000){
  #this needs to go to minimum number of features available across all 3 modalities or at least warn that number of features are inherently unbalanced
  stopifnot(
    "Feature names must end with a valid modality suffix (_m<modality>)." =
      all(grepl("_m[^_]+$", names(featureTotals)))
  )
  modalities = getModality(names(featureTotals))
  
  nModality = length(unique(modalities))
  
  featureTotals = split(featureTotals,modalities)
  
  selectedList <- lapply(
    featureTotals,
    function(x) selectTopFeatures(
      featureTotals = x,
      filterQuantile = filterQuantile,
      varFeatures = floor(varFeatures / nModality),
      totalFeatures = length(x)
    )
  )
  return(unlist(selectedList,use.names = F))
}

selectTopFeatures = function(featureTotals,filterQuantile=0.995,varFeatures=25000,totalFeatures=500000){
  origRowNames = names(featureTotals)
  ranked <- order(featureTotals, decreasing = TRUE)
  # Number of highest-count features to discard
  removeTop <- floor((1 - filterQuantile) * totalFeatures)
  if (sum(featureTotals > 0) > 2.25 * varFeatures) {
    candidateCount <- min(
      varFeatures + removeTop,
      length(ranked)
    )
    selected <- head(ranked, candidateCount)
    if (removeTop > 0) {
      selected <- selected[-seq_len(min(removeTop, length(selected)))]
    }
    selected <- head(selected, varFeatures)
  } else {
    message("Not enough non-zero features to apply upper-tail filtering.")
    selected <- head(ranked, varFeatures)
  }
  # # Remove any zero-count features
  selected <- selected[featureTotals[selected] > 0]
  # # Preserve original feature order
  # selected <- sort(selected)
  # mat[selected, , drop = FALSE]
  return(origRowNames[selected])
}