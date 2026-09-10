findTopFeatures = function(mat,quantile=0.995,varFeatures=25000,totalFeatures=500000){
  stopifnot(inherits(mat, "Matrix"))
  if (quantile < 0 || quantile > 1) {
    stop("Quantile must be between 0 and 1.")
  }
  if (varFeatures < 1000) {
    stop("varFeatures must be at least 1000.")
  }
  # Total counts per feature
  featureTotals <- Matrix::rowSums(mat)
  # Rank features from highest to lowest total count
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
  # Remove any zero-count features
  selected <- selected[featureTotals[selected] > 0]
  # Preserve original feature order
  selected <- sort(selected)
  mat[selected, , drop = FALSE]
}