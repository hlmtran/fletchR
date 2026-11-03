#' Detect multivariate outliers
#'
#' Identifies outliers using transformed variables and Mahalanobis distance.
#'
#' @param sce SingleCellExperiment object.
#' @param cols Columns used for outlier detection.
#' @param FN Column transformation function.
#' @param quantile_cutoff Chi-squared quantile threshold.
#'
#' @return SingleCellExperiment object with an `outliers` column.
#'
#' @export
detectOutliers = function(sce,
                          cols = c("nFrags_DEM", "nFrags_K4"),
                          FN = normit,
                          quantile_cutoff = 0.95) {
  
  variables = getVariables(sce, cols, FN)
  
  outliers = detectMahalanobisOutliers(
    mat = variables,
    quantile_cutoff = quantile_cutoff
  )
  
  sce$outliers = outliers
  
  return(sce)
}


#' Extract and transform variables
#'
#' Extracts selected colData columns and transforms each independently.
#'
#' @param sce SingleCellExperiment object.
#' @param cols Columns to extract.
#' @param FN Column transformation function.
#'
#' @return Matrix of transformed variables.
#'
#' @importFrom SummarizedExperiment colData
getVariables = function(sce,
                        cols = c("nFrags_DEM", "nFrags_K4"),
                        FN = normit) {
  
  stopifnot(
    is.character(cols),
    length(cols) >= 1,
    !anyNA(cols),
    !anyDuplicated(cols),
    all(cols %in% colnames(colData(sce))),
    is.function(FN)
  )
  
  nfrags_mat = colData(sce)
  nfrags_mat = nfrags_mat[, cols, drop = FALSE]
  
  stopifnot(all(vapply(nfrags_mat, is.numeric, logical(1))))
  
  nfrags_mat = apply(nfrags_mat, 2, FN)
  
  return(nfrags_mat)
}


#' Continuity-corrected inverse normal transformation
#'
#' Applies a rank-based inverse normal transformation.
#'
#' @param x Numeric vector to transform.
#'
#' @return Transformed numeric vector.
#'
#' @importFrom stats qnorm
normit = function(x) {
  
  notNA = which(!is.na(x))
  n = length(x[notNA])
  
  if (n == 0) {
    return(x)
  }
  
  x[notNA] = qnorm((rank(x[notNA]) - 0.5) / n)
  
  return(x)
}


#' Detect Mahalanobis outliers
#'
#' Flags observations exceeding a chi-squared Mahalanobis distance threshold.
#'
#' @param mat Numeric matrix with observations in rows.
#' @param quantile_cutoff Chi-squared quantile threshold.
#'
#' @return Logical vector indicating outliers.
#'
#' @importFrom stats mahalanobis cov qchisq
detectMahalanobisOutliers = function(mat, quantile_cutoff = 0.95) {
  
  stopifnot(
    is.matrix(mat),
    is.numeric(mat),
    !anyNA(mat),
    all(is.finite(mat)),
    length(quantile_cutoff) == 1,
    is.numeric(quantile_cutoff),
    is.finite(quantile_cutoff),
    quantile_cutoff > 0,
    quantile_cutoff < 1,
    nrow(mat) > ncol(mat) #invertible cov
  )
  
  D2 = mahalanobis(mat, colMeans(mat), cov(mat))
  
  cutoff = qchisq(quantile_cutoff, df = ncol(mat))
  
  outliers = D2 > cutoff
  
  return(outliers)
}