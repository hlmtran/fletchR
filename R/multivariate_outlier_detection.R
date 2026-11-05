#' Detect multivariate outliers
#'
#' Identifies outliers using transformed variables and Mahalanobis distance.
#'
#' @param sce SingleCellExperiment object.
#' @param cols Columns used for outlier detection.
#' @param FN Column transformation function.
#' @param quantile_cutoff Chi-squared quantile threshold.
#' @param outlierColname colData name to save.
#'
#' @return SingleCellExperiment object with an `outliers` column.
#'
#' @export
detectOutliers = function(sce,
                          cols = c("nFrags_DEM", "nFrags_K4"),
                          FN = normit,
                          quantile_cutoff = 0.95,
                          outlierColname = "outliers") {
  
  variables = getVariables(sce, cols)
  
  if(!is.null(FN)){  variables = apply(variables, 2, FN)}
  
  outliers = detectMahalanobisOutliers(
    mat = variables,
    quantile_cutoff = quantile_cutoff
  )
  
  stopifnot(identical(rownames(variables),names(outliers)))
  sce[[outlierColname]] = outliers
  
  return(sce)
}


#' Extract and transform variables
#'
#' Extracts selected colData columns and transforms each independently.
#'
#' @param sce SingleCellExperiment object.
#' @param cols Columns to extract.
#'
#' @return Matrix of transformed variables.
#'
#' @importFrom SummarizedExperiment colData
getVariables = function(sce, cols = c("nFrags_DEM", "nFrags_K4")) {
  stopifnot(
    is.character(cols),
    length(cols) >= 1,!anyNA(cols),!anyDuplicated(cols),
    all(cols %in% colnames(colData(sce)))
  )
  if (length(cols) == 1) {
    warning("Only 1 variable detected. This test is not robust to low end outliers in 1 dimension")
  }
  variables = colData(sce)
  variables = variables[, cols, drop = FALSE]
  
  stopifnot(all(vapply(variables, is.numeric, logical(1))))
  
  return(variables)
}


#' Rank-based inverse normal transformation
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

plotOutliers <- function(
    sceObj,
    dataCols = c("nFrags_K4", "nFrags_DEM"),
    outlierColname = "outliers",
    FN = normit,
    ...
) {
  
  stopifnot(length(dataCols) == 2L)
  stopifnot(outlierColname %in% colnames(colData(sceObj)))
  
  # Extract and transform data
  x <- FN(sceObj[[dataCols[1]]])
  y <- FN(sceObj[[dataCols[2]]])
  
  # Shared axis limits based on all cells
  xlim <- range(x, finite = TRUE)
  ylim <- range(y, finite = TRUE)
  
  # Identify valid cells
  idx <- is.finite(x) & is.finite(y)
  
  # Identify non-outlier cells
  idxNoOutliers <- idx &
    !is.na(sceObj[[outlierColname]]) &
    !sceObj[[outlierColname]]
  
  # Enable plot recording
  grDevices::dev.control(displaylist = "enable")
  
  # Plot all cells
  graphics::smoothScatter(
    x[idx], y[idx],
    xlim = xlim,
    ylim = ylim,
    xlab = dataCols[1],
    ylab = dataCols[2],
    main = "All cells",
    ...
  )
  
  pAll <- grDevices::recordPlot()
  
  # Plot without outliers
  graphics::smoothScatter(
    x[idxNoOutliers], y[idxNoOutliers],
    xlim = xlim,
    ylim = ylim,
    xlab = dataCols[1],
    ylab = dataCols[2],
    main = "Without outliers",
    nrpoints = 0,
    ...
  )
  
  pFiltered <- grDevices::recordPlot()
  
  return(list(
    all = pAll,
    filtered = pFiltered
  ))
}
