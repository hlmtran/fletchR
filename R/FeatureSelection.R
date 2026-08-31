# tfidfCalc = function(mat,scaleTo=10000){
#   stopifnot(inherits(mat, "sparseMatrix"))
#   # Compute term frequency (TF)
#   colSm <- Matrix::colSums(mat)
#   mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
#   
#   # Compute inverse document frequency (IDF)
#   rowSm <- Matrix::rowSums(mat)
#   idf <- as(ncol(mat) / rowSm, "sparseVector")
#   
#   # TF-IDF
#   mat <- as(Matrix::Diagonal(x = as.vector(idf)), "sparseMatrix") %*% mat
#   
#   # Log transform
#   mat@x <- log(mat@x * scaleTo + 1)
#   
#   return(mat)
# }
# 
# filterByQuantile = function(mat,quantile=0.995,varFeatures=25000,totalFeatures=500000){
#   stopifnot(inherits(mat, "Matrix"))
#   if (filterQuantile < 0 || filterQuantile > 1) {
#     stop("filterQuantile must be between 0 and 1.")
#   }
#   if (varFeatures < 1000) {
#     stop("varFeatures must be at least 1.")
#   }
#   # Total counts per feature
#   featureTotals <- Matrix::rowSums(mat)
#   # Rank features from highest to lowest total count
#   ranked <- order(featureTotals, decreasing = TRUE)
#   # Number of highest-count features to discard
#   removeTop <- floor((1 - filterQuantile) * totalFeatures)
#   if (sum(featureTotals > 0) > 2.25 * varFeatures) {
#     candidateCount <- min(
#       varFeatures + removeTop,
#       length(ranked)
#     )
#     selected <- head(ranked, candidateCount)
#     if (removeTop > 0) {
#       selected <- selected[-seq_len(min(removeTop, length(selected)))]
#     }
#     selected <- head(selected, varFeatures)
#   } else {
#     message("Not enough non-zero features to apply upper-tail filtering.")
#     selected <- head(ranked, varFeatures)
#   }
#   # Remove any zero-count features
#   selected <- selected[featureTotals[selected] > 0]
#   # Preserve original feature order
#   selected <- sort(selected)
#   mat[selected, , drop = FALSE]
# }
# 
# #from ArchR
# outlierByQuantile = function(mat,quantiles=c(0.02,0.98)){
#   stopifnot(inherits(mat, "Matrix"))
#   if (
#     length(quantiles) != 2 ||
#     any(quantiles < 0) ||
#     any(quantiles > 1)
#   ) {
#     stop("quantiles must contain two values between 0 and 1.")
#   }
# 
#   colSm = colSums(mat)
#   thresholds = stats::quantile(
#     colSm,
#     probs = sort(quantiles),
#     names = FALSE
#   )
# 
#   idxOutlier <- which(
#     colSm <= thresholds[1] |
#       colSm >= thresholds[2]
#   )
#   return(idxOutlier)
# }
# 
# remove0Rows = function(mat,rowSm=NULL){
#   if(is.null(rowSm)){
#     rowSm = rowSums(mat)
#   }
#   mat = mat[rowSm>0,]
#   return(mat)
# }
# 
# calcIDF = function(mat){
#   binMat = mat > 0
#   rSum = rowSums(binMat)
#   idf = as(ncol(mat)/(rSum), "sparseVector")
#   return(idf)
# }
# 
# testCorrelation <- function(pCheck, pCheck2) {
#   lapply(seq_len(ncol(pCheck)), function(x) {
#     cor(pCheck[, x], pCheck2[, x])
#   })
# }