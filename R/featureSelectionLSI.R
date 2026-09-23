featureSelectionLSI = function(x,
  useMatrix="counts",
  preTfIdf=NULL,
  binarize=TRUE,
  excludeChr=c("chrM","chrX","chrY"),
  filterQuantile=0.995,
  outlierQuantile=c(0.02,0.98),
  varFeatures=25000,
  totalFeatures=500000,
  scaleTo=10000){

  

    LSI = firstPass(x=x,quantile=quantile,varFeatures=varFeatures,totalFeatures=totalFeatures)



  }
firstPass = function(x,
                     useMatrix="counts",
                     preTfIdf=NULL,
                     binarize=TRUE,
                     excludeChr=c("chrM","chrX","chrY"),
                     filterQuantile=0.995,
                     outlierQuantile=c(0.02,0.98),
                     varFeatures=25000,
                     totalFeatures=500000,
                     scaleTo=10000){

  mat <- filterAndGetMat(x=x, useMatrix=useMatrix, excludeChr=excludeChr,binarize=binarize)

  selectedFeatures = findTopFeatures(mat,
    filterQuantile=filterQuantile,
    varFeatures=varFeatures,
    totalFeatures=totalFeatures)

  if(!is.null(preTfIdf)){
    tfidfRes <- list(tfidf = filterAndGetMat(x=x, 
                                          useMatrix=preTfIdf, 
                                          features=selectedFeatures,
                                          excludeChr=excludeChr,
                                          binarize=FALSE,
                                          replaceZeros=FALSE),
                  idf = metadata(x)[[preTfIdf]][["idf"]],
                  outliers = metadata(x)[[preTfIdf]][["outliers"]])
  } else {
    tfidfRes <- calcTfIdf(mat[rownames(mat) %in% selectedFeatures,], outlierQuantiles=outlierQuantile, excludeZeros=TRUE, scaleTo=scaleTo)
  }  
  return(tfidfRes) 
}

iterativePass = function(){

  return(selectedFeatures)
}
# testCorrelation <- function(pCheck, pCheck2) {
#   lapply(seq_len(ncol(pCheck)), function(x) {
#     cor(pCheck[, x], pCheck2[, x])
#   })
# }