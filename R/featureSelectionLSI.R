featureSelectionLSI = function(x,
  useMatrix="counts",
  preTfIdf=NULL,
  binarize=TRUE,
  excludeChr=c("chrM","chrX","chrY"),
  quantile=0.995,
  outlierQuantile=c(0.02,0.98),
  varFeatures=25000,
  totalFeatures=500000){

  

    LSI = firstPass(x=x,quantile=quantile,varFeatures=varFeatures,totalFeatures=totalFeatures)



  }
firstPass = function(x,preTfIdf=NULL,quantile=0.995,varFeatures=25000,totalFeatures=500000){

  mat <- filterAndGetMat(x=x, useMatrix=useMatrix, excludeChr=excludeChr,binarize=binarize)

  selectedFeatures = findTopFeatures(mat,
    quantile=quantile,
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
    tfidfRes <- calcTfIdf(mat, outlierQuantiles=outlierQuantiles, excludeZeros=TRUE, scaleTo=scaleTo)
  }  

}

iterativePass = function(){

  return(selectedFeatures)
}
# testCorrelation <- function(pCheck, pCheck2) {
#   lapply(seq_len(ncol(pCheck)), function(x) {
#     cor(pCheck[, x], pCheck2[, x])
#   })
# }