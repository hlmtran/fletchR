library(fletchR)

# test project
library(ArchR)
proj <- getTestProject()
proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)
#
# Checking Inputs...
# Detected less than 500 Cells.
# `filterBias` disabled.
# `outlierQuantiles` disabled
# `sampleCellsPre` disabled
# `testBias` in `addClusters` disabled
#
# proj <- addTileMatrix(proj) 
SEsmall <- ArchR::getMatrixFromProject(proj, "TileMatrix", binarize=TRUE)

# convert to SCE so that we can use ReducedDims
library(SingleCellExperiment) 
SCEsmall <- as(SEsmall, "SingleCellExperiment") 
metadata(SCEsmall)$LSI <- proj@reducedDims$IterativeLSI # just a List 

# test whether by naming the features we can robustify .projectLSI
rownames(SCEsmall) 
# NULL 
metadata(SCEsmall)$tileSize <- archRtileSize(proj)
stopifnot(metadata(SCEsmall)$tileSize == metadata(SCEsmall)$LSI$tileSize)
rowData(SCEsmall)$end <- rowData(SCEsmall)$start + metadata(SCEsmall)$tileSize
rowData(SCEsmall)$start <- rowData(SCEsmall)$start + 1
rowRanges(SCEsmall) <- as(rowData(SCEsmall), "GRanges") 
rownames(SCEsmall) <- as.character(rowRanges(SCEsmall))
head(rownames(SCEsmall))

# now do the same with LSI$LSIFeatures 
LSIfeats <- metadata(SCEsmall)$LSI$LSIFeatures
LSIfeats$end <- LSIfeats$start + metadata(SCEsmall)$tileSize
LSIfeats$start <- LSIfeats$start + 1 
LSIgr <- as(LSIfeats, "GRanges") 
names(LSIgr) <- as(LSIgr, "character") 

# are these the features we are looking for?
identical(unname(rowSums(assay(subsetByOverlaps(SCEsmall, LSIgr)))), 
          LSIgr$rowSums)
# [1] TRUE

ol <- findOverlaps(SCEsmall, LSIgr)
rowData(SCEsmall)$usedForLSI <- FALSE
rowData(SCEsmall)$usedForLSI[queryHits(ol)] <- TRUE 
table(rowData(SCEsmall)$usedForLSI)
# 
# FALSE  TRUE 
# 30593  1000 

# need to mask 0-sum rows 
rowData(SCEsmall)$rowSm <- rowSums(assay(SCEsmall))

# test projection and LSI mapping to archRiSEE
# 'TFIDF' not in names(assays(<RangedSummarizedExperiment>))
SCEsmall <- addTfIdf(SCEsmall, prune=1)

# try again
res <- try(addLSI(SCEsmall))

if (!inherits(res, "try-error")) {
  # subsample <- sample(colnames(SCEsmall), 100)
  # toProject <- assay(SCEsmall)[, subsample]
  # LSI <- metadata(SCEsmall)$LSI
  # projectedMatSVD <- archRiSEE::projectLSI(toProject, LSI)
  # # Subsetting TF-IDF matrix...
  # # Running SVD...
  # # Warning in irlba::irlba(mat, nDimensions + 5, nDimensions + 5) :
  # #   convergence criterion below machine epsilon
  # # Warning in irlba::irlba(mat, nDimensions + 5, nDimensions + 5) :
  # #   did not converge--results might be invalid!; try increasing work or maxit
  # # Scaling matSVD...
  # # Checking for depth-correlated columns...
  # # Kept 0% of columns.
  # # Error in matSVD[, toKeep][, seq_len(nDimensions)] : 
  # #   subscript out of bounds
  # testError <- projectedMatSVD - LSI$matSVD[subsample, ] 
  
  LSI <- metadata(SCEsmall)$LSI
  toProject <- assay(SCEsmall)[match(names(LSIgr),rownames(SCEsmall)),,drop=F]
  toProject_minus_outliers=toProject[,-match(LSI$outliers,colnames(toProject))]
  # 
  # identical(LSI$LSIFeatures,featureDF_iter2)
  # identical(LSI$idx,LSI_iter2$idx)
  # identical(colnames(toProject),colnames(mat_iter2_computeLSI))
  # 
  # all.equal(
  #   as.matrix(toProject),
  #   as.matrix(mat_iter2_computeLSI),
  #   check.attributes = FALSE
  # )
  # 
  projectedMatSVD <- archRiSEE::projectLSI(toProject_minus_outliers, LSI)
  
  # identical(projectedMatSVD,LSI$matSVD[rownames(projectedMatSVD),])
  # testCorrelation(projectedMatSVD,LSI$matSVD[rownames(projectedMatSVD),])
  
  cor(projectedMatSVD,LSI$matSVD[rownames(projectedMatSVD),]) |> diag()
  plot(projectedMatSVD,LSI$matSVD[rownames(projectedMatSVD),])
}

###testing tfidf internals
# SCEsmall2 = SCEsmall
# SCEsmall2 <- addTfIdf(SCEsmall2, prune=1,useMatrix = "TileMatrix",subsetLSI = TRUE,outlierQuantiles = c(0,1))
# 
# tfidf_n = SCEsmall2@metadata[["TFIDF"]]
# 
# identical(tfidf_n[,colnames(logTFIDF_iter2)]@x,logTFIDF_iter2@x)
# plot(tfidf_n[,colnames(logTFIDF_iter2)]@x,logTFIDF_iter2@x)
# cor(tfidf_n[,colnames(logTFIDF_iter2)]@x,logTFIDF_iter2@x)
# 
# SCEsmall2 = addLSI(SCEsmall2,scaleDims = FALSE,nDimensions = 5,corCutOff = 1)
# pmSVD = reducedDim(SCEsmall2,"LSI")
# cor(pmSVD,LSI$matSVD)



# #######
# test_mat = toProject
# #test outlier
# idxOutliers = outlierByQuantile(toProject,c(0,1))
# identical(colnames(toProject)[idxOutliers],LSI$outliers)
# 
# 
# #test log(TFIDF)
# toProject_minus_outliers_minus0Row = remove0Rows(toProject_minus_outliers)
# tfidf_test = logTfIdf(
#   toProject_minus_outliers_minus0Row,
#   idf = calcIDF(toProject_minus_outliers_minus0Row)
#   )
# 
# dim(tfidf_test)
# dim(logTFIDF_iter2)
# 
# all.equal(
#   tfidf_test |> as.matrix(),
#   logTFIDF_iter2 |> as.matrix()
#   )
# cor(
#   tfidf_test |> as.matrix(),
#   logTFIDF_iter2 |> as.matrix()
# ) |> diag()
# 
# tfidf_test2 = logTfIdf(
#   toProject_minus_outliers
# )
# cor(
#   tfidf_test2 |> as.matrix(),
#   logTFIDF_iter2 |> as.matrix()
# ) |> diag()
# plot(
#   tfidf_test2@x,
#   logTFIDF_iter2@x
# )
# 
# # test stacking with LSI mapping 
# 
# 
# # test projecting on stacked data
# 
# all.equal(
#   as.matrix(toProject),
#   as.matrix(mat_iter2_computeLSI),
#   check.attributes = FALSE
# )


# now try the actual archRiSEE function 

col_test <- matrix(
  c(
    1, 0, 0,
    1, 1, 0,
    1, 1, 1
  ),
  nrow = 3
)
col_test

testthat::expect_equal(
  pruneCols(col_test, prune = 2),
  c(TRUE, FALSE, FALSE)
)

SCE <- archRtoSCE(proj)
SCE <- addTfIdf(SCE, prune=1,subsetLSI = FALSE,outlierQuantiles = c(0,1),metadataSlot = "noFiltTfIdf")
# identical(rownames(tfidf_n),rownames(SCE@metadata$ArchRTfIdf))
SCE = addLSI(SCE,metadataSlot = "ArchRTfIdf",scaleDims = FALSE,nDimensions = 5,corCutOff = 0.75)
cor(reducedDim(SCE,"LSI"),SCE@metadata$LSI$matSVD) |> diag() 
cor(reducedDim(SCE,"LSI"),SCE@metadata$LSI$matSVD)

sum(getArchRLSIFeatures(SCE) %in% rownames(SCE@metadata$ArchRTfIdf)) == length(rownames(SCE@metadata$ArchRTfIdf))

SCE = addLSI(SCE,metadataSlot = "noFiltTfIdf",features = getArchRLSIFeatures(SCE),scaleDims = FALSE,nDimensions = 5,corCutOff = 0.75)


#####
SCE2 = addTfIdf(SCE,useMatrix = "counts",subsetLSI = T,assayName = "TfIdf",binarize=T,outlierQuantiles = c(0,1))
SCE2 = addLSI(SCE2,useMatrix = "TfIdf",subsetLSI=T,scaleDims = FALSE,nDimensions = 5, corCutOff = 0.75)
cor(reducedDim(SCE2,"LSI"),SCE2@metadata$LSI$matSVD) |> diag()

logtfidftest = fletchR:::filterAndGetMat(SCE2,useMatrix = "TfIdf",replaceZeros = F)

# cor(logtfidftest[rownames(logTFIDF_iter2),colnames(logTFIDF_iter2)],logTFIDF_iter2)
tt = logtfidftest[,colnames(logTFIDF_iter2)]
tt = tt[rowSums(tt)>0,]
cor(tt@x,logTFIDF_iter2@x)
