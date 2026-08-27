context("Testing ArchR consistency")
        
        
library(ArchR)
proj <- getTestProject()
proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)

SCE <- archRtoSCE(proj)
# SCE = addTfIdf(SCE,useMatrix = "counts",subsetLSI = T,assayName = "TfIdf",binarize=T,outlierQuantiles = c(0,1))
SCE = addLSI(SCE,useMatrix = "ArchRTfIdf",subsetLSI=T,scaleDims = FALSE,nDimensions = 5, corCutOff = 0.75)
c = cor(reducedDim(SCE,"LSI"),SCE@metadata$LSI$matSVD) |> diag()

testthat::expect_all_true(c >= 0.99999)