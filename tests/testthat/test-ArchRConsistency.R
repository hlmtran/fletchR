# context("Testing ArchR consistency")
        
        
library(ArchR)
proj <- getTestProject()
proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)

SCE <- archRtoSCE(proj)
# SCE = addTfIdf(SCE,useMatrix = "counts",subsetLSI = T,assayName = "TfIdf",binarize=T,outlierQuantiles = c(0,1))
SCE = addLSI(SCE,useMatrix = "ArchRTfIdf",subsetLSI=T,scaleDims = FALSE,nDimensions = 5, corCutOff = 0.75)
ArchRLSI = SCE@metadata$LSI

# c = cor(reducedDim(SCE,"LSI")[rownames(ArchRLSI$matSVD),],ArchRLSI$matSVD) |> diag()

logtfidftest = fletchR:::filterAndGetMat(SCE,useMatrix = "ArchRTfIdf",replaceZeros = F)

outliers = rownames(ArchRLSI$matSVD) %in% ArchRLSI$outliers
tt = logtfidftest[,!outliers]
tt = tt[rowSums(tt)>0,]

z = projectSVD(tt,u=ArchRLSI$svd$u,d=ArchRLSI$svd$d,v=ArchRLSI$svd$v,nDim=5)
d = projectSVD(tt,u=ArchRLSI$svd$u,d=ArchRLSI$svd$d,nDim=5)


describe("FletchR",{
  it("Projects RSV identical to ArchR",{testthat::expect_identical(z,ArchRLSI$matSVD[rownames(z),])})
  it("Projects RSV properly with calculated V",{testthat::expect_equal(z, d, tolerance = 1e-5)})
  it("Project RSV with outliers properly",{testthat::expect_equal(reducedDim(SCE,"LSI")[rownames(ArchRLSI$matSVD),],ArchRLSI$matSVD, tolerance = 1e-4)})
})


