describe("removeDepthCorrelatedDims", {
  
  lsi = readRDS(testthat::test_path("..","testdata","first_pass_LSI.rds"))
  nFrags = readRDS(testthat::test_path("..","testdata","nFrags.rds"))
  cutoff = 0.75
  # param = readRDS(testthat::test_path("..","testdata","first_pass_clustParams.rds"))
  
  it("removes the same dimensions for unscaled mat", {
    ArchRCor = lsi$corToDepth$none
    
    mat_a = lsi$matSVD[ ,ArchRCor <= cutoff ]
    mat_f = removeDepthCorrelatedDims(lsi$matSVD,nFrags,cutoff)

    expect_identical(mat_a,mat_f)
    
  })
  
  it("scales dims properly",{
    #fix sweepSparse in the future
    scaledA = ArchR:::.rowZscores(lsi$matSVD)
    scaledF = fletchR:::ArchRRowZscores(lsi$matSVD)
    
    expect_identical(scaledA,scaledF)
    })
  it("removes the same dimensions for scaled mat",{
    ArchRCor = lsi$corToDepth$scaled
    
    mat_a = ArchR:::.rowZscores(lsi$matSVD)[ ,ArchRCor <= cutoff ]
    mat_f = removeDepthCorrelatedDims(ArchRRowZscores(lsi$matSVD),nFrags,cutoff)
    
    expect_identical(mat_a,mat_f)
    })
})