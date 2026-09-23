describe("clusterSeurat", {
  
  rd = readRDS(testthat::test_path("..","testdata","first_pass_matSVD.rds"))
  param = readRDS(testthat::test_path("..","testdata","first_pass_clustParams.rds"))
  
  it("reproduces ArchR clustering", {
    
    clust1 = clusterSeurat(
      mat = rd,
      n.start = param$n.start,
      resolution = param$resolution
    )
    
    clust2 = ArchR:::.clustSeurat(
      mat = rd,
      clustParams = param
    )
    
    expect_identical(clust1,clust2)
    
  })
  
})