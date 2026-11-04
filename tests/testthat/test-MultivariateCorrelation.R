

describe("correlateDimensions", {
  it("returns a numeric matrix for univariate depth", {
    set.seed(123)
    
    projectedMat = matrix(rnorm(100 * 5), nrow = 100)
    depth = matrix(rpois(100, lambda = 100), ncol = 1)
    
    res = correlateDimensions(projectedMat = projectedMat, depth = depth)
    
    expect_true(is.matrix(res))
    expect_type(res, "double")
    expect_equal(dim(res), c(5L, 1L))
  })
  
  it("returns a numeric matrix for multivariate depth", {
    set.seed(123)
    
    projectedMat = matrix(rnorm(100 * 5), nrow = 100)
    depth = matrix(rpois(100 * 2, lambda = 100), ncol = 2)
    
    res = correlateDimensions(projectedMat = projectedMat, depth =depth)
    
    expect_true(is.matrix(res))
    expect_type(res, "double")
    expect_equal(dim(res), c(5L, 1L))
  })
})