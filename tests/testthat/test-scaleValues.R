library(testthat)

describe("Frobenius norm",{
  it("scaleValues() matches R's Frobenius norm calculation", {
    
    set.seed(123)
    mat = Matrix::rsparsematrix(20, 10, density = 0.5)
    
    expect_equal(
      scaleValues(mat@x),
      mat@x / norm(mat, type = "F"),
      tolerance = 1e-10
    )
  })
})

describe("scaleValues", {
  
  it("scales values using Frobenius normalization", {
    
    v = c(3, 4, 5, 8)
    
    expected = v / sqrt(sum(v^2))
    
    result = scaleValues(
      v = v,
      method = "Frobenius"
    )
    
    expect_equal(result, expected, tolerance = 1e-10)
  })
  
  it("scales values using min-max normalization", {
    
    v = c(2, 4, 6, 10)
    
    expected = (v - min(v)) / (max(v) - min(v))
    
    result = scaleValues(
      v = v,
      method = "minMax"
    )
    
    expect_equal(result, expected, tolerance = 1e-10)
  })
  
  it("rejects invalid scaling methods", {
    
    expect_error(
      scaleValues(
        v = c(1, 2, 3),
        method = "invalid"
      ),
      "Invalid method"
    )
  })
  
})