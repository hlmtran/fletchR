
library(testthat)

getTestMat = function(){
  
  set.seed(123)
  
  mat <- Matrix::rsparsematrix(
    nrow = 60,
    ncol = 20,
    density = 0.9
  )
  
  modalities <- rep(
    c("K4", "DEM", "K27"),
    times = c(50, 5, 5)
  )
  
  rownames(mat) <- paste0(modalities, "_", seq_len(60))
  colnames(mat) <- paste0("Cell_", seq_len(20))
  
  return(list(mat,modalities))
}

describe("scaleModalities", {
  
  it("preserves matrix dimensions and dimnames", {
    
    testData = getTestMat()
    
    mat = testData[[1]]
    modalities = testData[[2]]
    
    result = scaleModalities(mat, modalities)
    
    expect_equal(dim(result), dim(mat))
    expect_equal(dimnames(result), dimnames(mat))
    expect_s4_class(result, "dgCMatrix")
  })
  
  it("preserves the sparse matrix structure(only scales)", {
    
    testData = getTestMat()
    
    mat = testData[[1]]
    modalities = testData[[2]]
    
    result = scaleModalities(mat, modalities)
    
    expect_identical(result@i, mat@i)
    expect_identical(result@p, mat@p)
    expect_equal(length(result@x), length(mat@x))
  })

  it("matches independently scaled modality matrices using Frobenius", {
    
    testData = getTestMat()
    mat = testData[[1]]
    modalities = testData[[2]]
    
    # Split into three independent modality matrices
    matK4 = mat[modalities == "K4", , drop = FALSE]
    matDEM = mat[modalities == "DEM", , drop = FALSE]
    matK27 = mat[modalities == "K27", , drop = FALSE]
    
    # Scale each matrix independently
    matK4 = matK4 / sqrt(sum(matK4@x^2))
    matDEM = matDEM / sqrt(sum(matDEM@x^2))
    matK27 = matK27 / sqrt(sum(matK27@x^2))
    
    # Combine to create ground truth
    expected = rbind(matK4, matDEM, matK27)
    
    # Test implementation
    result = scaleModalities(
      mat,
      modalities,
      method = "Frobenius"
    )
    
    expect_equal(result, expected, tolerance = 1e-10)
  })
  
  it("matches independently scaled modality matrices using minMax", {
    
    testData = getTestMat()
    mat = testData[[1]]
    modalities = testData[[2]]
    
    matK4 = mat[modalities == "K4", , drop = FALSE]
    matDEM = mat[modalities == "DEM", , drop = FALSE]
    matK27 = mat[modalities == "K27", , drop = FALSE]
    
    matK4@x = (matK4@x - min(matK4@x)) /
      (max(matK4@x) - min(matK4@x))
    
    matDEM@x = (matDEM@x - min(matDEM@x)) /
      (max(matDEM@x) - min(matDEM@x))
    
    matK27@x = (matK27@x - min(matK27@x)) /
      (max(matK27@x) - min(matK27@x))
    
    expected = rbind(matK4, matDEM, matK27)
    
    result = scaleModalities(
      mat,
      modalities,
      method = "minMax"
    )
    
    expect_equal(result, expected, tolerance = 1e-10)
  })
})
