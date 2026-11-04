library(testthat)

describe("getModalityRowStartEnd", {
  
  it("returns correct zero-based indices", {
    
    modalities = rep(
      c("K4", "DEM", "K27"),
      times = c(50, 5, 5)
    )
    
    result = getModalityRowStartEnd(modalities)
    
    expect_equal(result["K4", "start"], 0)
    expect_equal(result["K4", "end"], 49)
    
    expect_equal(result["DEM", "start"], 50)
    expect_equal(result["DEM", "end"], 54)
    
    expect_equal(result["K27", "start"], 55)
    expect_equal(result["K27", "end"], 59)
  })
  
  it("rejects modalities occurring in nonconsecutive blocks", {
    
    modalities = c("K4", "DEM", "K4")
    
    expect_error(
      getModalityRowStartEnd(modalities),
      "Modalities must occur in consecutive blocks"
    )
  })
  
})
describe("getSparseIndices", {
  
  it("selects only indices belonging to the requested modality", {
    
    modalities = rep(
      c("K4", "DEM", "K27"),
      times = c(50, 5, 5)
    )
    
    coordinates = getModalityRowStartEnd(
      modalities = modalities
    )
    indices = c(0L, 49L, 50L, 54L, 55L, 59L)
    
    expect_equal(
      getSparseIndices(
        rowRanges = coordinates,
        sparseRowIndices = indices,
        modality = "K4"
      ),
      c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    )
    
    expect_equal(
      getSparseIndices(
        rowRanges = coordinates,
        sparseRowIndices = indices,
        modality = "DEM"
      ),
      c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
    )
    
    expect_equal(
      getSparseIndices(
        rowRanges = coordinates,
        sparseRowIndices = indices,
        modality = "K27"
      ),
      c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
    )
  })
  
  it("rejects unknown modalities", {
    testData = getTestMat()
    modalities = testData[[2]]
    
    coordinates = getModalityRowStartEnd(
      modalities = modalities
    )
    
    expect_error(
      getSparseIndices(
        rowRanges = coordinates,
        sparseRowIndices = testData[[1]]@i,
        modality = "unknown"
      )
    )
  })
})

describe("getModality", {
  
  it("extracts modality labels after _m", {
    
    v = c(
      "feature_1_some_name_mK4",
      "feature_2_mDEM",
      "feature_3_another_name_mK27"
    )
    
    expected = c("K4", "DEM", "K27")
    
    result = getModality(v = v)
    
    expect_equal(result, expected)
  })
  
  it("rejects row names without a modality delimiter", {
    
    v = c("feature_1_K4", "feature_2_mDEM")
    
    expect_error(
      getModality(v = v)
    )
  })
  
})