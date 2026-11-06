library(testthat)
library(Matrix)
library(sparseMatrixStats)

describe("findVariableFeaturesByCluster()", {
  
  it("rejects non-Matrix inputs and unsupported methods", {
    dense_mat <- matrix(1:100, nrow = 10, ncol = 10)
    clusters <- rep(c("A", "B"), each = 5)
    
    expect_error(
      findVariableFeaturesByCluster(dense_mat, clusters),
      "inherits\\(mat, \"Matrix\"\\)"
    )
    
    sparse_mat <- Matrix(dense_mat, sparse = TRUE)
    expect_error(
      findVariableFeaturesByCluster(sparse_mat, clusters, method = "cv"),
      "Method must be 'var' or 'vmr'"
    )
  })
  
  it("returns a character vector of top features with method = 'var'", {
    n_features <- 4
    n_cells <- 4
    
    counts <- matrix(c(
      0,  0, 100, 100,
      50, 50, 50,  50,
      0,  0,   0,   0,
      0,  0,   0,   0
    ), nrow = n_features, byrow = TRUE)
    
    mat <- Matrix(counts, sparse = TRUE)
    rownames(mat) <- paste0("feature_", seq_len(n_features))
    clusters <- c("A", "A", "B", "B")
    
    top_k <- 1
    res <- findVariableFeaturesByCluster(mat, clusters, nFeatures = top_k, method = "var")
    
    expect_type(res, "character")
    expect_length(res, top_k)
    expect_equal(res, "feature_1")
  })
  
  it("returns top features correctly when using method = 'vmr'", {
    counts <- matrix(c(
      0,  50,
      25, 25,
      0,   0
    ), nrow = 3, byrow = TRUE)
    
    mat <- Matrix(counts, sparse = TRUE)
    rownames(mat) <- paste0("feature_", seq_len(nrow(mat)))
    clusters <- c("Cluster1", "Cluster2")
    
    top_k <- 1
    res <- findVariableFeaturesByCluster(mat, clusters, nFeatures = top_k, method = "vmr")
    
    expect_type(res, "character")
    expect_length(res, top_k)
    expect_equal(res, "feature_1")
  })

  it("returns NULL if matrix has no row names", {
    mat <- Matrix(matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5), sparse = TRUE)
    clusters <- rep(c("A", "B"), length.out = ncol(mat))
    
    res <- findVariableFeaturesByCluster(mat, clusters, nFeatures = 2, method = "var")
    
    expect_null(res)
  })
  
  it("checks out with ArchR",{
    binMat = readRDS(testthat::test_path("..","testdata","first_pass_binMat.rds"))
    second_pass_variableFeatures = readRDS(testthat::test_path("..","testdata","second_pass_variableFeatures.rds")) 
    second_pass_variableFeatures$end = second_pass_variableFeatures$start + 10000 - 1 #tileSize = 10000 for this set
    ArchRVarFeat = as.character(as(second_pass_variableFeatures,"GRanges")) 
    
    rd = readRDS(testthat::test_path("..","testdata","first_pass_matSVD.rds"))
    param = readRDS(testthat::test_path("..","testdata","first_pass_clustParams.rds"))
    nFeat = 1000
    
    clust = clusterSeurat(
      mat = rd,
      n.start = param$n.start,
      resolution = param$resolution
    )
    
    varFeat = findVariableFeaturesByCluster(binMat,clust,nFeatures = nFeat)
    
    expect_true(
      identical(ArchRVarFeat,varFeat)
    )
    expect_false(
      identical(varFeat,head(rownames(binMat),nFeat))
      )
    })
})

describe("balanced variable feature selection", {
  
  it("applies balanced feature selection when balance = TRUE", {
    set.seed(123)
    
    modalities = rep(c("K4", "DEM", "K27"), times = c(100, 80, 60))
    
    counts = matrix(
      sample(1:100, length(modalities) * 12, replace = TRUE),
      nrow = length(modalities)
    )
    
    mat = Matrix::Matrix(counts, sparse = TRUE)
    rownames(mat) = paste0("feature_", seq_len(nrow(mat)), "_m", modalities)
    
    clusters = rep(c("A", "B", "C"), each = 4)
    
    res = findVariableFeaturesByCluster(
      mat, clusters, nFeatures = 90, balance = TRUE
    )
    
    expect_length(res, 90)
    expect_equal(
      rle(getModality(res))$lengths,
      c(30L, 30L, 30L)
    )
  })
  
  it("uses global feature selection when balance = FALSE", {
    counts = matrix(c(
      0, 100, 0, 100,
      50, 50, 50, 50,
      100, 0, 100, 0,
      20, 20, 20, 20
    ), nrow = 4, byrow = TRUE)
    
    mat = Matrix::Matrix(counts, sparse = TRUE)
    rownames(mat) = paste0("feature_", seq_len(nrow(mat)))
    
    clusters = c("A", "A", "B", "B")
    
    res = findVariableFeaturesByCluster(
      mat, clusters, nFeatures = 2, balance = FALSE
    )
    
    # No modality suffixes: this would fail if the balancing helper were called
    expect_length(res, 2)
    expect_true(all(res %in% rownames(mat)))
  })
  
})


describe("balanceTopVarFeatures()", {
  
  it("selects the highest-variance features independently within each modality", {
    vars = c(
      feature_1_mK4 = 2,
      feature_2_mK4 = 10,
      feature_3_mK4 = 8,
      feature_4_mDEM = 1,
      feature_5_mDEM = 9,
      feature_6_mDEM = 7,
      feature_7_mK27 = 5,
      feature_8_mK27 = 3,
      feature_9_mK27 = 6
    )
    
    res = balanceTopVarFeatures(vars, nFeatures = 6)
    
    expect_setequal(res, c(2L, 3L, 5L, 6L, 7L, 9L))
    expect_length(res, 6)
  })
  
  it("does not redistribute unused feature allocations across modalities", {
    vars = c(
      feature_1_mK4 = 10,
      feature_2_mK4 = 8,
      feature_3_mK4 = 6,
      feature_4_mK4 = 4,
      feature_5_mDEM = 9,
      feature_6_mK27 = 7,
      feature_7_mK27 = 5,
      feature_8_mK27 = 3
    )
    
    expect_warning(
      object = {res = balanceTopVarFeatures(vars, nFeatures = 6)},
      regexp = "Unequal numbers of features selected across modalities"
    )
    
    # Two allocated per modality, but DEM has only one feature
    expect_setequal(res, c(1L, 2L, 5L, 6L, 7L))
    expect_length(res, 5)
  })
  
  it("rounds down feature allocation when nFeatures is not divisible by the number of modalities", {
    vars = c(
      feature_1_mK4 = 10,
      feature_2_mK4 = 5,
      feature_3_mDEM = 9,
      feature_4_mDEM = 4,
      feature_5_mK27 = 8,
      feature_6_mK27 = 3
    )
    
    res = balanceTopVarFeatures(vars, nFeatures = 5)
    
    # floor(5 / 3) = 1 feature per modality
    expect_setequal(res, c(1L, 3L, 5L))
    expect_length(res, 3)
  })
  
  it("rejects feature names without valid modality suffixes", {
    vars = c(
      feature_1_mK4 = 10,
      feature_2 = 8,
      feature_3_mDEM = 6
    )
    
    expect_error(
      balanceTopVarFeatures(vars, nFeatures = 2),
      "Feature names must end with a valid modality suffix"
    )
  })
  
  it("selects the highest-ranked features when only one modality is present", {
    vars = c(
      feature_1_mK4 = 2,
      feature_2_mK4 = 10,
      feature_3_mK4 = 5,
      feature_4_mK4 = 8
    )
    
    res = balanceTopVarFeatures(vars, nFeatures = 2)
    
    expect_setequal(res, c(2L, 4L))
  })
  
})


