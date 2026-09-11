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
  
  it("preserves original feature row index ordering in the output vector", {
    # Rows 2 and 4 are variable, 1, 3, and 5 are constant
    counts <- matrix(c(
      10, 10,
      0,  50,
      10, 10,
      50,  0,
      10, 10
    ), nrow = 5, byrow = TRUE)
    
    mat <- Matrix(counts, sparse = TRUE)
    rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
    clusters <- c("C1", "C2")
    
    res <- findVariableFeaturesByCluster(mat, clusters, nFeatures = 2, method = "var")
    
    expect_equal(res, c("gene_2", "gene_4"))
  })
  
  it("returns NULL if matrix has no row names", {
    mat <- Matrix(matrix(rpois(20, lambda = 5), nrow = 4, ncol = 5), sparse = TRUE)
    clusters <- rep(c("A", "B"), length.out = ncol(mat))
    
    res <- findVariableFeaturesByCluster(mat, clusters, nFeatures = 2, method = "var")
    
    expect_null(res)
  })
  
})
