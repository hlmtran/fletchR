# library(ArchR)
# proj <- getTestProject()
# proj <- addIterativeLSI(proj, dimsToUse=1:5, varFeatures=1000, force=TRUE)

# SCE <- archRtoSCE(proj)

# mat = filterAndGetMat(SCE,
#     useMatrix = "counts",
#     replaceZeros = FALSE,
#     binarize = TRUE,
#     subsetLSI = FALSE,
#     prune = c(1,1))

# selectedFeatures = findTopFeatures(mat,
#     filterQuantile=0.995,
#     varFeatures=25000,
#     totalFeatures=500000)


library(testthat)
library(Matrix)

describe("findTopFeatures()", {
  
  it("throws an error when provided invalid inputs", {
    base_mat <- matrix(1:2000, nrow = 2000, ncol = 1)
    expect_error(
      findTopFeatures(base_mat),
      "inherits"
    )
    
    mat <- Matrix(1:2000, nrow = 2000, ncol = 1, sparse = TRUE)
    
    expect_error(
      findTopFeatures(mat, filterQuantile = 1.5),
      "Quantile must be between 0 and 1."
    )
    expect_error(
      findTopFeatures(mat, filterQuantile = -0.1),
      "Quantile must be between 0 and 1."
    )
    expect_error(
      findTopFeatures(mat, varFeatures = 999),
      "varFeatures must be at least 1000."
    )
  })
  
  it("returns names and applies upper-tail filtering when there are sufficient non-zero features", {
    # Parameters
    n_rows <- 3000
    filter_q <- 0.99
    var_features <- 1000
    total_features <- 10000
    
    # Construct matrix where row sum strictly equals the row index
    mat <- sparseMatrix(i = 1:n_rows, j = rep(1, n_rows), x = 1:n_rows, dims = c(n_rows, 1))
    rownames(mat) <- paste0("gene_", 1:n_rows)
    
    # Expected calculations
    expected_discard_count <- floor((1 - filter_q) * total_features)     # 100
    expected_max_sum <- n_rows - expected_discard_count                 # 3000 - 100 = 2900
    expected_min_sum <- expected_max_sum - var_features + 1             # 2900 - 1000 + 1 = 1901
    
    res <- findTopFeatures(
      mat, 
      filterQuantile = filter_q, 
      varFeatures = var_features, 
      totalFeatures = total_features
    )
    
    expect_type(res, "character")
    expect_length(res, var_features)
    
    kept_sums <- Matrix::rowSums(mat[res, , drop = FALSE])
    expect_equal(max(kept_sums), expected_max_sum)
    expect_equal(min(kept_sums), expected_min_sum)
  })
  
  it("skips upper-tail filtering and warns if non-zero features are insufficient", {
    # Parameters
    n_rows <- 3000
    var_features <- 1000
    threshold_multiplier <- 2.25
    # Keep non-zero count strictly below the threshold (2000 < 2250)
    n_nonzero <- floor(threshold_multiplier * var_features) - 250
    
    mat <- sparseMatrix(i = 1:n_nonzero, j = rep(1, n_nonzero), x = 1:n_nonzero, dims = c(n_rows, 1))
    rownames(mat) <- paste0("gene_", 1:n_rows)
    
    # Since upper-tail filtering is skipped, it takes the top `var_features` directly
    expected_max_sum <- n_nonzero                                        # 2000
    expected_min_sum <- n_nonzero - var_features + 1                     # 2000 - 1000 + 1 = 1001
    
    expect_message(
      res <- findTopFeatures(
        mat, 
        filterQuantile = 0.99, 
        varFeatures = var_features, 
        totalFeatures = 10000
      ),
      "Not enough non-zero features to apply upper-tail filtering."
    )
    
    expect_type(res, "character")
    expect_length(res, var_features)
    
    kept_sums <- Matrix::rowSums(mat[res, , drop = FALSE])
    expect_equal(max(kept_sums), expected_max_sum)
    expect_equal(min(kept_sums), expected_min_sum)
  })
  
  it("strictly removes all zero-count features from the final selection", {
    # Parameters
    n_rows <- 3000
    var_features <- 1000
    n_nonzero <- 500  # Less than var_features
    
    mat <- sparseMatrix(i = 1:n_nonzero, j = rep(1, n_nonzero), x = 1:n_nonzero, dims = c(n_rows, 1))
    rownames(mat) <- paste0("gene_", 1:n_rows)
    
    suppressMessages({
      res <- findTopFeatures(mat, varFeatures = var_features, totalFeatures = 10000)
    })
    
    # Must drop all zero-count rows and only return the available non-zero features
    expect_type(res, "character")
    expect_length(res, n_nonzero)
    
    kept_sums <- Matrix::rowSums(mat[res, , drop = FALSE])
    expect_true(all(kept_sums > 0))
  })
  
  it("returns NULL if the matrix does not have row names", {
    n_rows <- 3000
    mat <- sparseMatrix(i = 1:n_rows, j = rep(1, n_rows), x = 1:n_rows, dims = c(n_rows, 1))
    
    res <- findTopFeatures(mat, filterQuantile = 0.99, varFeatures = 1000, totalFeatures = 10000)
    
    expect_null(res)
  })
  
  it("checks out with ArchR results",{
    mat = readRDS(testthat::test_path("..","testdata","first_pass_binMat.rds"))
    archRTopFeats = readRDS(testthat::test_path("..","testdata","first_pass_topFeatures.rds")) 
    archRTopFeats$end = archRTopFeats$start + 10000 - 1 #tileSize = 10000 for this set
    archRTopFeats = archRTopFeats |> as("GRanges") |> as.character()
    topFeats = findTopFeatures(mat,filterQuantile = .995,varFeatures = 1000)
    expect_all_true(archRTopFeats %in% topFeats)
    })
})
