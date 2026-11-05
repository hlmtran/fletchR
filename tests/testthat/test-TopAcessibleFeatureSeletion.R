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
    expect_error(findTopFeatures(base_mat), "inherits")
    
    mat <- Matrix(1:2000,
                  nrow = 2000,
                  ncol = 1,
                  sparse = TRUE)
    
    expect_error(findTopFeatures(mat, filterQuantile = 1.5),
                 "Quantile must be between 0 and 1.")
    expect_error(findTopFeatures(mat, filterQuantile = -0.1),
                 "Quantile must be between 0 and 1.")
    expect_error(findTopFeatures(mat, varFeatures = 999),
                 "varFeatures must be at least 1000.")
  })
  
  it(
    "returns names and applies upper-tail filtering when there are sufficient non-zero features",
    {
      # Parameters
      n_rows <- 3000
      filter_q <- 0.99
      var_features <- 1000
      total_features <- 10000
      
      # Construct matrix where row sum strictly equals the row index
      mat <- sparseMatrix(
        i = 1:n_rows,
        j = rep(1, n_rows),
        x = 1:n_rows,
        dims = c(n_rows, 1)
      )
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
    }
  )
  
  it("skips upper-tail filtering and warns if non-zero features are insufficient",
     {
       # Parameters
       n_rows <- 3000
       var_features <- 1000
       threshold_multiplier <- 2.25
       # Keep non-zero count strictly below the threshold (2000 < 2250)
       n_nonzero <- floor(threshold_multiplier * var_features) - 250
       
       mat <- sparseMatrix(
         i = 1:n_nonzero,
         j = rep(1, n_nonzero),
         x = 1:n_nonzero,
         dims = c(n_rows, 1)
       )
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
  
  it("strictly removes all zero-count features from the final selection",
     {
       # Parameters
       n_rows <- 3000
       var_features <- 1000
       n_nonzero <- 500  # Less than var_features
       
       mat <- sparseMatrix(
         i = 1:n_nonzero,
         j = rep(1, n_nonzero),
         x = 1:n_nonzero,
         dims = c(n_rows, 1)
       )
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
    mat <- sparseMatrix(
      i = 1:n_rows,
      j = rep(1, n_rows),
      x = 1:n_rows,
      dims = c(n_rows, 1)
    )
    
    res <- findTopFeatures(
      mat,
      filterQuantile = 0.99,
      varFeatures = 1000,
      totalFeatures = 10000
    )
    
    expect_null(res)
  })
  
  it("checks out with ArchR results", {
    mat = readRDS(testthat::test_path("..", "testdata", "first_pass_binMat.rds"))
    archRTopFeats = readRDS(testthat::test_path("..", "testdata", "first_pass_topFeatures.rds"))
    archRTopFeats$end = archRTopFeats$start + 10000 - 1 #tileSize = 10000 for this set
    archRTopFeats = archRTopFeats |> as("GRanges") |> as.character()
    topFeats = findTopFeatures(mat, filterQuantile = .995, varFeatures = 1000)
    expect_all_true(archRTopFeats %in% topFeats)
  })
  
  
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
      expect_error(findTopFeatures(base_mat), "inherits")
      
      mat <- Matrix(1:2000,
                    nrow = 2000,
                    ncol = 1,
                    sparse = TRUE)
      
      expect_error(findTopFeatures(mat, filterQuantile = 1.5),
                   "Quantile must be between 0 and 1.")
      expect_error(findTopFeatures(mat, filterQuantile = -0.1),
                   "Quantile must be between 0 and 1.")
      expect_error(findTopFeatures(mat, varFeatures = 999),
                   "varFeatures must be at least 1000.")
    })
    
    it(
      "returns names and applies upper-tail filtering when there are sufficient non-zero features",
      {
        # Parameters
        n_rows <- 3000
        filter_q <- 0.99
        var_features <- 1000
        total_features <- 10000
        
        # Construct matrix where row sum strictly equals the row index
        mat <- sparseMatrix(
          i = 1:n_rows,
          j = rep(1, n_rows),
          x = 1:n_rows,
          dims = c(n_rows, 1)
        )
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
      }
    )
    
    it("skips upper-tail filtering and warns if non-zero features are insufficient",
       {
         # Parameters
         n_rows <- 3000
         var_features <- 1000
         threshold_multiplier <- 2.25
         # Keep non-zero count strictly below the threshold (2000 < 2250)
         n_nonzero <- floor(threshold_multiplier * var_features) - 250
         
         mat <- sparseMatrix(
           i = 1:n_nonzero,
           j = rep(1, n_nonzero),
           x = 1:n_nonzero,
           dims = c(n_rows, 1)
         )
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
    
    it("strictly removes all zero-count features from the final selection",
       {
         # Parameters
         n_rows <- 3000
         var_features <- 1000
         n_nonzero <- 500  # Less than var_features
         
         mat <- sparseMatrix(
           i = 1:n_nonzero,
           j = rep(1, n_nonzero),
           x = 1:n_nonzero,
           dims = c(n_rows, 1)
         )
         rownames(mat) <- paste0("gene_", 1:n_rows)
         
         suppressMessages({
           res <- findTopFeatures(mat,
                                  varFeatures = var_features,
                                  totalFeatures = 10000)
         })
         
         # Must drop all zero-count rows and only return the available non-zero features
         expect_type(res, "character")
         expect_length(res, n_nonzero)
         
         kept_sums <- Matrix::rowSums(mat[res, , drop = FALSE])
         expect_true(all(kept_sums > 0))
       })
    
    it("returns NULL if the matrix does not have row names", {
      n_rows <- 3000
      mat <- sparseMatrix(
        i = 1:n_rows,
        j = rep(1, n_rows),
        x = 1:n_rows,
        dims = c(n_rows, 1)
      )
      
      res <- findTopFeatures(
        mat,
        filterQuantile = 0.99,
        varFeatures = 1000,
        totalFeatures = 10000
      )
      
      expect_null(res)
    })
    
    it("checks out with ArchR results", {
      mat = readRDS(testthat::test_path("..", "testdata", "first_pass_binMat.rds"))
      archRTopFeats = readRDS(testthat::test_path("..", "testdata", "first_pass_topFeatures.rds"))
      archRTopFeats$end = archRTopFeats$start + 10000 - 1 #tileSize = 10000 for this set
      archRTopFeats = archRTopFeats |> as("GRanges") |> as.character()
      topFeats = findTopFeatures(mat, filterQuantile = .995, varFeatures = 1000)
      expect_all_true(archRTopFeats %in% topFeats)
    })
  })
  
  #multimodal balancing
  
  it("selects features equally across three modalities", {
    nFeatures = 3000
    
    modalities = rep(c("K4", "DEM", "K27"), times = c(1000, 1000, 1000))
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = 1:nFeatures,
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
    
    res = findTopFeatures(
      mat,
      filterQuantile = 0.995,
      varFeatures = 1500,
      totalFeatures = 3000,
      balance = TRUE
    )
    
    selectedModalities = getModality(res)
    
    expect_type(res, "character")
    expect_length(res, 1500)
    
    expect_equal(rle(selectedModalities)$lengths, c(500L, 500L, 500L))
  })
  
  
  it("handles modalities with unequal numbers of features", {
    modalities = rep(c("K4", "DEM", "K27"), times = c(2500, 250, 250))
    
    nFeatures = length(modalities)
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = 1:nFeatures,
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
    
    res = suppressMessages(
      findTopFeatures(
        mat,
        filterQuantile = 0.995,
        varFeatures = 1500,
        totalFeatures = 3000,
        balance = TRUE
      )
    )
    
    selectedModalities = getModality(res)
    
    expect_type(res, "character")
    
    expect_equal(length(res), 1000)
    
    expect_equal(rle(selectedModalities)$lengths, c(500L, 250L, 250L))
    
    expect_true(all(res %in% rownames(mat)))
  })
  
  
  it("applies upper-tail filtering independently to each modality",
     {
       nPerModality = 3000
       
       modalities = rep(c("K4", "DEM", "K27"), each = nPerModality)
       
       nFeatures = length(modalities)
       
       # Strictly increasing counts within each modality
       mat = sparseMatrix(
         i = 1:nFeatures,
         j = rep(1, nFeatures),
         x = rep(seq_len(nPerModality), 3),
         dims = c(nFeatures, 1)
       )
       
       rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
       
       res = findTopFeatures(
         mat,
         filterQuantile = 0.99,
         varFeatures = 3000,
         totalFeatures = 9000,
         balance = TRUE
       )
       
       # Each modality requests 1000 features.
       # Each modality independently discards its top 30.
       expectedMax = 2970
       expectedMin = 1971
       
       expect_length(res, 3000)
       
       for (modality in c("K4", "DEM", "K27")) {
         selected = res[getModality(res) == modality]
         
         values = Matrix::rowSums(mat[selected, , drop = FALSE])
         
         expect_length(selected, 1000)
         
         expect_equal(max(values), expectedMax)
         
         expect_equal(min(values), expectedMin)
       }
     })
  
  
  it("preserves original matrix feature order after selection", {
    nFeatures = 3000
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = rev(seq_len(nFeatures)),
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures))
    
    res = findTopFeatures(
      mat,
      filterQuantile = 1,
      varFeatures = 1000,
      totalFeatures = 3000
    )
    
    # Selected features should follow original row order
    expected = rownames(mat)[1:1000]
    
    expect_identical(res, expected)
  })
  
  
  it("preserves original feature order with balanced selection", {
    modalities = rep(c("K4", "DEM", "K27"), each = 1000)
    
    nFeatures = length(modalities)
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = rev(seq_len(nFeatures)),
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
    
    res = findTopFeatures(
      mat,
      filterQuantile = 1,
      varFeatures = 1500,
      totalFeatures = 3000,
      balance = TRUE
    )
    
    expected = rownames(mat)[rownames(mat) %in% res]
    
    expect_identical(res, expected)
    expect_length(res, 1500)
  })
  
  
  it("removes zero-count features independently within modalities",
     {
       nPerModality = 1000
       
       modalities = rep(c("K4", "DEM", "K27"), each = nPerModality)
       
       nFeatures = length(modalities)
       
       # K4: 1000 nonzero
       # DEM: 300 nonzero
       # K27: 100 nonzero
       nonzero = c(1:1000, 1001:1300, 2001:2100)
       
       mat = sparseMatrix(
         i = nonzero,
         j = rep(1, length(nonzero)),
         x = rep(1, length(nonzero)),
         dims = c(nFeatures, 1)
       )
       
       rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
       
       res = suppressMessages(
         findTopFeatures(
           mat,
           filterQuantile = 0.995,
           varFeatures = 1500,
           totalFeatures = 3000,
           balance = TRUE
         )
       )
       
       selectedModalities = getModality(res)
       
       expect_equal(rle(selectedModalities)$lengths, c(500L, 300L, 100L))
       
       expect_length(res, 900)
       
       expect_true(all(Matrix::rowSums(mat[res, , drop = FALSE]) > 0))
     })
  
  
  it("rejects invalid modality suffixes in balanced selection", {
    nFeatures = 3000
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = rep(1, nFeatures),
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures))
    
    expect_error(
      findTopFeatures(mat, varFeatures = 1500, balance = TRUE),
      "Feature names must end with a valid modality suffix"
    )
  })
  
  
  it("allows fewer than 1000 total features when balance is TRUE", {
    modalities = rep(c("K4", "DEM", "K27"), each = 1000)
    
    nFeatures = length(modalities)
    
    mat = sparseMatrix(
      i = 1:nFeatures,
      j = rep(1, nFeatures),
      x = rep(1, nFeatures),
      dims = c(nFeatures, 1)
    )
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
    
    res = suppressMessages(
      findTopFeatures(
        mat,
        filterQuantile = 1,
        varFeatures = 900,
        totalFeatures = 3000,
        balance = TRUE
      )
    )
    
    expect_length(res, 900)
    
    expect_equal(rle(getModality(res))$lengths, c(300L, 300L, 300L))
  })
  
  
  it("confirms that balancing does not alter feature selection when only one modality is present",
     {
       nFeatures = 3000
       
       mat = sparseMatrix(
         i = 1:nFeatures,
         j = rep(1, nFeatures),
         x = seq_len(nFeatures),
         dims = c(nFeatures, 1)
       )
       
       rownames(mat) = paste0("feature_", seq_len(nFeatures), "_mK4")
       
       unbalanced = findTopFeatures(
         mat,
         filterQuantile = 0.995,
         varFeatures = 1000,
         totalFeatures = 3000,
         balance = FALSE
       )
       
       balanced = findTopFeatures(
         mat,
         filterQuantile = 0.995,
         varFeatures = 1000,
         totalFeatures = 3000,
         balance = TRUE
       )
       
       expect_identical(balanced, unbalanced)
     })
  
  it("enforces equal feature allocation across modalities during the first LSI pass", {
    
    set.seed(123)
    
    nCells = 60
    modalities = rep(
      c("K4", "DEM", "K27"),
      times = c(2500, 1000, 1000)
    )
    
    nFeatures = length(modalities)
    
    mat = Matrix::rsparsematrix(
      nrow = nFeatures,
      ncol = nCells,
      density = 0.15
    )
    
    # Use nonnegative counts
    mat@x = rep(1, length(mat@x))
    
    rownames(mat) = paste0(
      "feature_", seq_len(nFeatures), "_m", modalities
    )
    colnames(mat) = paste0("cell_", seq_len(nCells))
    
    x = SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = mat)
    )
    
    res = .firstPass(
      x = x,
      excludeChr = character(0),
      nDimensions = 10,
      filterQuantile = 1,
      varFeatures = 1500,
      totalFeatures = nFeatures,
      multimodal = TRUE,
      scaleMethod = "Frobenius",
      balance = TRUE,
      seed = 123
    )
    
    selectedFeatures = res@metadata$selectedFeatures
    selectedModalities = getModality(selectedFeatures)
    
    expect_s4_class(res, "LinearEmbeddingMatrix")
    expect_length(selectedFeatures, 1500)
    
    expect_equal(
      rle(selectedModalities)$lengths,
      c(500L, 500L, 500L)
    )
    
    expect_true(all(selectedFeatures %in% rownames(mat)))
  })
})
