

describe("featureSelectionLSI multimodal", {
  it("runs full iterative LSI with Mahalanobis outlier detection", {
    set.seed(123)
    
    nFeatures = 3000
    nCells = 120
    
    modalities = rep(c("K4", "DEM", "K27"), times = c(2500, 250, 250))
    
    groups = rep(1:3, each = 40)
    
    # -----------------------------------------
    # Independent modality-specific depths
    # K4 ~ 10000, DEM ~ 5000, K27 ~ 10000
    # -----------------------------------------
    
    depth_K4 = round(rlnorm(nCells, meanlog = log(10000), sdlog = 0.25))
    
    depth_DEM = round(rlnorm(nCells, meanlog = log(5000), sdlog = 0.25))
    
    depth_K27 = round(rlnorm(nCells, meanlog = log(10000), sdlog = 0.25))
    
    # -----------------------------------------
    # Introduce high and low depth outliers
    # independently for each modality
    # -----------------------------------------
    
    # K4
    depth_K4[1:3] = depth_K4[1:3] * 20
    depth_K4[4:6] = pmax(1, round(depth_K4[4:6] * 0.05))
    
    # DEM
    depth_DEM[7:9] = depth_DEM[7:9] * 15
    depth_DEM[10:12] = pmax(1, round(depth_DEM[10:12] * 0.05))
    
    # K27
    depth_K27[13:15] = depth_K27[13:15] * 18
    depth_K27[16:18] = pmax(1, round(depth_K27[16:18] * 0.05))
    
    highOutlierCells = c(1:3, 7:9, 13:15)
    lowOutlierCells = c(4:6, 10:12, 16:18)
    normalCells = 19:nCells
    
    # -----------------------------------------
    # Simulate three cell populations
    # -----------------------------------------
    
    lambda = matrix(1, nrow = nFeatures, ncol = nCells)
    
    for (group in 1:3) {
      cells = which(groups == group)
      
      features = c(((group - 1) * 300 + 1):(group * 300),
                   2500 + ((group - 1) * 60 + 1):(group * 60),
                   2750 + ((group - 1) * 60 + 1):(group * 60))
      
      lambda[features, cells] =
        lambda[features, cells] + 5
    }
    
    # -----------------------------------------
    # Apply independent sequencing processes
    # Normalize expected counts within each
    # modality to its own sequencing depth
    # -----------------------------------------
    
    for (modality in c("K4", "DEM", "K27")) {
      idx = which(modalities == modality)
      
      depth = switch(modality,
                     K4 = depth_K4,
                     DEM = depth_DEM,
                     K27 = depth_K27)
      
      expected = colSums(lambda[idx, , drop = FALSE])
      
      lambda[idx, ] = sweep(lambda[idx, , drop = FALSE], 2, depth / expected, "*")
    }
    
    # -----------------------------------------
    # Sample Poisson fragment counts
    # -----------------------------------------
    
    mat = matrix(rpois(nFeatures * nCells, lambda = as.vector(lambda)),
                 nrow = nFeatures,
                 ncol = nCells)
    
    mat = as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
    
    rownames(mat) = paste0("feature_", seq_len(nFeatures), "_m", modalities)
    
    colnames(mat) = paste0("Cell_", seq_len(nCells))
    
    # -----------------------------------------
    # Construct SingleCellExperiment
    # -----------------------------------------
    
    x = SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))
    
    # Observed modality-specific sequencing depths
    x$nFrags_K4 = Matrix::colSums(mat[modalities == "K4", , drop = FALSE])
    
    x$nFrags_DEM = Matrix::colSums(mat[modalities == "DEM", , drop = FALSE])
    
    x$nFrags_K27 = Matrix::colSums(mat[modalities == "K27", , drop = FALSE])
    
    x$nFrags = Matrix::colSums(mat)
    
    # -----------------------------------------
    # Verify sequencing depths
    # -----------------------------------------
    
    expect_equal(as.numeric(x$nFrags),
                 as.numeric(x$nFrags_K4 + x$nFrags_DEM + x$nFrags_K27))
    
    # All observed fragment counts are integers
    for (depthCol in c("nFrags_K4", "nFrags_DEM", "nFrags_K27")) {
      expect_true(all(x[[depthCol]] >= 0))
      expect_true(all(x[[depthCol]] %% 1 == 0))
    }
    
    # DEM should have approximately half
    # the sequencing depth of K4 and K27
    # Exclude intentionally introduced outliers
    
    expect_equal(mean(x$nFrags_DEM[normalCells]) /
                   mean(x$nFrags_K4[normalCells]), 0.5, tolerance = 0.15)
    
    expect_equal(mean(x$nFrags_DEM[normalCells]) /
                   mean(x$nFrags_K27[normalCells]), 0.5, tolerance = 0.15)
    
    # -----------------------------------------
    # Mahalanobis outlier detection
    # using all three modalities
    # -----------------------------------------
    
    x = detectOutliers(
      sce = x,
      cols = c("nFrags_K4", "nFrags_DEM", "nFrags_K27"),
      FN = normit,
      quantile_cutoff = 0.95,
      outlierColname = "outliers"
    )
    
    expect_type(x$outliers, "logical")
    expect_length(x$outliers, nCells)
    expect_false(anyNA(x$outliers))
    expect_true(any(x$outliers))
    
    # Verify detection of both high and low
    # depth outliers
    expect_true(any(x$outliers[highOutlierCells]))
    expect_true(any(x$outliers[lowOutlierCells]))
    
    # -----------------------------------------
    # Full iterative multimodal LSI
    # -----------------------------------------
    
    result = featureSelectionLSI(
      x = x,
      iter = 2,
      multimodal = TRUE,
      scaleMethod = "Frobenius",
      outlierCellsCol = "outliers",
      depthCol = c("nFrags_K4", "nFrags_DEM", "nFrags_K27"),
      excludeChr = character(0),
      varFeatures = 1500,
      totalFeatures = 3000,
      nDimensions = 10,
      seed = 123
    )
    
    # -----------------------------------------
    # Verify LSI output
    # -----------------------------------------
    
    expect_s4_class(result, "LinearEmbeddingMatrix")
    
    tfidfRes = result@metadata$tfidf
    tfidf = tfidfRes$tfidf
    
    # Verify Mahalanobis outlier propagation
    expect_setequal(tfidfRes$outliers, colnames(x)[x$outliers])
    
    # Verify feature selection produced results
    expect_gt(length(result@metadata$selectedFeatures), 0)
    
    # Verify all cells are retained
    expect_equal(ncol(tfidf), nCells)
    
    expect_identical(colnames(tfidf), colnames(x))
    
    # Verify finite TF-IDF values
    expect_true(all(is.finite(tfidf@x)))
    
    # -----------------------------------------
    # Verify Frobenius normalization
    # for whichever modalities were retained
    # -----------------------------------------
    
    selectedModalities = getModality(v = rownames(tfidf))
    
    for (modality in unique(selectedModalities)) {
      values = tfidf[selectedModalities == modality, , drop = FALSE]@x
      
      expect_equal(sqrt(sum(values^2)), 1, tolerance = 1e-8)
    }
  })
})
