
describe("featureSelectionLSI multimodal", {
  
  it("runs full iterative LSI with Mahalanobis outlier detection", {
    set.seed(123)
    
    nFeatures = 3000
    nCells = 120
    
    modalities = rep(
      c("K4", "DEM", "K27"),
      times = c(2500, 250, 250)
    )
    
    # Simulate three cell populations
    mat = matrix(
      rpois(nFeatures * nCells, lambda = 1),
      nrow = nFeatures,
      ncol = nCells
    )
    
    groups = rep(1:3, each = 40)
    
    for (group in 1:3) {
      cells = which(groups == group)
      
      features = c(
        ((group - 1) * 300 + 1):(group * 300),
        2500 + ((group - 1) * 60 + 1):(group * 60),
        2750 + ((group - 1) * 60 + 1):(group * 60)
      )
      
      mat[features, cells] =
        mat[features, cells] +
        matrix(
          rpois(length(features) * length(cells), lambda = 5),
          nrow = length(features)
        )
    }
    
    mat = as(Matrix::Matrix(mat, sparse = TRUE), "dgCMatrix")
    
    rownames(mat) = paste0(
      "feature_", seq_len(nFeatures), "_m", modalities
    )
    colnames(mat) = paste0("Cell_", seq_len(nCells))
    
    # Introduce high-depth outliers
    mat[, 1:3] = mat[, 1:3] * 20
    
    x = SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = mat)
    )
    
    # Modality-specific sequencing depths
    x$nFrags_K4 = Matrix::colSums(
      mat[modalities == "K4", , drop = FALSE]
    )
    
    x$nFrags_DEM = Matrix::colSums(
      mat[modalities == "DEM", , drop = FALSE]
    )
    
    x$nFrags = Matrix::colSums(mat)
    
    # Detect Mahalanobis outliers
    x = detectOutliers(
      sce = x,
      cols = c("nFrags_K4", "nFrags_DEM"),
      FN = normit,
      quantile_cutoff = 0.95,
      outlierColname = "outliers"
    )
    
    expect_type(x$outliers, "logical")
    expect_length(x$outliers, nCells)
    expect_true(any(x$outliers))
    
    # Run full iterative multimodal LSI
    result = featureSelectionLSI(
      x = x,
      iter = 2,
      multimodal = TRUE,
      scaleMethod = "Frobenius",
      outlierCellsCol = "outliers",
      excludeChr = character(0),
      varFeatures = 1500,
      totalFeatures = 3000,
      nDimensions = 10,
      seed = 123
    )
    
    expect_s4_class(result, "LinearEmbeddingMatrix")
    
    tfidfRes = result@metadata$tfidf
    tfidf = tfidfRes$tfidf
    
    # Verify outlier propagation
    expect_setequal(
      tfidfRes$outliers,
      colnames(x)[x$outliers]
    )
    
    # Verify all three modalities are retained
    selectedModalities = getModality(
      v = rownames(tfidf)
    )
    
    expect_setequal(
      unique(selectedModalities),
      c("K4", "DEM", "K27")
    )
    
    # Verify modality-specific Frobenius normalization
    for (modality in unique(selectedModalities)) {
      values = tfidf[
        selectedModalities == modality, ,
        drop = FALSE
      ]@x
      
      expect_equal(
        sqrt(sum(values^2)),
        1,
        tolerance = 1e-8
      )
    }
    
    expect_true(
      length(result@metadata$selectedFeatures) > 0
    )
  })
  
})
