#' turn an ArchR project into a SingleCellExperiment 
#'
#' @param proj          an ArchRproject 
#' @param ...           additional arguments for ArchR::getMatrixFromProject()
#'
#' @return              a SingleCellExperiment 
#'
#' @details Presumes LSI has been done (this will likely be relaxed soon). 
#'          If you are working with single-cell histone or transcription factor
#'          data, you will most likely want to use larger (5kb-50kb) tiles (see 
#'          Janssens et al, Nature Protocols, 2024 for sciCUT&Tag benchmarks).
#'          Recent (2026) versions of this function simply cram a lot of ArchR
#'          outputs into metadata()$output for simplicity.
#'
#' @import scuttle
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, ...){

  if (!require(ArchR)) stop("This function won't work without an ArchR install")
  g <- archRgenome(proj)
  tileSize <- archRtileSize(proj)
  warning("tileSize set to ", tileSize) 
        
  # grab everything (after possibly warning the user about the above) 
  if (!"TileMatrix" %in% getAvailableMatrices(proj)) {
    stop("You need to add a TileMatrix to your project to export it")
  } else { 
    b <- unique(sapply(getArrowFiles(proj), h5read, "TileMatrix/Info/Class"))
    bb <- (b == "Sparse.Binary.Matrix")
    SCE <- as(getMatrixFromProject(proj, "TileMatrix", binarize=bb), 
              "SingleCellExperiment")
    assayNames(SCE) <- "counts"
    metadata(SCE)$tileSize <- tileSize
    rowData(SCE)$start <- rowData(SCE)$start + 1
    rowData(SCE)$end <- rowData(SCE)$start + tileSize - 1
    rowRanges(SCE) <- as(rowData(SCE), "GRanges")
    rownames(SCE) <- as.character(rowRanges(SCE))
    genome(SCE) <- g

  }

  # pull in embeddings etc. 
  keep <- colnames(SCE)

  if ("IterativeLSI" %in% names(proj@reducedDims)) {
   
    message("Copying UMAP parameters to metadata(SCE)$UMAP...")
    metadata(SCE)$LSI <- proj@reducedDims$IterativeLSI 
    reducedDim(SCE, "LSI") <- getReducedDims(proj, "IterativeLSI")[keep, ]
    colnames(reducedDim(SCE, "LSI")) <- 
      paste0("LSI", seq_len(ncol(reducedDim(SCE, "LSI"))))
    
    LSIfeats <- metadata(SCE)$LSI$LSIFeatures
    LSIfeats$end <- LSIfeats$start + metadata(SCE)$tileSize
    LSIfeats$start <- LSIfeats$start + 1
    LSIgr <- as(LSIfeats, "GRanges")
    names(LSIgr) <- as(LSIgr, "character")
    genome(LSIgr) <- g 
    stopifnot(identical(unname(rowSums(assay(subsetByOverlaps(SCE, LSIgr)))),
                        LSIgr$rowSums))
    ol <- findOverlaps(SCE, LSIgr)
    rowData(SCE)$usedForLSI <- FALSE
    rowData(SCE)$usedForLSI[queryHits(ol)] <- TRUE
    # res <- try(addLSI(SCE))
  
  }

  if ("UMAP" %in% names(proj@embeddings)) {
    message("Copying UMAP to reducedDim(SCE, 'UMAP')...")
    metadata(SCE)$UMAP <- proj@embeddings$UMAP
    reducedDim(SCE, "UMAP") <- getEmbedding(proj, "UMAP")[keep, ]
    colnames(reducedDim(SCE, "UMAP")) <- 
      paste0("UMAP", seq_len(ncol(reducedDim(SCE, "UMAP"))))
  }

  # since it's feasible to stack (e.g. DEM + H3K4me + H3K27me)
  mainExpName(SCE) <- "FragmentCounts"
  rowRanges(SCE)$idx <- NULL # irrelevant to us and gets in the way
  rowData(SCE)$assay <- "FragmentCounts" # for binding to any other altExps
  message("You may want to update mcols(SCE)$assay to be more specific.")
  message("(For example, 'DEM' or 'H3K27me3' or 'H3K4me3' or 'ATAC'...)")

  # could also add others if it makes sense here 
  mats <- setdiff(getAvailableMatrices(proj), "TileMatrix")
  for (mat in mats) {
    message("Attempting to add ", mat, " to altExp(SCE, '", mat, "')...")
    if (mat == "GeneScoreMatrix") { 
      altExp(SCE, mat) <- archRgeneScoreMatrix(proj, genome=g)[, keep]
    } else if (mat == "PeakMatrix") { 
      altExp(SCE, mat) <- archRpeakMatrix(proj, genome=g)[, keep]
    } else { 
      altExp(SCE, mat) <- getMatrixFromProject(proj, mat)[, keep]
    }
  }

  message("Copying cellColData...")
  colData(SCE) <- getCellColData(proj)
  SCE <- addTfIdf(SCE, prune=1)
  
  message("Copying metadata...")
  md <- archRmetadata(proj)
  for (i in names(md)) metadata(SCE)[[i]] <- md[[i]]

  message("Done.")
  return(SCE)

}
