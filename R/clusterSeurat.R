#' Cluster a Matrix Using Seurat
#'
#' Performs graph-based clustering on a matrix using Seurat's
#' \code{FindNeighbors} and \code{FindClusters} functions. The rows of
#' \code{mat} are treated as observations and the columns as dimensions
#' of a PCA embedding.
#'
#' Singleton observations in the shared nearest neighbor (SNN) graph are
#' handled by reordering the input matrix, reconstructing the graph, and
#' restoring the original order of the observations after clustering.
#'
#' @param mat A numeric matrix containing the embedding to cluster. Rows
#'   correspond to observations and columns correspond to dimensions.
#'   Row names are required.
#' @param n.start Integer specifying the number of random starts used by
#'   \code{\link[Seurat]{FindClusters}}. Default is \code{10}.
#' @param resolution Numeric value specifying the clustering resolution
#'   passed to \code{\link[Seurat]{FindClusters}}. Default is \code{2}.
#' @param dims Integer vector specifying the dimensions of \code{mat} to
#'   use when constructing the neighbor graph. If \code{NULL}, all
#'   dimensions are used.
#'
#' @import Seurat
#' @return A named character vector containing cluster assignments. Names
#'   correspond to the row names of \code{mat}, and clusters are labeled
#'   \code{"Cluster1"}, \code{"Cluster2"}, and so on.
#'
#' @export
clusterSeurat = function(mat,n.start = 10,resolution=2,dims=NULL){
  
  if(is.null(dims)){dims=seq_len(ncol(mat))}
  
  seu <- .createSeurat(mat)
  
  seu = Seurat::FindNeighbors(object = seu,reduction="pca",dims = dims)
  
  cS <- Matrix::colSums(seu@graphs$RNA_snn)
  if(cS[length(cS)] == 1){
    
    #Error Handling with Singletons
    idxSingles <- which(cS == 1)
    idxNonSingles <- which(cS != 1)
    
    rn <- rownames(mat) #original order
    mat <- mat[c(idxSingles, idxNonSingles), ,drop = FALSE]
    
    seu <- .createSeurat(mat)
    seu = Seurat::FindNeighbors(object = seu,reduction="pca",dims = dims)
    seu = Seurat::FindClusters(object = seu,n.start=n.start,resolution=resolution)
    
    clust <- .getClusters(seu)
    clust <- clust[rn]
    
  }else{
    seu = Seurat::FindClusters(object = seu,n.start=n.start,resolution=resolution)
    
    clust <- .getClusters(seu)
  }
  
  clust
  
}

#' Create a Temporary Seurat Object
#'
#' Creates a temporary Seurat object and stores the supplied matrix as
#' a PCA dimensional reduction.
#'
#' @param mat A numeric matrix with observations as rows and dimensions
#'   as columns. Row names are used as cell names.
#'
#' @return A \code{Seurat} object containing \code{mat} as the
#'   \code{"pca"} dimensional reduction.
#'
#' @keywords internal
.createSeurat = function(mat){
  
  tmp <- matrix(rnorm(nrow(mat) * 3, 10), ncol = nrow(mat), nrow = 3)
  colnames(tmp) <- rownames(mat)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  
  seu <- suppressWarnings(Seurat::CreateSeuratObject(tmp, project='temp', min.cells=0, min.features=0))
  seu[['pca']] <- Seurat::CreateDimReducObject(embeddings=mat, key='PC_', assay='RNA')
  
  seu
  
}


#' Extract Cluster Assignments from a Seurat Object
#'
#' Extracts the final metadata column from a clustered Seurat object and
#' converts the cluster identifiers to sequential labels of the form
#' \code{"Cluster1"}, \code{"Cluster2"}, and so on.
#'
#' @param seu A clustered \code{Seurat} object.
#'
#' @return A named character vector containing cluster assignments, with
#'   cell names as names.
#'
#' @keywords internal
.getClusters = function(seu){
  
  clust <- seu@meta.data[,ncol(seu@meta.data)]
  clust <- paste0("Cluster",match(clust, unique(clust)))
  names(clust) <- rownames(seu@meta.data)
  
  clust
  
}