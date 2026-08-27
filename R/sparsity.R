#' compute sparsity for a Matrix, SummarizedExperiment, or SingleCellExperiment
#' 
#' @param object  the object whose sparsity we would like to evaluate
#' 
#' @return        1 - (length(object@x)/prod(dim(object))), i.e. (1 - density)
#' 
#' @details       just adding a method for Matrix objects to RcppML's generic.
#'                if a SummarizedExperiment is passed, the first assay is used.
#'
#' @import        RcppML
#' @import        Matrix
#' @import        SummarizedExperiment
#'
#' @rdname        sparsity 
#'
#' @export
#'
#' @name          sparsity
# NULL
# 
# if (!isGeneric("sparsity")) { 
#   setGeneric("sparsity", 
#              function(object) 1 - (length(object@x) / prod(dim(object))))
# }
# 
# setMethod("sparsity", "Matrix", 
#           function(object) 1 - (length(object@x) / prod(dim(object))))
# setMethod("sparsity", "SummarizedExperiment", 
#           function(object) 1 - (length(assay(object)@x) / prod(dim(object))))
setGeneric(
  "sparsity",
  function(object) standardGeneric("sparsity")
)

#' @rdname sparsity
setMethod(
  "sparsity",
  signature(object = "Matrix"),
  function(object) {
    1 - length(object@x) / prod(dim(object))
  }
)

#' @rdname sparsity
setMethod(
  "sparsity",
  signature(object = "SummarizedExperiment"),
  function(object) {
    x <- SummarizedExperiment::assay(object)
    sparsity(x)
  }
)