#' Extract data matrix from Seurat object or expression matrix
#'
#' Get the gene expression matrix from a Seurat object or an expression matrix, optionally centered
#' and/or subset on highly variable genes
#'
#' @param obj Seurat object or a \code{\link[base]{matrix}}-like object
#' @param genes List of variable genes to subset the matrix. If NULL, uses all genes
#' @param do_filtering Whether filter the gene by `min_expr` and `min_pct`.
#' @param min_expr The minum expression of gene expression
#' @param min_pct The minum precent of cell expression gene greater than `min_expr`.
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param non_negative Whether fill the negative value with zero
#' @param verbose Show more message
#' @param \dots Further arguments passed to each method
#'
#' @return Returns a sparse data matrix (cells per genes), subset according to the given parameters
#'
#' @rdname getDataMatrix
#'
#' @export
#'
getDataMatrix <- function(
  obj,
  ...
) {
  UseMethod(generic = 'getDataMatrix', object = obj)
}
