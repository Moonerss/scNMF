#'
#' @param min_expr The min expression value of a gene
#' @param min_pct The min percent of cells expressed a gene
#' @param verbose Print more info
#' @param ... Arguments passed to other methods
#'
#' @importFrom Seurat RelativeCounts
#'
#' @rdname pre_filter
#' @export
#' @method pre_filter default
#'
pre_filter.default <- function(object, min_expr = 3.5, min_pct = 0.02, verbose = TRUE, ...) {

  if (verbose) cli_alert_info_with_time("convert counts to {cli::col_blue('cpm')}")
  cpm <- RelativeCounts(object, scale.factor = 1e6)

  if (verbose) cli_alert_info_with_time("normalize {cli::col_blue('cpm')} to {cli::col_blue('log2((cpm/10) + 1)')}")
  cpm_log <- log2((cpm / 10) + 1)

  if (verbose) cli_alert_info_with_time("filter gene with {cli::col_blue('min_expr')} and {cli::col_blue('min_pct')}")
  cpm_log <- cpm_log[apply(cpm_log, 1, function(x) length(which(x > min_expr)) > ncol(cpm_log) * min_pct), ]

  if (verbose) cli_alert_info_with_time("scale the data with {cli::col_blue('x - rowMeans(x)')}")
  cpm_log <- cpm_log - rowMeans(cpm_log)

  if (verbose) cli_alert_info_with_time("set negative value to {cli::col_blue('0')}")
  cpm_log[cpm_log < 0] <- 0

  return(cpm_log)
}

#'
#' @import Matrix Matrix
#'
#' @rdname pre_filter
#' @export
#' @method pre_filter Matrix
#'
pre_filter.Matrix <- pre_filter.default

#'
#' @param assay Specific assay to get data from or set data for
#'
#' @importFrom SeuratObject %||% GetAssayData DefaultAssay
#'
#' @rdname pre_filter
#' @export
#' @method pre_filter Seurat
#'
pre_filter.Seurat <- function(object, assay = 'RNA',
                              min_expr = 3.5, min_pct = 0.02,
                              verbose = TRUE,
                              ...) {
  assay <- assay %||% DefaultAssay(object)
  if (!is.element('counts', slotNames(object@assays[[assay]]))) {
    cli::cli_abort('There is no {cli::col_red(\'counts\')} layer in object, please check')
  }

  counts_data <- GetAssayData(object, assay = assay, layer = 'counts')
  final_res <- pre_filter(counts_data, min_expr = min_expr, min_pct = min_pct, verbose = verbose, ...)

  n_genes <- nrow(counts_data) - nrow(final_res)
  if (verbose) cli_alert_info_with_time('Filter {cli::col_blue(n_genes)} gene')

  return(final_res)
}


