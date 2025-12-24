#' Pre normalize and filter data
#'
#' Pre normalize and filter data
#'
#' @param object An object
#'
#' @return Return a expression matrix after normalization and filtered for \code{nmf_programs}
#'
#' @rdname pre_filter
#' @export pre_filter
#'
pre_filter <- function(object, ...) {
  UseMethod(generic = 'pre_filter', object = object)
}
