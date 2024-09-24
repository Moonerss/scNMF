## add name for list element
add_list_name <- function(list_obj) {
  if (is.null(names(list_obj))) {
    cli::cli_alert_warning('There is no name for each object in the list')
    cli::cli_alert_info('Set name with `sample_{.emph index}`')
    names(list_obj) <- paste0('sample_', 1:length(list_obj))
  }
  res <- Map(function(x, y) {
    if (is.vector(x)) {
      if (length(x == 1) & is.na(x)) {
        cli::cli_alert_warning('Element {.val {y}} only have one {.val {NA}}, skip ...')
        x <- NULL
      } else if (is.null(names(x))) {
        cli::cli_alert_warning('There is no name for element {.val {y}}, skip ...')
        x <- NULL
      } else {
        if (grepl(y, names(x))) {
          cli::cli_alert_warning('The name contain element name {.val {y}}, skip ...')
        } else {
          names(x) <- paste0(y, '_', names(x))
        }
      }
    } else if (is.data.frame(x) | is.matrix(x)) {
      if (is.null(colnames(x))) {
        cli::cli_alert_warning('There is no column name for element {.val {y}}, skip ...')
        x <- NULL
      } else {
        if (any(grepl(y, colnames(x)))) {
          cli::cli_alert_info('The name contain element name {.val {y}}, skip ...')
        } else {
          colnames(x) <- paste0(y, '_', colnames(x))
        }
      }
    } else {
      cli::cli_alert_warning('The element {.val {y}} is not vector\\matrix\\data.frame, skip ...')
      x <- NULL
    }
    return(x)
  }, list_obj, names(list_obj))
  res[sapply(res, is.null)] <- NULL
  return(res)
}

# filter low expr data
filter_genes <- function(mat, min_expr = 3.5, min_pct = 0.2) {
  mat[apply(mat, 1, function(x) length(which(x > min_expr)) > ncol(mat)*min_pct),]
}

# scale matrix
scale_data <- function(mat, center = TRUE, scale = FALSE, non_negative = TRUE) {
  #Center and rescale
  mat <- t(scale(Matrix::t(mat), center=center, scale=scale))
  #check scaling with rowsum=0 (gives NaN)
  if (scale) {
    mat[is.na(mat)] <- 0
  }
  if (non_negative) {
    mat[mat<0] <- 0
  }
  return(mat)
}

#Calculate Jaccard Index
jaccardIndex <- function(a, b) {
  index <- length(intersect(a, b))/length(union(a, b))
  return (index)
}

jaccardSimilarity <- function(gene_lists) {
  nprogs <- length(gene_lists)
  J <- matrix(data=0, ncol=nprogs, nrow = nprogs)
  colnames(J) <- names(gene_lists)
  rownames(J) <- names(gene_lists)
  for (i in 1:nprogs) {
    for (j in 1:nprogs) {
      J[i,j] <- jaccardIndex(gene_lists[[i]], gene_lists[[j]])
    }
  }
  return(J)
}

# color
custom_magma <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))

# calculate cell complexity
get_cell_complexity <- function(mat) {
  complexity <- apply(mat, 2, function(y) length(which(y != 0)))
  return(complexity)
}

# use pbapply
my.pbapply <- ifelse(require(pbapply, quietly = T), pbapply, apply)
