#' @title Plot similarity matrix
#' @description Plot the Jaccard similarity heatmap.
#'
#' @param top_program_genes a list of object get use \code{select_top_program_genes()}
#' @param robust_program The result of \code{robust_program()}
#' @param hclust.method Method to build similarity tree between individual programs
#'
#' @importFrom stats as.dendrogram as.dist hclust order.dendrogram reorder
#' @importFrom ggplot2 ggplot aes geom_tile scale_color_gradient2 scale_x_discrete scale_y_discrete scale_fill_gradient2 labs guides guide_colourbar theme element_blank element_rect element_text
#' @importFrom scales squish
#' @importFrom purrr flatten
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#'
#' @return return a ggplot object
#' @export
#'
plot_similarity_matrix <- function(top_program_genes, robust_program, hclust.method = "average") {
  ## check
  if (!identical(names(top_program_genes), names(robust_program))) {
    cli::cli_abort('The `robust_program` don not match with `top_program_genes`, please check')
  }
  ## select robust genes
  nmf_programs_genes_selected <- Map(function(x, y) {
    x[, is.element(colnames(x), y), drop = F]
  }, top_program_genes, robust_program)
  nmf_programs_genes_selected <- add_list_name(nmf_programs_genes_selected)
  nmf_programs_genes_selected <- lapply(nmf_programs_genes_selected, function(x) {as.list(as.data.frame(x))})
  nmf_programs_genes_selected <- purrr::flatten(nmf_programs_genes_selected)
  ## Jaccard index
  nmf_jaccard <- jaccardSimilarity(nmf_programs_genes_selected)
  # Cluster programs by gene overlap
  tree <- hclust(as.dist(1 - nmf_jaccard), method = hclust.method)
  reorder_tree <- reorder(as.dendrogram(tree), colMeans(nmf_jaccard))
  # order jaccard matrix
  nmf_jaccard_ordered <- nmf_jaccard[order.dendrogram(reorder_tree), order.dendrogram(reorder_tree)]

  # plot similarity matrix heatmap
  plot_dat <- reshape2::melt(nmf_jaccard_ordered)
  if (length(unique(plot_dat$Var1)) > 200) {
    x_y_unit <- 100
  } else {
    x_y_unit <- 10
  }
  p <- ggplot(data = plot_dat, aes(x = .data$Var1, y = .data$Var2, fill = .data$value, color = .data$value)) +
    geom_tile() +
    scale_color_gradient2(limits = c(0.02, 0.25), low = custom_magma[1:111], mid = custom_magma[112:222], high = custom_magma[223:333], midpoint = 0.135, oob = squish, name = "Similarity\n(Jaccard index)") +
    scale_fill_gradient2(limits = c(0.02, 0.25), low = custom_magma[1:111], mid = custom_magma[112:222], high = custom_magma[223:333], midpoint = 0.135, oob = squish, name = "Similarity\n(Jaccard index)") +
    scale_x_discrete(name="\nPrograms", breaks=unique(plot_dat$Var1)[seq(x_y_unit, length(unique(plot_dat$Var1)), by=x_y_unit)], labels= seq(x_y_unit, length(unique(plot_dat$Var1)), by=x_y_unit)) +
    scale_y_discrete(name="\nPrograms", breaks=unique(plot_dat$Var2)[seq(x_y_unit, length(unique(plot_dat$Var2)), by=x_y_unit)], labels= seq(x_y_unit, length(unique(plot_dat$Var2)), by=x_y_unit)) +
    labs(x = "\nPrograms", y = "\nPrograms") +
    guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
    theme(axis.ticks = element_blank(),
          panel.border = element_rect(fill = F),
          panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10, hjust = 0.5),
          legend.justification = "bottom")

  return(p)
}


#' @title select robust program
#' @description Select robust rank program from the NMF result
#' @param top_program_genes the object obtained from \code{select_top_program_genes()}
#' @param intra_min minimum percents overlap with a program from the same sample (for selecting robust programs), default 0.7
#' @param intra_max maximum percents overlap with a program from the same sample (for removing redundant programs), default 0.2
#' @param inter_filter logical; indicates whether programs should be filtered based on their similarity to programs of other sample, default TRUE
#' @param inter_min minimum overlap with a program from another sample, default 0.2
#'
#' @importFrom stringr str_replace
#' @return return a list of select nmf rank program for each sample
#'
#' @export
#'
select_robust_nmf_program <- function(top_program_genes, intra_min = 0.7, intra_max = 0.2, inter_filter = T, inter_min = 0.2) {
  # convert data
  number_of_gene <- sapply(top_program_genes, nrow)
  if (length(unique(number_of_gene)) == 1) {
    number_of_gene <- unique(number_of_gene)
  } else {
    number_of_gene <- max(number_of_gene, na.rm = T)
  }
  intra_min <- intra_min * number_of_gene
  intra_max <- intra_max * number_of_gene
  inter_min <- inter_min * number_of_gene

  # Select NMF programs based on the minimum overlap with other NMF programs from the same sample
  intra_intersect <- lapply(top_program_genes, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y)))))
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))
  nmf_sel <- Map(function(x, y) {
    x[, y >= intra_min]
  }, top_program_genes, intra_intersect_max)

  # Select NMF programs based on
  # i) the maximum overlap with other NMF programs from the same sample
  # ii) the minimum overlap with programs from another sample
  nmf_sel <- add_list_name(nmf_sel)

  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  ## calculating intersection between all programs
  inter_intersect <- apply(nmf_sel_unlist, 2, function(x) apply(nmf_sel_unlist, 2, function(y) length(intersect(x, y))))

  final_filter <- NULL
  for (i in names(nmf_sel)) {
    a <- inter_intersect[grep(i, colnames(inter_intersect), invert = T), grep(i, colnames(inter_intersect))]
    # for each sample, ranks programs based on their maximum overlap with programs of other sample
    b <- sort(apply(a, 2, max), decreasing = T)
    if (isTRUE(inter_filter)) {
      # selects programs with a maximum intersection of at least 10
      b <- b[b >= inter_min]
    }
    if (length(b) > 1) {
      d <- names(b[1])
      for (y in 2:length(b)) {
        if (max(inter_intersect[d, names(b[y])]) <= intra_max) {
          # selects programs iteratively from top-down.
          # Only selects programs that have a intersection smaller than 10 with a previously selected programs
          d <- c(d, names(b[y]))
        }
      }
      final_filter <- c(final_filter, d)
    } else {
      final_filter <- c(final_filter, names(b))
    }
  }

  ## make the result same with input
  final_filter_list <- lapply(names(nmf_sel), function(x) {
    final_filter[grepl(x, final_filter)]
    })
  names(final_filter_list) <- names(nmf_sel)

  final_filter_list <- Map(function(x, y) {
    str_replace(x, pattern = '.*rank_', 'rank_')
  }, final_filter_list, names(final_filter_list))

  return(final_filter_list)
}


#' @title Select top n gene feature
#' @description
#' Select top n gene for each NMF program
#'
#' @param nmf_programs_list a list of object get use \code{nmf_programs()}
#' @param top_n the number of gene to select from each NMF program
#'
#' @return return a list for each object with the selected top n gene of each program
#'
#' @export
#'
select_top_program_genes <- function(nmf_programs_list, top_n = 50) {
  genes_w <- lapply(nmf_programs_list, function(x) {
    x$w_basis
  })
  ## get gene programs by NMF score
  top_n_genes <- lapply(genes_w, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:top_n]))

  return(top_n_genes)
}



#' @title Run NMF
#' @description Run nonnegative matrix factorization (NMF) with matrix
#' @param mat a numeric matrix with genes in row and cells in column.
#' @param lognormalize whether log transformed the data, if TRUE, the data will transformed with `log2(x+1)`
#' @param do_filtering Whether filter the gene by `min_expr` and `min_pct`.
#' @param min_expr The minum expression of gene expression
#' @param min_pct The minum precent of cell expression gene greater than `min_expr`.
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param non_negative Whether fill the negative value with zero
#' @param rank specification of the factorization rank. Here it must be a single numeric value.
#' @param method specification of the NMF algorithm.
#' @param seed specification of the starting point or seeding method.
#' @param ... Other argument in `nmf` function.
#'
#' @importFrom NMF nmf basis coef
#' @return return a list contain w and h matrix of NMF.
#'
#' @export
#'
nmf_programs <- function(mat, lognormalize = T,
                         do_filtering = TRUE, min_expr = 3.5, min_pct = 0.2,
                         center = TRUE, scale = FALSE, non_negative = TRUE,
                         rank = 2, method = "snmf/r", seed = 1, ...) {
  ## log tranformed
  if (lognormalize) {
    cli::cli_alert_info("Log tranformed data by `log2(x + 1)`")
    new_mat <- log2(mat + 1)
  } else {
    new_mat <- mat
  }
  ## filter genes with low expression
  if (do_filtering) {
    cli::cli_alert_info("Filter genes expression with minum value {.val {min_expr}} less than {.val {min_pct * 100}}% cells")
    new_mat <- filter_genes(new_mat, min_expr = min_expr, min_pct = min_pct)
  }
  ## center data
  new_mat <- scale_data(new_mat, center = center, scale = scale, non_negative = non_negative)
  ## run nmf
  nmf_res <- nmf(new_mat, rank = rank, method = method, seed = seed, ...)
  if (length(rank) == 1) {
    ## run in single rank
    w_basis <- basis(nmf_res)
    colnames(w_basis) <- paste0("rank_", rank, ".", 1:rank)
    h_coef <- t(coef(nmf_res))
    colnames(h_coef) <- paste0("rank_", rank, ".", 1:rank)
    nmf_programs_scores <- list(w_basis = w_basis, h_coef = h_coef)
  } else {
    ## run in multiple rank
    w <- NULL
    h <- NULL
    for (i in rank) {
      w_basis <- basis(nmf_res$fit[[as.character(i)]])
      colnames(w_basis) <- paste0("rank_", i, ".", 1:i)
      h_coef <- t(coef(nmf_res$fit[[as.character(i)]]))
      colnames(h_coef) <- paste0("rank_", i, ".", 1:i)
      w <- cbind(w, w_basis)
      h <- cbind(h, h_coef)
    }
    nmf_programs_scores <- list(w_basis = w, h_coef = h)
  }

  return(nmf_programs_scores)
}


#' Extract data matrix from Seurat object
#'
#' Get the gene expression matrix from a Seurat object, optionally centered
#' and/or subset on highly variable genes
#'
#' @param obj Seurat object
#' @param assay Get data matrix from this assay
#' @param slot Get data matrix from this slot (=layer)
#' @param genes List of variable genes to subset the matrix. If NULL, uses all genes
#' @param center Whether to center the data matrix
#' @param scale Whether to scale the data matrix
#' @param non_negative Whether fill the negative value with zero
#'
#' @return Returns a sparse data matrix (cells per genes), subset
#' according to the given parameters
#'
#' @importFrom Seurat GetAssayData
#' @export
#'
getDataMatrix <- function(obj, assay = "RNA", slot = "data", genes = NULL, center = TRUE, scale = FALSE, non_negative = TRUE) {
  mat <- GetAssayData(obj, assay = assay, layer = slot)

  # subset on
  if (!is.null(genes)) mat <- mat[genes, ]

  # scale data
  mat <- scale_data(mat, center = center, scale = scale, non_negative = non_negative)

  return(mat)
}

