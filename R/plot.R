#'
#' @title Plot similarity matrix
#' @description Plot the Jaccard similarity heatmap.
#'
#' @param similarity_matrix jaccard similarity matrix get use \code{getSimilarityMatrix}
#' @param hclust.method cluster method for \code{\link[stats]{hclust}}
#' @param reverse_x reverse the order of x axis.
#' @param reverse_y reverse the order of y axis.
#' @param MP two cloumns data.frame with nmf program and meta-program cluster info,
#' you can get cluster info use \code{getMP}. If you supply this, the heatmap will order by meta-program cluster
#'
#' @importFrom stats as.dendrogram as.dist hclust order.dendrogram reorder
#' @importFrom ggplot2 ggplot aes geom_tile scale_color_gradient2 scale_x_discrete scale_y_discrete scale_fill_gradient2 labs guides guide_colourbar theme element_blank element_rect element_text
#' @importFrom scales squish
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#'
#' @return return a ggplot object
#' @export
#'
plot_similarity_matrix <- function(similarity_matrix, hclust.method = "average", reverse_x = FALSE, reverse_y = FALSE, MP = NULL) {

  ## check
  if (!is.null(MP)) {
    MP <- as.data.frame(MP)
    if (!all(MP[,1] %in% colnames(similarity_matrix))) {
      cli::cli_abort('Not all program in MP current in similarity_matrix, please check')
    }
  }
  # Cluster programs by gene overlap
  tree <- hclust(as.dist(1 - similarity_matrix), method = hclust.method)
  reorder_tree <- reorder(as.dendrogram(tree), rev(colMeans(similarity_matrix)))
  # order jaccard matrix
  nmf_jaccard_ordered <- similarity_matrix[order.dendrogram(reorder_tree), order.dendrogram(reorder_tree)]

  # order by cluster
  if (!is.null(MP)) {
    inds_sorted <- c()
    for (j in unique(MP[,2])){
      cluster_program <- MP[MP[,2] == j, 1]
      inds_sorted <- c(inds_sorted , match(cluster_program , colnames(similarity_matrix)))
    }
    ### clustered NMFs will appear first, and the latter are the NMFs that were not clustered
    inds_new <- c(inds_sorted, which(is.na(match(1:dim(similarity_matrix)[2],inds_sorted))))
    nmf_jaccard_ordered <- nmf_jaccard_ordered[inds_new, inds_new]
  }

  # plot similarity matrix heatmap
  plot_dat <- reshape2::melt(nmf_jaccard_ordered)
  if (reverse_x) plot_dat$Var1 <- factor(plot_dat$Var1, levels = rev(levels(plot_dat$Var1)))
  if (reverse_y) plot_dat$Var2 <- factor(plot_dat$Var2, levels = rev(levels(plot_dat$Var2)))
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
