#' @title Run NMF
#' @description Run nonnegative matrix factorization (NMF) with matrix
#' @param mat a properly matrix data as input for \code{\link[NMF]{nmf}}, you can get it through \code{\link[scNMF]{pre_filter}}.
#' @param rank specification of the factorization rank. Here it must be a single numeric value.
#' @param method specification of the NMF algorithm.
#' @param seed specification of the starting point or seeding method.
#' @param ... Other argument in \code{\link[NMF]{nmf}} function.
#'
#' @importFrom NMF nmf basis coef
#' @return return a list contain w and h matrix of NMF.
#'
#' @export
#'
nmf_programs <- function(mat, rank = 2, method = "snmf/r", seed = 1, ...) {
  ## run nmf
  nmf_res <- nmf(mat, rank = rank, method = method, seed = seed, ...)
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

#' @title Select NMF program features
#' @description Select top n important features for each NMF program
#'
#' @param nmf_programs_list a list of object get use \code{nmf_programs}
#' @param top_n the number of features to select from each NMF program by importance
#'
#' @return return a list for each object with the selected top n gene of each program
#'
#' @export
#'
select_program_features <- function(nmf_programs_list, top_n = 50) {
  genes_w <- lapply(nmf_programs_list, function(x) {
    x$w_basis
  })
  ## get gene programs by NMF score
  top_n_genes <- lapply(genes_w, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:top_n]))

  return(top_n_genes)
}

#' @title select robust program
#' @description Select robust rank program from the NMF result
#' @param top_program_features the object obtained from \code{select_program_features}
#' @param intra_min minimum percents overlap with a program from the same sample (for selecting robust programs), default 0.7
#' @param intra_max maximum percents overlap with a program from the same sample (for removing redundant programs), default 0.2
#' @param inter_filter logical; indicates whether programs should be filtered based on their similarity to programs of other sample, default TRUE
#' @param inter_min minimum overlap with a program from another sample, default 0.2
#'
#' @importFrom stringr str_replace
#'
#' @return return a list of select nmf rank program for each sample. If one sample have
#' no robust program, it will be filtered.
#'
#' @export
#'
select_robust_nmf_program <- function(top_program_features, intra_min = 0.7, intra_max = 0.2, inter_filter = T, inter_min = 0.2) {
  # convert data
  number_of_gene <- sapply(top_program_features, nrow)
  if (length(unique(number_of_gene)) == 1) {
    number_of_gene <- unique(number_of_gene)
  } else {
    number_of_gene <- max(number_of_gene, na.rm = T)
  }
  intra_min <- intra_min * number_of_gene
  intra_max <- intra_max * number_of_gene
  inter_min <- inter_min * number_of_gene

  # Select NMF programs based on the minimum overlap with other NMF programs from the same sample
  intra_intersect <- lapply(top_program_features, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x, y)))))
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))
  nmf_sel <- Map(function(x, y) {
    x[, y >= intra_min]
  }, top_program_features, intra_intersect_max)

  # Select NMF programs based on
  # i) the maximum overlap with other NMF programs from the same sample
  # ii) the minimum overlap with programs from another sample
  nmf_sel <- add_list_name(nmf_sel)

  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  ## calculating intersection between all programs
  inter_intersect <- my.pbapply(nmf_sel_unlist, 2, function(x) apply(nmf_sel_unlist, 2, function(y) length(intersect(x, y))))

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

  ## if one sample have no robust program, set it to NULL
  if (any(sapply(final_filter_list, function(x) {length(x) == 0}))) {
    cli::cli_alert_warning('Have sample with no robust program, set the selected program as {.val {NA}}')
  }
  final_filter_list <- lapply(final_filter_list, function(x) {
    if (length(x)  == 0) {
      NA
    } else {
      x
    }
  })

  return(final_filter_list)
}


#' @title get similarity matrix
#' @description Calculate the jaccard similarity matrix among robust program
#'
#' @param top_program_genes a list of object get use \code{select_top_program_genes}
#' @param robust_program The result of \code{robust_program}
#' @param hclust.method Method to build similarity tree between individual programs
#'
#'
#' @importFrom stats as.dendrogram as.dist hclust order.dendrogram reorder
#'
#' @return return a ordered jaccard similarity matrix
#'
#' @export
#'
getSimilarityMatrix <- function(top_program_genes, robust_program, hclust.method = "average") {
  ## check
  if (!identical(names(top_program_genes), names(robust_program))) {
    cli::cli_abort('The `robust_program` don not match with `top_program_genes`, please check')
  }
  ## select robust genes
  nmf_programs_genes_selected <- Map(function(x, y) {
    res <- x[, is.element(colnames(x), y), drop = F]
    if (ncol(res) == 0) {
      res <- NA
    }
    return(res)
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

  return(nmf_jaccard_ordered)
}


#' Generate meta-program
#'
#' @description
#' Generate meta-program cluster
#'
#' @param similarity_matrix jaccard similarity matrix get use \code{getSimilarityMatrix}
#' @param nmf_programs result get use \code{select_program_features}
#' @param Genes_nmf_w_basis a list of w matrix get use \code{nmf_programs}
#' @param initial_score the minus cutoff jaccard index for defining the first NMF program in a cluster
#' @param min_add_score the minus cutoff jaccard index for adding a new NMF to the forming cluster
#' @param min_group_size the minimal program number to consider for defining a cluster
#'
#' @return return a list with meta-program cluster member and meta-program feature set
#'
#' @export
#'
getMP <- function(similarity_matrix, nmf_programs, Genes_nmf_w_basis,
                  initial_score = 0.2, min_add_score = 0.2,
                  min_group_size = 5) {

  ## overlap programs with similarity_matrix
  inter_program <- intersect(colnames(nmf_programs), colnames(similarity_matrix))
  if (length(inter_program) != ncol(similarity_matrix)) {
    cli::cli_abort('The program in `similarity_matrix` must all in `nmf_programs`, please check')
  }
  nmf_programs <- nmf_programs[, colnames(nmf_programs) %in% colnames(similarity_matrix)]

  ## set w matrix name
  Genes_nmf_w_basis <- add_list_name(Genes_nmf_w_basis)

  ## order programs by similarity with other programs
  Sorted_intersection <- sort(apply(similarity_matrix, 2, function(x) (length(which(x >= initial_score)) - 1)), decreasing = TRUE)

  ## set default param
  Cluster_list <- MP_list <- list()
  k <- 1
  Curr_cluster <- c()
  n_feature <- nrow(nmf_programs)

  ## select program cluster
  while (Sorted_intersection[1] > min_group_size) {
    cli::cli_alert_info('Clculating cluster{.val {k}}')
    ## select one base MP
    Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))

    # ------------ intersection between all remaining NMFs and Genes in MP ------------------ ##
    ## 1. Genes in the forming MP are first chosen to be those in the first NMF.
    ##    Genes_MP always has only 50 genes and evolves during the formation of the cluster
    Genes_MP <- nmf_programs[,names(Sorted_intersection[1])]
    ## remove selected NMF
    nmf_programs <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))]
    ## intersection between all other NMFs and Genes_MP, ordered by jaccard index
    # Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) jaccardIndex(Genes_MP,x)) , decreasing = TRUE)
    Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)

    NMF_history <- Genes_MP

    ## 2. Create gene list is composed of intersecting genes (in descending order by frequency).
    ##    When the number of genes with a given frequency span bewond the 50th genes,
    ##    they are sorted according to their NMF score.
    while (Intersection_with_Genes_MP[1] >= min_add_score) {
      ## add new program to cluster feature
      Curr_cluster <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
      ## overlap with added program
      Genes_MP_temp <- sort(table(c(NMF_history, nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)
      ## genes with overlap equal to the last gene
      Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[n_feature])]

      if (length(Genes_at_border)>1){
        ### Sort last genes in Genes_at_border according to maximal NMF gene scores
        ### Run across all NMF programs in Curr_cluster and extract NMF scores for each gene
        Genes_curr_NMF_score <- c()
        for (i in Curr_cluster) {
          splited_name <- strsplit(i, '_rank_')[[1]]
          curr_study <- splited_name[1]
          nmf_w_basis_idx <- match(names(Genes_at_border), rownames(Genes_nmf_w_basis[[curr_study]]))
          Q <- Genes_nmf_w_basis[[curr_study]][nmf_w_basis_idx[!is.na(nmf_w_basis_idx)] ,i]
          names(Q) <- names(Genes_at_border[!is.na(nmf_w_basis_idx)])  ### sometimes when adding genes the names do not appear
          Genes_curr_NMF_score <- c(Genes_curr_NMF_score, Q)
        }
        ## sort by NMF w score
        Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
        Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]

        Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[n_feature])]), names(Genes_curr_NMF_score_sort))
      } else {
        Genes_MP_temp <- names(Genes_MP_temp)[1:n_feature]
      }

      NMF_history <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])
      Genes_MP <- Genes_MP_temp[1:n_feature]

      # remove selected NMF
      nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]
      # intersection between all other NMFs and Genes_MP
      Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
    }

    ## select one MP
    Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
    MP_list[[paste0("MP_",k)]] <- Genes_MP
    k <- k+1

    # Remove current chosen cluster
    similarity_matrix <- similarity_matrix[-match(Curr_cluster,rownames(similarity_matrix)) , -match(Curr_cluster,colnames(similarity_matrix))]

    # Sort intersection of remaining NMFs not included in any of the previous clusters
    Sorted_intersection <- sort(apply(similarity_matrix, 2, function(x) (length(which(x>=initial_score))-1)) , decreasing = TRUE)
    Curr_cluster <- c()

  }

  if (length(Cluster_list) == 0) {
    cli::cli_abort('Can\'t find cluster')
  } else {
    list(MP = Cluster_list, MP_features = MP_list)
  }

}

