#' Generate meta-program
#'
#' @description
#' Generate meta-program cluster
#'
#' @param similarity_matrix jaccard similarity matrix get use \code{getSimilarity}
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
