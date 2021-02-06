#' processBranchData
#'
#' @param tree a phylogenetic tree in the treedata format
#' @param df data frame with branch specific data
#' @param burnin fraction of the markov-chain to discard
#' @param parnames defaults to c("avg_lambda", "avg_mu", "num_shifts")
#' @param summary function to summarize the continuous parameter. Typically mean or median
#'
#' @return a treedata file with attached branch-specific data
#' @export
#'
#' @examples
#' \dontrun{
#'
#' treefile <- system.file("extdata", "bds/primates.tre", package="RevGadgets")
#' logfile <- system.file("extdata", "bds/primates_BDS_rates_truncated.p", package="RevGadgets")
#'
#' branch_data <- readTrace(logfile)[[1]]
#' tree <- readTrees(paths = treefile)[[1]][[1]]
#'
#' annotated_tree <- processBranchData(tree, branch_data, summary = "median")
#'
#' p <- ggtree(annotated_tree) +
#'   aes(colour = avg_lambda) +
#'   theme(legend.position=c(0.2,0.80), legend.background=element_blank()) +
#'   scale_color_continuous("Posterior median \nspeciation rate",
#'                          low="blue", high="green")
#'
#' # OR
#'
#' p <- plotTree(tree = list(list(annotated_tree)),
#'               node_age_bars = FALSE,
#'               node_pp = F,
#'               tip_labels = FALSE,
#'               color_branch_by = "avg_lambda",
#'               line_width = 0.8) +
#'      ggplot2::theme(legend.position=c(.1, .9))
#'
#' }
processBranchData <- function(tree, df, burnin = 0.25,
                              parnames = c("avg_lambda", "avg_mu", "num_shifts"),
                              summary = "median"){

  df <- df[floor(nrow(df)*burnin):nrow(df),]
  tree_tbl <- tibble::as_tibble(tree)
  map <- matchNodes(tree@phylo)

  for (item in parnames){
    parameter <- unname(sapply(df[,grepl(item, colnames(df))], summary))[map$Rev]
    tree_tbl[[item]] <- parameter
  }
  tree2 <- tidytree::as.treedata(tree_tbl)

  return(tree2)
}

