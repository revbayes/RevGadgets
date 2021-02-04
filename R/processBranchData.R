#' Title
#'
#' @param phy tree in ape format
#'
#' @return a data frame that translates ape node numbers to RevBayes node numbers
#' @export
#'
#' @examples
matchNodes = function(phy) {

  # get some useful info
  num_tips = length(phy$tip.label)
  num_nodes = phy$Nnode
  tip_indexes = 1:num_tips
  node_indexes = num_tips + num_nodes:1

  node_map = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  current_node = phy$Nnode + 2
  k = 1
  t = 1

  while(TRUE) {

    if ( current_node <= num_tips ) {
      node_map$Rev[node_map$R == current_node] = t
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1
    } else {

      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1

      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go right
        current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }
    }

  }

  return(node_map[,1:2])

}


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
#' library(RevGadgets)
#' library(tidytree)
#' library(ggplot2)
#' library(ape)
#'
#' treefile <- system.file("extdata", "bds/primates.tre", package="RevGadgets")
#' logfile <- system.file("extdata", "bds/primates_BDS_rates_truncated.log", package="RevGadgets")
#'
#' branch_data <- read.table(logfile, header = TRUE, sep = "\t")
#' tree <- as.treedata(read.tree(treefile))
#'
#' annotated_tree <- processBranchData(tree, branch_data, summary = "median")
#'
#' p <- ggtree(annotated_tree) +
#'   aes(colour = avg_lambda) +
#'   theme(legend.position=c(0.2,0.80), legend.background=element_blank()) +
#'   scale_color_continuous("Posterior median \nspeciation rate",
#'                          low="blue", high="green")
processBranchData <- function(tree, df, burnin = 0.25,
                              parnames = c("avg_lambda", "avg_mu", "num_shifts"),
                              summary = "median"){

  df <- df[floor(nrow(df)*burnin):nrow(df),]
  tree_tbl <- as_tibble(tree)
  map <- matchNodes(tree@phylo)

  for (item in parnames){
    parameter <- unname(sapply(df[,grepl(item, colnames(df))], summary))[map$Rev]
    tree_tbl[[item]] <- parameter
  }
  tree2 <- as.treedata(tree_tbl)

  return(tree2)
}

