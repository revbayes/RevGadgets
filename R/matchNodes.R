#' match Nodes
#'
#' @param phy (tree in ape format; no default) Tree on which to match nodes
#'
#' @return a data frame that translates ape node numbers to RevBayes node
#' numbers
#'
#' @examples
#'
#' treefile <- system.file("extdata", "bds/primates.tre", package="RevGadgets")
#' tree <- readTrees(treefile)
#' map <- matchNodes(tree[[1]][[1]]@phylo)
#'
#' @export

matchNodes <- function(phy) {
  # get some useful info
  num_tips <- length(phy$tip.label)
  num_nodes <- phy$Nnode
  tip_indexes <- 1:num_tips
  node_indexes <- num_tips + num_nodes:1

  node_map <-
    data.frame(R = 1:(num_tips + num_nodes),
               Rev = NA,
               visits = 0)
  current_node <- phy$Nnode + 2
  k <- 1
  t <- 1

  while (TRUE) {
    if (current_node <= num_tips) {
      node_map$Rev[node_map$R == current_node] <- t
      current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
      t <- t + 1
    } else {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1

      if (node_map$visits[node_map$R == current_node] == 1) {
        # go right
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][2]
      } else if (node_map$visits[node_map$R == current_node] == 2) {
        # go left
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][1]
      } else if (node_map$visits[node_map$R == current_node] == 3) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
    }
  }

  return(node_map[, 1:2])

}


