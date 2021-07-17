#' dropTip
#'
#' Drop one or multiple tips from your tree
#'
#' Modifies a tree object (in RevGadget's format) by dropping one or more tips
#' from the tree and from any associated data. Wrapper for treeio::drop.tip().
#'
#'
#' @param tree (list of lists of treedata objects; no default) Name of a list of
#' lists of treedata objects, such as produced by readTrees().
#'
#' @param tips (character or numeric, no default) The tips(s) to drop. Either a
#' single taxon name or node number or vector of such.
#'
#' @return returns a list of list of treedata objects, with the modified tips.
#'
#' @seealso treeio: \link[treeio]{drop.tip} and ape: \link[ape]{drop.tip}.
#'
#' @examples
#'
#' file <- system.file("extdata",
#'                     "sub_models/primates_cytb_GTR_MAP.tre",
#'                     package="RevGadgets")
#' tree <- readTrees(paths = file)
#' # reroot tree, then drop the tip
#' tree <- rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")
#' tree_dropped <- dropTip(tree, "Otolemur_crassicaudatus")
#'
#'
#' @export

dropTip <- function(tree, tips) {
  if (!is.list(tree))
    stop("tree should be a list of lists of treedata objects")
  if (class(tree[[1]][[1]]) != "treedata")
    stop("tree should be a list of lists of treedata objects")
  if (class(tips) != "character" & class(tips) != "numeric")
    stop("tips should be of class character or numeric")
  if (length(tips) > length(tree[[1]][[1]]@phylo$tip.label))
    stop("number of tips to drop larger than the number of tips in the tree")
  missing_tips <- tips[ !tips %in% tree[[1]][[1]]@phylo$tip.label ]
  if (length(missing_tips > 0))
    stop(paste0("Tips not found in tree object: ",
                paste0(missing_tips, collapse = ", ")))

  for (i in seq_len(length(tree))) {
    for (j in seq_len(length(tree[[i]]))) {
      t <- tree[[i]][[j]]
      t_dropped <- treeio::drop.tip(t, tip = tips)
      # replace old treedata object with new
      tree[[i]][[j]] <- t_dropped
    }
  }
  return(tree)
}
