#' Read trees
#'
#' Reads in a tree file containing one or multiple trees
#'
#' Reads in a tree file in either nexus or newick format, and containing a single tree
#' or multiple trees (as in the results of a Bayesian analysis).
#'
#' @param path (single character string; no default) File path to tree.
#' @param format (single character string; nexus) Format of tree file: nexus or newick
#'
#' @return Object of type multiPhylo, with length one if only one tree provided
#'
#' @examples
#' single_tree <- readTrees(path = "data/sub_models/primates_cytb_covariotide_MAP.tre", format = "nexus")
#' multi_trees <- readTrees(path = "data/sub_models/primates_cytb_covariotide.trees", format = "newick")

readTrees <- function(path, format){

  # enforce argument matching

  if (is.character(path) == FALSE) stop("path must be a single character string")
  if (file.exists(path) == FALSE) stop("file does not exist")
  format <- match.arg(format, choices = c("nexus", "newick"))

  # read in tree(s) of type nexus or newick

  if (format == "nexus") {
    tree <- read.nexus(file = path)
  } else if (format == "newick") {
    tree <- read.tree(file = path)
  }

  # convert to type multiPhylo for consistency

  if (class(tree) == "phylo") {
    tree <- c.phylo(tree)
  } else if (class(tree) != "multiPhylo") {
    stop ("tree(s) not of type phylo or multiPhylo")
  }

  return(tree)

}