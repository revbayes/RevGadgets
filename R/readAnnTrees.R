#' Read annotated trees
#'
#' Reads in a tree files containing annotated trees
#'
#' Reads in a tree file in nexus format containing a single or multiple
#' annotated trees. Relies on TreeIO.
#'
#' @param path (single character string; no default) File path to tree(s)
#'
#' @return List of treedata objects ("beastList"), of length one if only one tree provided
#'
#' @examples
#'
#' \dontrun{
#' tree_single_file <- system.file("extdata",
#'     "comp_method_disc/ase_freeK.tree", package="RevGadgets")
#' tree_single <- readAnnTrees(path = tree_single_file)
#'
#' tree_multi_file <- system.file("extdata",
#'     "nexus_multi_ann.nex", package="RevGadgets")
#' tree_multi <- readAnnTrees(path = tree_multi_file)
#'
#' @export
#' @importFrom treeio read.beast


readAnnTrees <- function(path){

  # enforce argument matching

  if (is.character(path) == FALSE) stop("path must be a single character string")
  if (file.exists(path) == FALSE) stop("file does not exist")

  # read in tree(s) using TreeIO

  t <- treeio::read.beast(file = path)

  if (class(t) == "treedata") {
    t <- list(t)
    class(t) <- "beastList"
    }

  return(t)
}
