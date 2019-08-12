#' Read trees
#'
#' Reads in a tree file containing one or multiple trees
#'
#' Reads in a tree file in either nexus or newick format, and containing a single tree
#' or multiple trees (as in the results of a Bayesian analysis).
#'
#' @param paths (vector of character strings; no default) File path(s) to tree(s).
#' @param tree_name (character string; default psi) Name of the tree variable.
#' @param burnin (single numeric value; default = 0.1) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1).
#' @param verbose (logical; default true) Display a status bar?
#'
#' @return A list (across runs) of lists (across samples) of treedata objects.
#'
#' @examples
#'
#' \dontrun{
#'}
#'
#' @export

readTrees <- function(paths, tree_name =  "psi", burnin = 0.1, verbose =  TRUE, ...) {

  # enforce argument matching
  character_paths_are_strings <- is.character(paths)
  if ( any(character_paths_are_strings == FALSE) == TRUE ) {
    # print out the ones that are not character strings
    cat( "Some paths are not character strings:",
         paste0("\t",paths[character_paths_are_strings == FALSE]), sep="\n")
    stop()
  }

  do_files_exist <- file.exists(paths)
  if ( any(do_files_exist == FALSE) == TRUE ) {
    # print out paths to files that don't exist
    cat( "Some files do not exist:",
         paste0("\t",paths[do_files_exist == FALSE]), sep="\n")
    stop()
  }

  all_nexus <- sapply(paths, isNexusFile)
  if ( all(all_nexus == TRUE) ) {
    trees <- lapply(paths, readNexusTrees, ...)
  } else if ( all(all_nexus == FALSE) ) {
    n_paths  <- length(paths)
    trees    <- vector("list", n_paths)
    for(i in 1:n_paths) {
      cat("Reading trees in file: ", paths[i], "\n", sep="")
      trees[[i]] <- readTreeLogs(paths[i], tree_name = tree_name, burnin = burnin, verbose = verbose, ...)
    }
  } else {
    stop("Write a good error message here.")
  }

  recover()

  return(trees)

#  # read in tree(s) of type nexus or newick
#
#  if (format == "nexus") {
#    tree <- ape::read.nexus(file = path)
#  } else if (format == "newick") {
#    tree <- ape::read.tree(file = path)
#  }
#
#  # convert to type multiPhylo for consistency
#
#  if (class(tree) == "phylo") {
#    tree <- c(tree)
#  } else if (class(tree) != "multiPhylo") {
#    stop ("tree(s) not of type phylo or multiPhylo")
#  }
#
#  return(tree)
#
}
