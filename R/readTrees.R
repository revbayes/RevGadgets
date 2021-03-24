#' Read trees
#'
#' Reads in a tree file containing one or multiple trees
#'
#' Reads in a tree file in either nexus or newick format, and containing a single tree
#' or multiple trees (as in the results of a Bayesian analysis). For reading in annotated
#' tree files of continuous character evolution, the parameter must be considered a node
#' parameter rather than branch parameter. Set isNodeParameter = TRUE in the extended
#' newick monitor (mnExtNewick)
#'
#' @param paths (vector of character strings; no default) File path(s) to tree(s).
#' @param tree_name (character string; default psi) Name of the tree variable.
#' @param burnin (single numeric value; default = 0.1) Fraction of generations to
#' discard (if value provided is between 0 and 1) or number of generations (if
#' value provided is greater than 1).
#' @param n_cores (integer; default 1) Number of cores for parallelizing.
#' @param verbose (logical; default true) Display a status bar?
#'
#' @return A list (across runs) of lists (across samples) of treedata objects.
#'
#' @examples
#'
#' \dontrun{
#'
#' # read in a single nexus file
#' file <- system.file("extdata",
#'     "sub_models/primates_cytb_GTR_MAP.tre", package="RevGadgets")
#' tree_single <- readTrees(paths = file)
#'
#' # read in a single newick string
#' file <- system.file("extdata", "bds/primates.tre", package="RevGadgets")
#' tree_new <- readTrees(path = file)
#'
#' # read in a tree trace (may take a few seconds)
#' file <- system.file("extdata", "sub_models/primates_cytb_GTR.trees", package="RevGadgets")
#' tree_multi <- readTrees(path = file)
#'
#' }
#' @export

readTrees <- function(paths, tree_name =  "psi", burnin = 0, n_cores = 1L, verbose = TRUE) {
  # enforce argument matching
  if (is.character(tree_name) == FALSE) stop("tree_name should be a single character")
  if (is.numeric(burnin) == FALSE) stop("burnin should be a number")
  if (is.logical(verbose) == FALSE) stop("verbose should be TRUE or FALSE")

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

  n_paths  <- length(paths)
  trees    <- vector("list", n_paths)
  for (i in 1:length(paths)){
    nexus <- .isNexusFile(paths[i])
    newick_single <- .isSingleNewick(paths[i])
    if (!newick_single & !nexus){
      newick_trace <- T
    } else {newick_trace <- F}

    if (nexus & !newick_single & !newick_trace){

      trees[[i]] <- .readNexusTrees(path = paths[i], burnin = burnin, verbose = verbose)

    } else if (!nexus & newick_single & !newick_trace) {

      tree_string <- readLines(paths[i], n=1)
      trees[[i]] <- list(.parseTreeString(tree_string))

    } else if (!nexus & !newick_single & newick_trace) {

      trees[[i]] <- .readTreeLogs(path = paths[i], tree_name = tree_name, burnin = burnin,
                                  verbose = verbose)

    } else {stop("tree file format unrecognized")}

  }
  return(trees)
}

