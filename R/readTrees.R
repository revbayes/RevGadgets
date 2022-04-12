#' Read trees
#'
#' Reads in a tree file containing one or multiple trees
#'
#' Reads in a tree file in either nexus or newick format, and containing a
#' single tree or multiple trees (as in the results of a Bayesian analysis).
#' For reading in annotated tree files of continuous character evolution,
#' the parameter must be considered a node parameter rather than branch
#' parameter. Set isNodeParameter = TRUE in the extended newick monitor
#' (mnExtNewick)
#'
#' @param paths (vector of character strings; no default) File path(s) to
#' tree(s).
#' @param tree_name (character string; default psi) Name of the tree variable.
#' @param burnin (single numeric value; default = 0.1) Fraction of generations
#' to discard (if value provided is between 0 and 1) or number of generations
#' (if value provided is greater than 1).
#' @param n_cores (integer; default 1) Number of cores for parallelizing.
#' @param verbose (logical; default true) Display a status bar?
#'
#' @return A list (across runs) of lists (across samples) of treedata objects.
#'
#' @examples
#'
#' \donttest{
#' # read in a single nexus file
#'
#' # download the example dataset to working directory
#' url_nex <-
#'  "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR_MAP.tre"
#' dest_path_nex <- "primates_cytb_GTR_MAP.tre"
#' download.file(url_nex, dest_path_nex)
#'
#' # to run on your own data, change this to the path to your data file
#' file <- dest_path_nex
#' tree_single_old <- readTrees(paths = file)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_nex)
#'
#' # read in a single newick string
#'
#' # download the example dataset to working directory
#' url_new <-
#'  "https://revbayes.github.io/tutorials/intro/data/primates.tre"
#' dest_path_new <- "primates.tre"
#' download.file(url_new, dest_path_new)
#'
#' # to run on your own data, change this to the path to your data file
#' file_new <- dest_path_new
#' tree_new <- readTrees(paths = file_new)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_new)
#'
#'
#' # read in a tree trace (may take a few seconds)
#'
#' # download the example dataset to working directory
#' url_multi <-
#'  "https://revbayes.github.io/tutorials/intro/data/primates_cytb_GTR.trees"
#' dest_path_multi <- "primates_cytb_GTR.trees"
#' download.file(url_multi, dest_path_multi)
#'
#' # to run on your own data, change this to the path to your data file
#' file_multi <- dest_path_multi
#' tree_multi <- readTrees(paths = file_multi)
#'
#' # remove file
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_multi)
#' }
#'
#' @export

readTrees <-
  function(paths,
           tree_name =  "psi",
           burnin = 0,
           n_cores = 1L,
           verbose = TRUE) {
    # enforce argument matching
    if (is.character(tree_name) == FALSE)
      stop("tree_name should be a single character")
    if (is.numeric(burnin) == FALSE)
      stop("burnin should be a number")
    if (is.logical(verbose) == FALSE)
      stop("verbose should be TRUE or FALSE")

    character_paths_are_strings <- is.character(paths)
    if (any(character_paths_are_strings == FALSE) == TRUE) {
      # print out the ones that are not character strings
      stop(
        paste0("Some paths are not character strings:",
          paste0("\t", paths[character_paths_are_strings == FALSE]),
          sep = "\n")
        )
    }

    do_files_exist <- file.exists(paths)
    if (any(do_files_exist == FALSE) == TRUE) {
      # print out paths to files that don't exist
      stop(
        paste0("Some files do not exist:",
          paste0("\t", paths[do_files_exist == FALSE]), sep = "\n")
      )
    }

    n_paths  <- length(paths)
    trees    <- vector("list", n_paths)
    for (i in seq_len(length(paths))) {
      nexus <- .isNexusFile(paths[i])
      newick_single <- .isSingleNewick(paths[i])
      if (!newick_single & !nexus) {
        newick_trace <- TRUE
      } else {
        newick_trace <- FALSE
      }

      if (nexus & !newick_single & !newick_trace) {
        trees[[i]] <-
          .readNexusTrees(path = paths[i],
                          burnin = burnin,
                          verbose = verbose)

      } else if (!nexus & newick_single & !newick_trace) {
        tree_string <- readLines(paths[i], n = 1)
        trees[[i]] <- list(.parseTreeString(tree_string))

      } else if (!nexus & !newick_single & newick_trace) {
        trees[[i]] <-
          .readTreeLogs(
            path = paths[i],
            tree_name = tree_name,
            burnin = burnin,
            verbose = verbose
          )

      } else {
        stop("tree file format unrecognized")
      }

      # add index if missing (for trees not output by RevBayes)
      for (j in seq_len(length(trees[[i]]))) {
        if (!"index" %in% colnames(trees[[i]][[j]]@data)) {
          t <- trees[[i]][[j]]
          if (!"node" %in% colnames(t@data) ||
              length(t@data$node) == 0) {
            trees[[i]][[j]]@data <-
              dplyr::tibble(node = trees[[i]][[j]]@phylo$edge[, 2])
          }
          node_matches <- dplyr::as_tibble(matchNodes(t@phylo))
          colnames(node_matches) <- c("node", "index")
          class(node_matches$node) <-
            class(node_matches$index) <- "character"
          trees[[i]][[j]]@data <-
            dplyr::left_join(node_matches, t@data)
        }
      }
    }
    return(trees)
  }
