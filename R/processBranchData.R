#' processBranchData
#'
#' @param tree (treedata object; no default) a phylogenetic tree in the
#' treedata format, or a list of lists of a single tree data object, such as the
#' output of readTrees().
#' @param dat (data.frame or list; no default) a data frame, or a list
#' (of length 1) of a data frame, with branch specific data, such as the output
#' of readTrace().
#' @param burnin (numeric; 0.25) fraction of the markov-chain to discard
#' @param parnames (character vector; c("avg_lambda", "avg_mu", "num_shifts"))
#' Names of parameters to process
#' @param summary (character; "median") function to summarize the continuous
#' parameter. Typically mean or median
#' @param net_div (logical; FALSE) Calculate net diversification?
#'
#' @return a treedata file with attached branch-specific data
#' @export
#'
#' @examples
#' \donttest{
#'
#' # download the example dataset to working directory
#' url_rates <-
#'   "https://revbayes.github.io/tutorials/intro/data/primates_BDS_rates.log"
#' dest_path_rates <- "primates_BDS_rates.log"
#' download.file(url_rates, dest_path_rates)
#'
#' url_tree <-
#'   "https://revbayes.github.io/tutorials/divrate/data/primates_tree.nex"
#' dest_path_tree <- "primates_tree.nex"
#' download.file(url_tree, dest_path_tree)
#'
#' # to run on your own data, change this to the path to your data file
#' treefile <- dest_path_tree
#' logfile <- dest_path_rates
#'
#' branch_data <- readTrace(logfile)
#' tree <- readTrees(paths = treefile)
#'
#' annotated_tree <- processBranchData(tree, branch_data, summary = "median")
#'
#' # you can plot this output
#' p <- plotTree(tree = annotated_tree,
#'               node_age_bars = FALSE,
#'               node_pp = FALSE,
#'               tip_labels = FALSE,
#'               color_branch_by = "avg_lambda",
#'               line_width = 0.8,
#'               branch_color = c("blue","green")) +
#'      ggplot2::theme(legend.position=c(.1, .9));p
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(dest_path_tree, dest_path_rates)
#' }
processBranchData <- function(tree,
                              dat,
                              burnin = 0.25,
                              parnames = c("avg_lambda",
                                           "avg_mu",
                                           "num_shifts"),
                              summary = "median",
                              net_div = FALSE) {
  if (class(dat) == "list") {
    dat <- dat[[1]]
  }
  if (class(tree) == "list") {
    tree <- tree[[1]][[1]]
  }

  if (!"data.frame" %in% class(dat))
    stop("dat must be a data.frame or a single list of a data.frame")
  if (!"treedata" %in% class(tree))
    stop("tree must be a treedata object or a list of
         lists of treedata objects")

  dat <- dat[floor(nrow(dat) * burnin):nrow(dat), ]
  tree_tbl <- tibble::as_tibble(tree)
  map <- matchNodes(tree@phylo)

  for (item in parnames) {
    parameter <-
      unname(unlist(lapply(dat[, grepl(item,
                                       colnames(dat))],
                           summary)))[map$Rev]
    tree_tbl[[item]] <- parameter
  }

  if (net_div) {
    if ("avg_lambda" %in% parnames & "avg_mu" %in% parnames) {
      lambdas <- as.matrix(dat[, grepl("avg_lambda", colnames(dat))])
      mus <- as.matrix(dat[, grepl("avg_mu", colnames(dat))])
      net_divs <- as.data.frame(lambdas - mus)
      tree_tbl[["net_div"]] <-
        unname(unlist(lapply(net_divs, summary)))[map$Rev]
    } else {
      stop(
        "You set net_div = TRUE. Cannot calculate net_div without
        'avg_lambda and avg_mu' in parnames"
      )
    }
  }
  tree2 <- tidytree::as.treedata(tree_tbl)

  return(list(list(tree2)))
}
