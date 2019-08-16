# Non-exported utility functions for RevGadgets

isNexusFile <- function(file) readLines(file, n=1) == "#NEXUS"

parseTreeString <- function(string) {
  text <- sub("[^(]*", "", string)
  stats <- treeio:::read.stats_beast_internal( "", text )
  tree <- ape::read.tree(text = text)
  obj <- treeio:::BEAST("", text, stats, tree )
  return(obj)
}

readNexusTrees <- function(path, burnin, verbose, ...) {

  # read the lines
  lines <- readLines(path)

  # the line with a tree
  tree_strings <- lines[grep("tree ", lines, ignore.case=FALSE)]

  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin+1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin*length(tree_strings))
    tree_strings <- tree_strings[(discard+1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }

  # get the trees
  n_trees <- length(tree_strings)
  if ( verbose == TRUE ) {
    bar <- txtProgressBar(style=3, width=40)
  }
  trees <- vector("list", n_trees)
  for(i in 1:n_trees) {
    trees[[i]] <- parseTreeString( tree_strings[i] )
    if ( verbose == TRUE ) { setTxtProgressBar(bar, i / n_trees)  }
  }
  if ( verbose == TRUE ) {
    cat("\n")
  }

  # return the trees
  return(trees)

}

readTreeLogs <- function(path, tree_name, burnin, verbose, ...) {

  # read the samples
  samples <- read.table(path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

  # check that there is a column with the given name
  if (tree_name %in% colnames(samples) == FALSE ) {
    stop(paste0("No column named ", tree_name, " found."))
  }

  # get the tree strings
  tree_strings <- samples[,tree_name]

  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin+1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin*length(tree_strings))
    tree_strings <- tree_strings[(discard+1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }

  # get the trees
  n_trees <- length(tree_strings)
  if ( verbose == TRUE ) {
    bar <- txtProgressBar(style=3, width=40)
  }
  trees <- vector("list", n_trees)
  for(i in 1:n_trees) {
    trees[[i]] <- parseTreeString( tree_strings[i] )
    if ( verbose == TRUE ) { setTxtProgressBar(bar, i / n_trees)  }
  }
  if ( verbose == TRUE ) {
    cat("\n")
  }

  # return the trees
  return(trees)

}
