# Non-exported utility functions for RevGadgets

buildTranslateDictionary <- function(lines) {

  start_tree_block <- grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <- grep("end;", lines[start_tree_block:length(lines)], ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block : end_tree_block]

  # look for translate block

  start_translations <- grep("translate", tree_block, ignore.case = TRUE)
  end_translations <- grep(";", tree_block[start_translations:length(tree_block)])[1] + start_translations -  1
  translations <- tree_block[start_translations : end_translations]
  translations <- gsub("\\\t","", translations)
  translations <- gsub(",", "", translations)

  # grab only the numbers and taxa names
  translations <- translations[grep("[1-9]", translations)]
  translations_split <- lapply(translations, strsplit, split = " ")
  dictionary <- as.data.frame(matrix(unlist(translations_split),
                                     ncol = 2, byrow = TRUE),
                              stringsAsFactors = FALSE)
  colnames(dictionary) <- c("number","taxon")

  return(dictionary)
}

colFun <- function(n) {
  if (n == 1) {return("#005ac8")}
  if (n == 2) {return(c("#005ac8","#fa7850"))}
  if (n == 3) {return(c("#14d2dc","#005ac8","#fa7850"))}
  if (n == 4) {return(c("#14d2dc","#005ac8","#fa7850", "#aa0a3c"))}
  if (n == 5) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                        "#0ab45a"))}
  if (n == 6) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                        "#0ab45a","#006e82"))}
  if (n == 7) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                        "#0ab45a","#006e82", "#fa78fa"))}
  if (n == 8) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                        "#0ab45a","#006e82", "#fa78fa", "#8214a0"))}
  if (n == 9) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                        "#0ab45a","#006e82", "#fa78fa", "#8214a0",
                        "#fae6be"))}
  if (n == 10) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                         "#0ab45a","#006e82", "#fa78fa", "#8214a0",
                         "#fae6be", "#00a0fa"))}
  if (n == 11) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                         "#0ab45a","#006e82", "#fa78fa", "#8214a0",
                         "#fae6be", "#00a0fa","#f0f032"))}

  if (n == 12) {return(c("#14d2dc","#005ac8","#fa7850","#aa0a3c",
                         "#0ab45a","#006e82", "#fa78fa", "#8214a0",
                         "#fae6be", "#00a0fa","#f0f032", "#a0fa82"))}
  if (n >= 13 ) {stop("more than 12 colors is not supported")}
}

findTreeLines <- function(lines) {

  # pull out tree block only
  start_tree_block <- grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <- grep("end;", lines[start_tree_block:length(lines)], ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block : end_tree_block]

  # pull out trees

  # find all starting lines by searching for "tree"
  trees_start <- grep("tree ", tree_block, ignore.case=TRUE)
  # find all ending lines by searching for ";" except for the last line of the tree block
  semicols <- grep("\\;", tree_block)
  semicols <- semicols[ semicols >= trees_start[1] ]
  trees_end <- semicols[ 1 : (length(semicols) - 1) ]
  # if tree are each on one line, return tree strings, else concatenate multiple lines
  if (all(trees_start  == trees_end)) {
    tree_strings <- tree_block[grep("tree ", tree_block, ignore.case=TRUE)]
  } else {stop("RevGadgets currently doesn't support line breaks in trees in nexus files")}
  #  search  for  semicolon to signal end of line

  # return tree strings
  return(tree_strings)

}

isNexusFile <- function(file) readLines(file, n=1) == "#NEXUS"

parseTreeString <- function(string) {
  text <- sub("[^(]*", "", string)
  stats <- treeio:::read.stats_beast_internal( "", text )
  tree <- ape::read.tree(text = text)
  obj <- treeio:::BEAST("", text, stats, tree )
  return(obj)
}

# Right tail probability of the horseshoe
# Integrates the density function via grid
pRightTailHorseshoeGrid <- function(x, gamma=1, grid.size=5000) {
  quants <- seq(1e-10,1-1e-10,length.out=grid.size)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants)/2.0
  probs <- 1/length(quants) # we're using quantiles, each gamma is equally likely
  sigmas <- qcauchy(quants,0,gamma)
  sum(pnorm(x,0,sigmas,lower.tail=FALSE) * probs)
}

readNexusTrees <- function(path, burnin, verbose, ...) {

  # read the lines
  lines <- readLines(path)

  # the line with a tree
  tree_strings <- findTreeLines(lines)

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

  # translate using dictionary if translate block present in file
  if (length(grep("translate", lines, ignore.case = TRUE)) >= 1) {
    dictionary <- buildTranslateDictionary(lines = lines)
    for (i in 1:n_trees) {
      n_tips <- length(trees[[i]]@phylo$tip.label)
      for (j in 1:n_tips) {
        ind <- which(trees[[i]]@phylo$tip.label[j] == dictionary[, 1])
        trees[[i]]@phylo$tip.label[j] <- dictionary[ind, 2]
      }
    }
  }

  # return the trees
  return(trees)

}

readTreeLogs <- function(path, tree_name, burnin, verbose, ...) {

  # read the samples
  samples <- utils::read.table(path, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)

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

