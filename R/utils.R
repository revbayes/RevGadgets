# Non-exported utility functions for RevGadgets
.add_epoch_times <- function( p, max_age, dy_bars, dy_text ) {

  max_x = max(p$data$x)
  max_y = max(p$data$y)

  epoch_names <- rev(c("Holocene", "Pleistocene", "Pliocene",
                   "Miocene", "Oligocene", "Eocene", "Paleocene",
                   "Upper\nCretaceous", "Lower\nCretaceous"))
  epoch_ages <- rev(c(0, 0.117, 2.588, 5.332, 23.03,
                      33.9, 55.8, 65.5, 99.6))
  period_names <- rev(c("Quaternary", "Neogene", "Paleogene",
                    "Cretaceous", "Jurassic", "Triassic",
                    "Permian", "Carboniferous", "Devonian",
                    "Silurian", "Ordovician", "Cambrian"))

  period_ages <- rev(c(0, 2.588, 23.03, 65.5, 145.5,
                       199.6, 251, 299, 359.2, 416,
                       443.7, 488.3))
  if (max_age > 140) {
    x_pos <- max_x - c(max_age, period_ages)
  } else x_pos <- max_x - c(max_age, epoch_ages)

  y_pos = rep(max_y, length(x_pos))
  x_pos_mid = ( x_pos[1:(length(x_pos)-1)] + x_pos[2:length(x_pos)] ) / 2

  for (k in 2:(length(x_pos))) {
    box_col = "gray92"
    if (k %% 2 == 0) box_col = "white"
    box = ggplot2::geom_rect( xmin=x_pos[k-1], xmax=x_pos[k], ymin=dy_bars, ymax=y_pos[k], fill=box_col )
    p = gginnards::append_layers(p, box, position = "bottom")
  }
  if (max_age > 140) {
    for (k in 1:length(period_names)) {
      p <- p + ggplot2::annotate( geom="text", label=period_names[k], angle = 90,
                                  x=x_pos_mid[k], y=dy_text, hjust=0, size=3.25)
    }
  } else {
    for (k in 1:length(epoch_names)) {
      p <- p + ggplot2::annotate( geom="text", label=epoch_names[k], angle = 90,
                                  x=x_pos_mid[k], y=dy_text, hjust=0, size=3.25)
    }
  }
  return(p)
}

.build_state_probs <- function(t, state_labels, include_start_states, p_threshold = 0.01) {

  n_states = length(state_labels)
  n_tips = length(attributes(t)$phylo$tip.label)
  n_node = 2 * n_tips - 1

  dat = list()

  if (include_start_states) {
    state_tags = c("start","end")
  } else {
    state_tags = c("anc")
  }

  for (s in state_tags) {
    dat[[s]] = data.frame( matrix(0, nrow=n_node, ncol=n_states) )
    #dat[[s]] = cbind(node=1:n_node, dat[[s]])

    for (i in 1:3)
    {
      m = paste(s,"_state_",i,sep="")
      pp_str = paste(m,"_pp",sep="")
      n_tmp = as.numeric(as.vector(attributes(t)$data$node)) # node index
      x_tmp = as.vector(attributes(t)$data[[m]])
      pp_tmp = as.numeric(as.vector(attributes(t)$data[[pp_str]]))

      for (j in 1:length(x_tmp))
      {
        if (!is.na(x_tmp[j])) {

          if (pp_tmp[j] > p_threshold) {
            k = which(x_tmp[j]==state_labels)
            dat[[s]][n_tmp[j], k] = pp_tmp[j]
          }
        }
      }
    }

    # format column names
    colnames(dat[[s]])=as.vector(unlist(state_labels))

    # add probs for >3rd state under ... label
    rem_prob = c()
    for (i in 1:nrow(dat[[s]])) {
      rem_prob[i] = 1
      for (j in 1:length(dat[[s]][i,])) {
        rem_prob[i] = rem_prob[i] - dat[[s]][i,j]
      }
    }
    dat[[s]]$`...` = rem_prob
    dat[[s]]$node = 1:n_node
    #print(dat[[s]][250:260,])
  }

  return(dat)
}

.buildTranslateDictionary <- function(lines) {

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

.colFun <- function(n) {
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

.collect_probable_states <- function(p, p_threshold=0.005) {
  labels = c("end_state", "start_state")
  index = c(1,2,3)

  codes = c()
  labels_pp = c()
  for (l in labels) {
    for (i in index) {
      label_index = paste(l,"_",i,sep="")
      label_index_pp = paste(l,"_",i,"_pp",sep="")
      index_threshold = p$data[[ label_index_pp ]] > p_threshold
      codes = c(codes, unique( p$data[[label_index]][ index_threshold ] ))
    }
  }
  codes = unique(codes)
  codes = c(codes, "...")
  return(codes)
}

.convertAndRound <- function(L) {
  #sometimes there will be NAs before forcing to convert - got to remove nas before doing this test!
  k <- L[!is.na(L)]
  if (any(is.na(as.numeric(k))) == FALSE) { # if integer or numeric
    if (sum(as.numeric(L) %% 1, na.rm = T) == 0) { # if integer
      labs <- L
    } else { # if numeric
      labs <- sprintf("%.3f",as.numeric(L)) # round nicely
      labs[labs == "1.000"] <- "1"
    }
  } else { # if character
    labs <- L
  }
  return(labs)
}

.findTreeLines <- function(lines) {

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

.isColor <- function(var) {
  if (is.null(var)) {
    return(FALSE) } else {
  t <- try(col2rgb(var), silent = TRUE)
  if (class(t) == "try-error") {return(FALSE)}
  else return(TRUE)
    }
}

.isNexusFile <- function(file) readLines(file, n=1) == "#NEXUS"

.parseTreeString <- function(string) {
  text <- sub("[^(]*", "", string)
  # stats <- treeio:::read.stats_beast_internal( "", text )
  stats <- .read.stats_revbayes_internal( "", text )
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

.read.stats_revbayes_internal <- function(beast, tree) {

  phylo <- ape::read.tree(text = tree) # read the tree
  tree2 <- treeio:::add_pseudo_nodelabel(phylo) # add nodelabels (if there aren't already any)

  ## node name corresponding to stats
  nn <- unlist(strsplit(unlist(strsplit(tree2, split=",")), "\\)"))
  nn <- gsub("\\(*", "", nn)
  nn <- gsub("[:;].*", "", nn)
  nn <- gsub(" ", "", nn)
  nn <- gsub("'", "", nn)
  nn <- gsub('"', "", nn)

  # nn <- strsplit(tree2, split=",") %>% unlist %>%
  #   strsplit(., split="\\)") %>% unlist %>%
  #   gsub("\\(*", "", .) %>%
  #   gsub("[:;].*", "", .) %>%
  #   gsub(" ", "", .) %>%
  #   gsub("'", "", .) %>%
  #   gsub('"', "", .)
  # get the label of each node (internal and external) in the order
  # they appear in the newick string

  phylo <- ape::read.tree(text = tree2)
  root <- tidytree:::rootnode(phylo)
  nnode <- phylo$Nnode

  tree_label <- c(phylo$tip.label, phylo$node.label)
  ii <- match(nn, tree_label)

  if ( any(grepl("TRANSLATE", beast, ignore.case = TRUE)) ) {
    label2 <- c(phylo$tip.label, root:treeio:::getNodeNum(phylo))
  } else {
    label2 <- as.character(1:treeio:::getNodeNum(phylo))
  }
  node <- label2[match(nn, tree_label)]

  ## stats <- unlist(strsplit(tree, "\\["))[-1]
  ## stats <- sub(":.+$", "", stats

  ## BEAST1 edge stat fix
  tree <- gsub("\\]:\\[&(.+?\\])", ",\\1:", tree)
  tree <- gsub(":(\\[.+?\\])", "\\1:", tree)

  if (grepl("\\]:[0-9\\.eE+\\-]*\\[", tree) || grepl("\\]\\[", tree)) {
    ## MrBayes output
    # stats <- strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\[") %>% unlist
    stats <- unlist(strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\["))
    lstats <- lapply(stats, function(x) {
      unlist(strsplit(x, split="\\][,\\)]"))
    })

    for (i in seq_along(stats)) {
      n <- length(lstats[[i]])
      if (i == length(stats)) {
        stats[i] <- lstats[[i]][n]
      } else {
        stats[i] <- paste0(lstats[[i]][n],
                           sub("&", ",", lstats[[i+1]][1])
        )
      }
    }
    stats <- gsub("\\]\\[&", ",", stats)
  } else {
    ## BEAST output
    stats <- unlist(strsplit(tree, ":"))
  }

  names(stats) <- node

  stats <- stats[grep("\\[", stats)]
  stats <- sub("[^\\[]*\\[", "", stats)

  stats <- sub("^&", "", stats)
  stats <- sub("];*$", "", stats)
  stats <- gsub("\"", "", stats)

  #this is what is breaking readTrees for the OU output
  stats2 <- lapply(seq_along(stats), function(i) {

    x <- stats[[i]]
    y <- unlist(strsplit(x, ","))
    sidx <- grep("=\\{", y)
    eidx <- grep("\\}$", y)

    flag <- FALSE
    if (length(sidx) > 0) {
      flag <- TRUE
      # SETS <- lapply(seq_along(sidx), function(k) {
      #   p <- y[sidx[k]:eidx[k]]
      #   gsub(".*=\\{", "", p) %>% gsub("\\}$", "", .)
      # })
      SETS <- lapply(seq_along(sidx), function(k) {
        p <- y[sidx[k]:eidx[k]]
        gsub("\\}$", "",gsub(".*=\\{", "", p))
      })
      names(SETS) <- gsub("=.*", "", y[sidx])

      kk <- unlist(lapply(seq_along(sidx), function(k) {
        sidx[k]:eidx[k]
      }))
      y <- y[-kk]
    }

    if (length(y) == 0)
      return(SETS)

    name <- gsub("=.*", "", y)
    val <- gsub(".*=", "", y)
    val <- gsub("^\\{", "", val)
    val <- gsub("\\}$", "", val)

    # %>%
    #   gsub("^\\{", "", .) %>%
    #   gsub("\\}$", "", .)

    if (flag) {
      nn <- c(name, names(SETS))
    } else {
      nn <- name
    }

    res <- rep(NA, length(nn))
    names(res) <- nn

    for (i in seq_along(name)) {
      # res[i] <- if(treeio:::is.numeric(val[i])) as.numeric(val[i]) else val[i]
      res[i] <- if(is.numeric(val[i])) as.numeric(val[i]) else val[i]
    }
    if (flag) {
      j <- i
      for (i in seq_along(SETS)) {
        if(is.numeric(SETS[[i]])) {
          res[i+j] <- list(as.numeric(SETS[[i]]))
        } else {
          res[i+j] <- SETS[i]
        }
      }
    }

    return(res)
  })

  nn <- sort(unique(unlist(lapply(stats2, names))))

  # nn <- lapply(stats2, names) %>% unlist %>%
  #   unique %>% sort

  stats2 <- lapply(stats2, function(x) {
    y <- x[nn]
    names(y) <- nn
    y[vapply(y, is.null, logical(1))] <- NA
    y
  })

  stats3 <- do.call(rbind, stats2)
  stats3 <- tibble::as_tibble(stats3)

  ## no need to extract sd from prob+-sd
  ## as the sd is stored in prob_stddev
  ##
  ## "prob_stddev"   "prob(percent)" "prob+-sd"
  ##
  ##
  ##
  ## idx <- grep("\\+-", colnames(stats3))
  ## if (length(idx)) {
  ##     for (i in idx) {
  ##         stats3[,i] <- as.numeric(gsub("\\d+\\+-", "", stats3[,i]))
  ##     }
  ## }

  cn <- gsub("(\\d+)%", "0.\\1", colnames(stats3))
  cn <- gsub("\\(([^\\)]+)\\)", "_\\1", cn)
  ## cn <- gsub("\\+-", "_", cn)

  colnames(stats3) <- cn
  stats3$node <- names(stats)

  i <- vapply(stats3,
              function(x) max(vapply(x, length, numeric(1))),
              numeric(1))

  for (j in which(i==1)) {
    stats3[,j] <- unlist(stats3[,j])
  }
  stats3$node <- as.integer(stats3$node)
  return(stats3)
}

.readNexusTrees <- function(path, burnin, verbose, ...) {
  # read the lines
  lines <- readLines(path)

  # the line with a tree
  tree_strings <- .findTreeLines(lines)

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
    trees[[i]] <- .parseTreeString( tree_strings[i] )
    if ( verbose == TRUE ) { setTxtProgressBar(bar, i / n_trees)  }
  }
  if ( verbose == TRUE ) {
    cat("\n")
  }

  # translate using dictionary if translate block present in file
  if (length(grep("translate", lines, ignore.case = TRUE)) >= 1) {
    dictionary <- .buildTranslateDictionary(lines = lines)
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

.readTreeLogs <- function(path, tree_name, burnin, verbose, ...) {

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
    trees[[i]] <- .parseTreeString( tree_strings[i] )
    if ( verbose == TRUE ) { setTxtProgressBar(bar, i / n_trees)  }
  }
  if ( verbose == TRUE ) {
    cat("\n")
  }

  # return the trees
  return(trees)

}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}
