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

# set custom state labels
.assign_state_labels <- function(t, state_labels, include_start_states, labels_as_numbers, n_states=3) {
  #recover()
  # what is the ancestral state name tag?
  if (include_start_states) {
    state_pos_str_base = c("start_state_", "end_state_")
  } else {
    state_pos_str_base = c("anc_state_")
  }

  # send error if state_labels are provided without names
  if (!is.null(state_labels) && is.null(names(state_labels))) {
    error("names(state_labels) must identify all unlabeled state names in attributes(t)$data")
  }

  #make matrix of all anc state values
  col_num <- grep(state_pos_str_base[1], colnames(t@data))
  if (length(state_pos_str_base) > 1){
    col_num2 <- grep(state_pos_str_base[2], colnames(t@data))
    col_num <- c(col_num, col_num2)
  }
  pps <- grep("_pp", colnames(t@data))
  columns <- col_num[!col_num %in% pps]

  all_anc_states <- unique(c(as.matrix(t@data[, columns])))

  # send error if state labels are provided but there are any
  # states without a corresponding state label
  if (!is.null(state_labels) && any(all_anc_states %in% c("NA", names(state_labels)) == FALSE)) {
    stop(paste0("names(state_labels): ",
                paste0(names(state_labels), collapse = ", "),
                " do not match data in tree file: ",
                paste0(sort(all_anc_states[all_anc_states != "NA"]), collapse = ", ")))
  }

  # generate state labels if none provided and not a chromosome analysis
  if ( is.null(state_labels) == TRUE & labels_as_numbers == FALSE) {
    warning("State labels not provided by user. Will be generated automatically.")
    states <- unique(unlist(attributes(t)$data[grepl(paste0("state_","[0-9]$"),names(attributes(t)$data))]))
    states <- states[!states == "NA"]
    states <- states[order(states)]
    state_labels <- list()
    for(i in 1:length(states) ) {
      state_labels[as.character(states[i])] = LETTERS[i]
    }
    state_labels["other"] <- "other"
  }

  # for chromosome analyses, just keep the names as is (numbers of chromos)
  if ( is.null(state_labels) == TRUE & labels_as_numbers == TRUE) {
    state_labels <- unique(unlist(attributes(t)$data[grepl(paste0("state_","[0-9]$"),names(attributes(t)$data))]))
    state_labels <- state_labels[-which(state_labels == "NA")]
    names(state_labels) <- state_labels
  }

  # create list of ancestral state name tags
  state_pos_str_to_update = c(sapply(1:n_states, function(x) { paste(state_pos_str_base,x,sep="")}))

  # overwrite state labels
  for (m in state_pos_str_to_update)
  {
    # get the states
    x_state = attributes(t)$data[[m]]
    x_state = as.vector(x_state)
    x_state_valid = which( x_state != "NA" )
    x_state_invalid = which( x_state == "NA" )
    x_state_tmp = unlist(sapply(x_state, function(z) { state_labels[ names(state_labels)==z ] }))
    x_state[x_state_valid] = x_state_tmp
    x_state[x_state_invalid] = NA
    attributes(t)$data[[m]] = x_state
  }

  unique(c(as.matrix(t@data[, columns])))
  # Just add the USED state_labels here
  used_state_labels <-  na.omit(unique(c(as.matrix(t@data[, columns]))))
  attributes(t)$state_labels <- as.character(used_state_labels)

  return(t)
}

.build_state_probs <- function(t, state_labels, include_start_states, p_threshold = 0) {
  n_states = length(state_labels)
  n_tips = length(attributes(t)$phylo$tip.label)
  n_node = 2 * n_tips - 1

  dat = list()

  if (include_start_states == TRUE) {
    state_tags = c("end","start")
  } else if (include_start_states == FALSE & "anc_state_1" %in% colnames(t@data)) {
    state_tags = c("anc")
  } else if (include_start_states == FALSE & "end_state_1" %in% colnames(t@data)) {
    state_tags = c("end")
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

    # add probs for >3rd state under "other" label
    rem_prob = c()
    for (i in 1:nrow(dat[[s]])) {
      rem_prob[i] = 1
      for (j in 1:length(dat[[s]][i,])) {
        rem_prob[i] = rem_prob[i] - dat[[s]][i,j]
      }
    }
    dat[[s]]$"other" = rem_prob
    dat[[s]]$node = 1:n_node
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

.collect_probable_states <- function(p, p_threshold = 0.005) {
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
  codes = c(codes, "other")
  return(codes)
}

.computeMeanInterval <- function(item, rates, probs){
  interval_times <- unlist(rates[["speciation time"]][1,grepl("interval_times", names(rates$`speciation time`))])
  interval_times <- sort(interval_times) # For some reason these are ordered differently than rate vectors

  rate <- rates[[item]]
  rate <- rate[, grep("[0-9]", colnames(rate))]

  mean_rate <- colMeans(rate)
  quantiles <- apply(rate, 2,
                     quantile,
                     probs = probs)

  df <- dplyr::tibble(.rows = length(mean_rate))
  df["mean"] <- mean_rate
  df["lower"] <- quantiles[1,]
  df["upper"] <- quantiles[2,]
  df$time <- interval_times
  df$item <- item

  return(df)
}

.convertAndRound <- function(L) {
  #sometimes there will be NAs before forcing to convert - got to remove nas before doing this test!
  k <- L[!is.na(L)]
  if (any(is.na(as.numeric(k))) == FALSE) { # if integer or numeric
    if (sum(as.numeric(L) %% 1, na.rm = T) == 0) { # if integer
      labs <- L
      labs[labs == "1.000000"] <- "1" # catch case of all posterios of 1
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

.geom_subview_revgadgets <- function (mapping = NULL, data = NULL, width = 0.1, height = 0.1,
                                      x = NULL, y = NULL, subview = NULL) {
  if (is.null(data)) {
    data <- dplyr::tibble(x = x, y = y)
  } else if (!inherits(data, "tbl")) {
    data <- dplyr::as_tibble(data)
  }
  if (is.null(mapping)) {
    mapping <- ggplot2::aes_(x = ~x, y = ~y)
  }

  mapping <- as.list(mapping)

  if (is.null(mapping$x)) {
    stop("x aesthetic mapping should be provided")
  }

  if (is.null(mapping$y)) {
    stop("y aesthetic mapping should be provided")
  }

  if (is.null(mapping$subview) && is.null(subview)) {
    stop("subview must be provided")
  }

  if (is.null(mapping$subview)) {
    if (!inherits(subview, "list")) {
      subview <- list(subview)
    }
    data$subview <- subview
  } else {
    sv_var <- rvcheck::get_aes_var(mapping, "subview")
    data$subview <- data[[sv_var]]
  }

  xvar <- rvcheck::get_aes_var(mapping, "x")
  yvar <- rvcheck::get_aes_var(mapping, "y")

  if (is.null(mapping$width)) {
    data$width <- width
  } else {
    width_var <- rvcheck::get_aes_var(mapping, "width")
    data$width <- data[[width_var]]
  }

  if (is.null(mapping$height)) {
    data$height <- height
  } else {
    height_var <- rvcheck::get_aes_var(mapping, "height")
    data$height <- data[[height_var]]
  }

  data$xmin <- data[[xvar]] - data$width/(2*max(data[[yvar]]))
  data$xmax <- data[[xvar]] + data$width/(2*max(data[[yvar]]))
  data$ymin <- data[[yvar]] - data$width/(2*max(data[[xvar]]))
  data$ymax <- data[[yvar]] + data$width/(2*max(data[[xvar]]))

  # save pies as images and plot as raster grobs - this plots centered but maybe not
  # a good long-term solution

   results <- list()
     for (i in 1:nrow(data)){
       ggplot2::ggsave(".temp.png",
                       plot = data$subview[[i]],
                       bg = "transparent",
                       width = 3, height = 3,
                       units = "cm", dpi = 200)
       pie <- png::readPNG(".temp.png")
       g <- grid::rasterGrob(pie, interpolate=TRUE)
       results[[i]] <- ggplot2::annotation_custom(ggplotify::as.grob(g),
                                  xmin = data$xmin[i],  xmax = data$xmax[i],
                                  ymin = data$ymin[i], ymax = data$ymax[i])
     }
   file.remove(".temp.png")
   return(results)

   # old way of plotting pies - won't plot centered on nodes
  # lapply(1:nrow(data), function(i) {
  #   ggplot2::annotation_custom(ggplotify::as.grob(data$subview[[i]]), xmin = data$xmin[i],
  #                              xmax = data$xmax[i], ymin = data$ymin[i], ymax = data$ymax[i])
  #})
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getParent <- function(tr, node) {
  if ( node == ggtree:::getRoot(tr) )
    return(0)
  edge <- tr[["edge"]]
  parent <- edge[,1]
  child <- edge[,2]
  res <- parent[child == node]
  if (length(res) == 0) {
    stop("cannot found parent node...")
  }
  if (length(res) > 1) {
    stop("multiple parent found...")
  }
  return(res)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getXcoord <- function(tr) {
  edge <- tr$edge
  parent <- edge[,1]
  child <- edge[,2]
  root <- ggtree:::getRoot(tr)

  len <- tr$edge.length

  N <- ggtree:::getNodeNum(tr)
  x <- numeric(N)
  x <- .getXcoord2(x, root, parent, child, len)
  return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getXcoord2 <- function(x, root, parent, child, len, start=0, rev=FALSE) {
  x[root] <- start
  x[-root] <- NA  ## only root is set to start, by default 0

  currentNode <- root
  direction <- 1
  if (rev == TRUE) {
    direction <- -1
  }
  while(anyNA(x)) {
    idx <- which(parent %in% currentNode)
    newNode <- child[idx]
    x[newNode] <- x[parent[idx]]+len[idx] * direction
    currentNode <- newNode
  }

  return(x)
}

# modified from https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getYcoord <- function(tr, step=1) {
  Ntip <- length(tr[["tip.label"]])
  N <- ggtree:::getNodeNum(tr)

  edge <- tr[["edge"]]
  parent <- edge[,1]
  child <- edge[,2]

  cl <- split(child, parent)
  child_list <- list()
  child_list[as.numeric(names(cl))] <- cl

  y <- numeric(N)
  tip.idx <- child[child <= Ntip]
  y[tip.idx] <- 1:Ntip * step
  y[-tip.idx] <- NA

  currentNode <- 1:Ntip
  while(anyNA(y)) {
    pNode <- unique(parent[child %in% currentNode])
    ## piping of magrittr is slower than nested function call.
    ## pipeR is fastest, may consider to use pipeR
    ##
    ## child %in% currentNode %>% which %>% parent[.] %>% unique
    ## idx <- sapply(pNode, function(i) all(child[parent == i] %in% currentNode))
    idx <- sapply(pNode, function(i) all(child_list[[i]] %in% currentNode))
    newNode <- pNode[idx]

    y[newNode] <- sapply(newNode, function(i) {
      mean(y[child_list[[i]]], na.rm=TRUE)
      ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
    })

    currentNode <- c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
    ## currentNode <- c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
    ## parent %in% newNode %>% child[.] %>%
    ##     `%in%`(currentNode, .) %>% `!` %>%
    ##         currentNode[.] %>% c(., newNode)
  }

  return(y)
}

.inset.revgadgets <- function (tree_view, insets, width = 0.1, height = 0.1, hjust = 0,
                               vjust = 0, x = "node", pos = 0.5) {
  df <- tree_view$data[as.numeric(names(insets)), ]

  # position subviews based on tree part
  x <- match.arg(x, c("node", "branch", "edge", "parent_shoulder"))
  if (x == "node") {
    xx <- df$x
  } else if (x == "parent_shoulder") {
    xx <- df$x[ match(df$parent, df$node) ]
  } else {
    xx <- df$branch
  }
  yy <- df$y
  xx <- xx - hjust # x-coordinates for nodes
  yy <- yy - vjust # y-coordinates for nodes

  if (length(width) == 1) { width <- rep(width, length(insets)) }
  if (length(height) == 1) { height <- rep(height, length(insets)) }

    # old way
  tree_view <- tree_view +
    .geom_subview_revgadgets(subview = insets, width = width,
                             height = height, x = xx, y = yy)
  # return treeview with subviews
  return(tree_view)
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

.makeNodeNames <- function(tree) {
  pr <- ape::prop.part(tree)
  labels <- attributes(pr)$labels
  names(labels) <- 1:length(labels)
  nodes <- lapply(pr[1:length(pr)], dplyr::recode, !!!labels)
  nodes <- append(attributes(pr)$labels, nodes)

  node_names <- numeric()
  node_names_op <- numeric()
  for (i in 1:length(nodes)) {
    node_names[i] <- paste(as.numeric(sort(tree$tip.label) %in% nodes[[i]]),
                           sep = "", collapse = "")
    node_names_op[i] <- paste(as.numeric(!sort(tree$tip.label) %in% nodes[[i]]),
                              sep = "", collapse = "")
  }
  return(data.frame(node_names = node_names,
                    node_names_op = node_names_op))
}

.makePlotData <- function(rates, probs){
  rates <- .removeNull(rates)
  res <- lapply(names(rates), function(e) .computeMeanInterval(e, rates = rates, probs = probs))
  plotdata <- do.call(rbind, res)
  plotdata$item <- factor(plotdata$item,
                          levels = c("speciation rate", "extinction rate", "speciation time", "extinction time",
                                     "net-diversification rate", "relative-extinction rate"))
  return(plotdata)
}

.makeStates <- function(label_fn, color_fn) {

  # generate colors for ranges
  range_color_list <- read.csv(color_fn, header=T, sep=",", colClasses="character")

  # get area names
  area_names <- unlist(sapply(range_color_list$range, function(y) { if (nchar(y)==1) { return(y) } }))

  # get state labels
  state_descriptions <- read.csv(label_fn, header=T, sep=",", colClasses="character")

  # map presence-absence ranges to area names
  range_labels <- sapply(state_descriptions$range[2:nrow(state_descriptions)],
                         function(x) {
                           present = as.vector(gregexpr(pattern="1", x)[[1]])
                           paste( area_names[present], collapse="")
                         })

  # map labels to colors
  range_colors <- range_color_list$color[ match(range_labels, range_color_list$range) ]

  # generate state/color labels
  idx <- 1
  st_lbl <- list()
  st_colors <- c()
  for (j in 1:(nrow(state_descriptions)-1)) {
    st_lbl[[ as.character(j) ]] <- range_labels[j]
    st_colors[j] <- range_colors[j]
  }
  st_colors[ length(st_colors)+1 ] <- "lightgray"
  st_lbl[["other"]] <- "other"

  return( list(state_labels = st_lbl, state_color = st_colors) )
}

# Fast data.frame constructor and indexing
# No checking, recycling etc. unless asked for
# Stolen from ggplot2
new_data_frame <- function(x = list(), n = NULL) {
  if (length(x) != 0 && is.null(names(x))) {
    abort("Elements must be named")
  }
  lengths <- vapply(x, length, integer(1))
  if (is.null(n)) {
    n <- if (length(x) == 0 || min(lengths) == 0) 0 else max(lengths)
  }
  for (i in seq_along(x)) {
    if (lengths[i] == n) next
    if (lengths[i] != 1) {
      abort("Elements must equal the number of rows or 1")
    }
    x[[i]] <- rep(x[[i]], n)
  }

  class(x) <- "data.frame"

  attr(x, "row.names") <- .set_row_names(n)
  x
}

.parseTreeString <- function(string) {
  text <- sub("[^(]*", "", string)
  # stats <- treeio:::read.stats_beast_internal( "", text )
  stats <- .read.stats_revbayes_internal( "", text )
  tree <- ape::read.tree(text = text)
  obj <- treeio:::BEAST("", text, stats, tree )
  return(obj)
}

# Right tail probability of the horseshoe: integrates the density function via grid
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

### Functions required by densiTreeWithBranchData
# attribute colors to a vector based the value in a range
color_gradient <- function(x, intervals = seq(0,11,0.1), colors = c("red","yellow","green"), bias = 1) {
  colfun <- grDevices::colorRampPalette(colors, bias = bias)
  return(  colfun(length(intervals)) [ findInterval(x, intervals, all.inside = TRUE) ] )
}

# function to sort a treedata
sort_tips <- function(x) {
  x <- reorder_treedata(x)
  nTip <- as.integer(length(x@phylo$tip.label))
  e2 <- x@phylo$edge[, 2]
  x@data <- x@data[c(e2[e2 <= nTip], (nTip+1):(nTip + x@phylo$Nnode)),]
  x@phylo$tip.label <- x@phylo$tip.label[e2[e2 <= nTip]]
  x@phylo$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# idem but with phylo
sort_tips_phylo <- function(x) {
  x <- ape::reorder.phylo(x)
  nTip <- as.integer(length(x$tip.label))
  e2 <- x$edge[, 2]
  x$tip.label <- x$tip.label[e2[e2 <= nTip]]
  x$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# get MRCA height from tree(s)
get_MRCA_heights <- function(x) {
  fun <- function(t) max(ape::node.depth.edgelength(t))
  height <- NULL
  if (inherits(x, "phylo")) height <- fun(x)
  if (inherits(x, "multiPhylo")) {
    if (!is.null(attr(x, "TipLabel"))) {
      x <- ape::.uncompressTipLabel(x)
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
    else {
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
  }
  else {
    height <- vapply(x, fun, 0)
  }
  height
}

# add tip labels to a tree plot - copied from phangorn
add_tiplabels <- function(xy, tip.label, direction, adj, font, srt = 0, cex = 1,
                          col = 1, label_offset = 0) {
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards",
                                      "downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")
  nTips <- length(tip.label)
  xx <- rep(1, nrow(xy))
  yy <- xy[, 2 ]
  if (direction == "leftwards" | direction == "downwards") xx <- xx * 0
  if (!horizontal) {
    #    tmp <- yy
    yy <- xx
    xx <- xy[, 1]
  }
  MAXSTRING <- max(strwidth(tip.label, cex = cex))
  loy <- 0
  if (direction == "rightwards") lox <- label_offset + MAXSTRING * 1.05 * adj
  if (direction == "leftwards")
    lox <- -label_offset - MAXSTRING * 1.05 * (1 - adj)
  if (!horizontal) {
    psr <- par("usr")
    MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
    loy <- label_offset + MAXSTRING * 1.05 * adj
    lox <- 0
    srt <- 90 + srt
    if (direction == "downwards") {
      loy <- -loy
      srt <- 180 + srt
    }
  }
  text(xx[1:nTips] + lox, yy[1:nTips] + loy, tip.label, adj = adj,
       font = font, srt = srt, cex = cex, col = col)
}

# adapted from treeplyr (package no longer available on CRAN)
reorder_treedata <- function(tdObject, order = "postorder") {
  dat.attr <- attributes(tdObject@data)
  phy <- tdObject@phylo
  ntips <- length(phy$tip.label)
  phy$node.label <- (ntips+1):(ntips+phy$Nnode)
  phy <- ape::reorder.phylo(phy, order)
  index <- match(tdObject@phylo$tip.label, phy$tip.label)
  index.node <- match((ntips+1):(ntips+phy$Nnode), phy$node.label)
  
  tdObject@data <- tdObject@data[c(index,index.node),]
  attributes(tdObject@data) <-dat.attr
  attributes(tdObject)$tip.label <- phy$tip.label
  tdObject@phylo <- phy
  
  tdObject
}

## End functions required by densiTreeWithBranchData

.removeNull <- function(x){
  res <- x[which(!sapply(x, is.null))]
}

# set prob factors
.set_pp_factor_range <- function(t, include_start_states, n_states=1) {

  # what is the ancestral state name tag?
  if (include_start_states) {
    state_pos_str_base = c("start_state_", "end_state_")
  } else {
    state_pos_str_base = c("anc_state_")
  }

  # create list of ancestral state name tags
  state_pos_str_to_update = c(sapply(1:n_states, function(x) { paste(state_pos_str_base,x,"_pp",sep="")}))

  # overwrite state labels
  for (m in state_pos_str_to_update)
  {
    x_state = attributes(t)$data[[m]]
    #levels(x_state) = c(levels(x_state))
    attributes(t)$data[[m]] = x_state
  }
  return(t)
}

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

