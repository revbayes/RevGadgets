.add_epoch_times <- function(p, max_age, dy_bars, dy_text) {
  max_x <- max(p$data$x)
  max_y <- max(p$data$y)

  epoch_names <- rev(
    c(
      "Holocene",
      "Pleistocene",
      "Pliocene",
      "Miocene",
      "Oligocene",
      "Eocene",
      "Paleocene",
      "Upper\nCretaceous",
      "Lower\nCretaceous"
    )
  )
  epoch_ages <- rev(c(0, 0.117, 2.588, 5.332, 23.03,
                      33.9, 55.8, 65.5, 99.6))
  period_names <- rev(
    c(
      "Quaternary",
      "Neogene",
      "Paleogene",
      "Cretaceous",
      "Jurassic",
      "Triassic",
      "Permian",
      "Carboniferous",
      "Devonian",
      "Silurian",
      "Ordovician",
      "Cambrian"
    )
  )

  period_ages <- rev(c(
    0,
    2.588,
    23.03,
    65.5,
    145.5,
    199.6,
    251,
    299,
    359.2,
    416,
    443.7,
    488.3
  ))
  if (max_age > 140) {
    x_pos <- max_x - c(max_age, period_ages)
  } else
    x_pos <- max_x - c(max_age, epoch_ages)

  y_pos <- rep(max_y, length(x_pos))
  x_pos_mid <-
    (x_pos[1:(length(x_pos) - 1)] + x_pos[2:length(x_pos)]) / 2

  for (k in 2:(length(x_pos))) {
    box_col <- "gray92"
    if (k %% 2 == 0)
      box_col <- "white"
    box <-
      ggplot2::geom_rect(
        xmin = x_pos[k - 1],
        xmax = x_pos[k],
        ymin = dy_bars,
        ymax = y_pos[k],
        fill = box_col
      )
    p <- gginnards::append_layers(p, box, position = "bottom")
  }
  if (max_age > 140) {
    for (k in seq_len(length(period_names))) {
      p <-
        p + ggplot2::annotate(
          geom = "text",
          label = period_names[k],
          angle = 90,
          x = x_pos_mid[k],
          y = dy_text,
          hjust = 0,
          size = 3.25
        )
    }
  } else {
    for (k in seq_len(length(epoch_names))) {
      p <-
        p + ggplot2::annotate(
          geom = "text",
          label = epoch_names[k],
          angle = 90,
          x = x_pos_mid[k],
          y = dy_text,
          hjust = 0,
          size = 3.25
        )
    }
  }
  return(p)
}

# stolen from treeio: https://github.com/YuLab-SMU/treeio
.validTblTree <-
  function(object, cols = c("parent", "node", "label")) {
    cc <- cols[!cols %in% colnames(object)]
    if (length(cc) > 0) {
      msg <-
        paste0("invalid tbl_tree object.\n  missing column:\n    ",
               paste(cc, collapse = ","),
               ".")
    }
  }

# unexported ggtree function getXcoord:
# https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getXcoord <- function(tr) {
  edge <- tr$edge
  parent <- edge[, 1]
  child <- edge[, 2]
  root <- .rootNode(tr)

  len <- tr$edge.length

  N <- ape::Nnode(tr, internal.only = FALSE)
  x <- numeric(N)
  x <- .getXcoord2(x, root, parent, child, len)
  return(x)
}

# unexported ggtree function getXcoord2:
# https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getXcoord2 <-
  function(x,
           root,
           parent,
           child,
           len,
           start = 0,
           rev = FALSE) {
    x[root] <- start
    x[-root] <- NA  ## only root is set to start, by default 0

    currentNode <- root
    direction <- 1
    if (rev == TRUE) {
      direction <- -1
    }
    while (anyNA(x)) {
      idx <- which(parent %in% currentNode)
      newNode <- child[idx]
      x[newNode] <- x[parent[idx]] + len[idx] * direction
      currentNode <- newNode
    }

    return(x)
  }

# unexported ggtree function getYcoord:
# https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
.getYcoord <- function(tr, step = 1) {
  Ntip <- length(tr[["tip.label"]])
  N <- ape::Nnode(tr, internal.only = FALSE)

  edge <- tr[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]

  cl <- split(child, parent)
  child_list <- list()
  child_list[as.numeric(names(cl))] <- cl

  y <- numeric(N)
  tip.idx <- child[child <= Ntip]
  y[tip.idx] <- 1:Ntip * step
  y[-tip.idx] <- NA

  currentNode <- 1:Ntip
  while (anyNA(y)) {
    pNode <- unique(parent[child %in% currentNode])
    ## piping of magrittr is slower than nested function call.
    ## pipeR is fastest, may consider to use pipeR
    ##
    ## child %in% currentNode %>% which %>% parent[.] %>% unique
    ## idx <- sapply(pNode,
    ##              function(i) all(child[parent == i] %in% currentNode))
    idx <-
      unlist(lapply(pNode, function(i)
        all(child_list[[i]] %in% currentNode)))
    newNode <- pNode[idx]

    y[newNode] <- unlist(lapply(newNode, function(i) {
      mean(y[child_list[[i]]], na.rm = TRUE)
      ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
    }))

    currentNode <-
      c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
    ## currentNode <-
    ##    c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
    ## parent %in% newNode %>% child[.] %>%
    ##     `%in%`(currentNode, .) %>% `!` %>%
    ##         currentNode[.] %>% c(., newNode)
  }

  return(y)
}
