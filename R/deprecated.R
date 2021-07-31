# .add_epoch_times <- function(p, max_age, dy_bars, dy_text) {
#   max_x <- max(p$data$x)
#   max_y <- max(p$data$y)
#
#   epoch_names <- rev(
#     c(
#       "Holocene",
#       "Pleistocene",
#       "Pliocene",
#       "Miocene",
#       "Oligocene",
#       "Eocene",
#       "Paleocene",
#       "Upper\nCretaceous",
#       "Lower\nCretaceous"
#     )
#   )
#   epoch_ages <- rev(c(0, 0.117, 2.588, 5.332, 23.03,
#                       33.9, 55.8, 65.5, 99.6))
#   period_names <- rev(
#     c(
#       "Quaternary",
#       "Neogene",
#       "Paleogene",
#       "Cretaceous",
#       "Jurassic",
#       "Triassic",
#       "Permian",
#       "Carboniferous",
#       "Devonian",
#       "Silurian",
#       "Ordovician",
#       "Cambrian"
#     )
#   )
#
#   period_ages <- rev(c(
#     0,
#     2.588,
#     23.03,
#     65.5,
#     145.5,
#     199.6,
#     251,
#     299,
#     359.2,
#     416,
#     443.7,
#     488.3
#   ))
#   if (max_age > 140) {
#     x_pos <- max_x - c(max_age, period_ages)
#   } else
#     x_pos <- max_x - c(max_age, epoch_ages)
#
#   y_pos <- rep(max_y, length(x_pos))
#   x_pos_mid <-
#     (x_pos[1:(length(x_pos) - 1)] + x_pos[2:length(x_pos)]) / 2
#
#   for (k in 2:(length(x_pos))) {
#     box_col <- "gray92"
#     if (k %% 2 == 0)
#       box_col <- "white"
#     box <-
#       ggplot2::geom_rect(
#         xmin = x_pos[k - 1],
#         xmax = x_pos[k],
#         ymin = dy_bars,
#         ymax = y_pos[k],
#         fill = box_col
#       )
#     p <- gginnards::append_layers(p, box, position = "bottom")
#   }
#   if (max_age > 140) {
#     for (k in seq_len(length(period_names))) {
#       p <-
#         p + ggplot2::annotate(
#           geom = "text",
#           label = period_names[k],
#           angle = 90,
#           x = x_pos_mid[k],
#           y = dy_text,
#           hjust = 0,
#           size = 3.25
#         )
#     }
#   } else {
#     for (k in seq_len(length(epoch_names))) {
#       p <-
#         p + ggplot2::annotate(
#           geom = "text",
#           label = epoch_names[k],
#           angle = 90,
#           x = x_pos_mid[k],
#           y = dy_text,
#           hjust = 0,
#           size = 3.25
#         )
#     }
#   }
#   return(p)
# }
#
# # stolen from treeio: https://github.com/YuLab-SMU/treeio
# .validTblTree <-
#   function(object, cols = c("parent", "node", "label")) {
#     cc <- cols[!cols %in% colnames(object)]
#     if (length(cc) > 0) {
#       msg <-
#         paste0("invalid tbl_tree object.\n  missing column:\n    ",
#                paste(cc, collapse = ","),
#                ".")
#     }
#   }
#
# # unexported ggtree function getXcoord:
# # https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
# .getXcoord <- function(tr) {
#   edge <- tr$edge
#   parent <- edge[, 1]
#   child <- edge[, 2]
#   root <- .rootNode(tr)
#
#   len <- tr$edge.length
#
#   N <- ape::Nnode(tr, internal.only = FALSE)
#   x <- numeric(N)
#   x <- .getXcoord2(x, root, parent, child, len)
#   return(x)
# }
#
# # unexported ggtree function getXcoord2:
# # https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
# .getXcoord2 <-
#   function(x,
#            root,
#            parent,
#            child,
#            len,
#            start = 0,
#            rev = FALSE) {
#     x[root] <- start
#     x[-root] <- NA  ## only root is set to start, by default 0
#
#     currentNode <- root
#     direction <- 1
#     if (rev == TRUE) {
#       direction <- -1
#     }
#     while (anyNA(x)) {
#       idx <- which(parent %in% currentNode)
#       newNode <- child[idx]
#       x[newNode] <- x[parent[idx]] + len[idx] * direction
#       currentNode <- newNode
#     }
#
#     return(x)
#   }
#
# # unexported ggtree function getYcoord:
# # https://github.com/GuangchuangYu/ggtree/blob/master/R/tree-utilities.R
# .getYcoord <- function(tr, step = 1) {
#   Ntip <- length(tr[["tip.label"]])
#   N <- ape::Nnode(tr, internal.only = FALSE)
#
#   edge <- tr[["edge"]]
#   parent <- edge[, 1]
#   child <- edge[, 2]
#
#   cl <- split(child, parent)
#   child_list <- list()
#   child_list[as.numeric(names(cl))] <- cl
#
#   y <- numeric(N)
#   tip.idx <- child[child <= Ntip]
#   y[tip.idx] <- 1:Ntip * step
#   y[-tip.idx] <- NA
#
#   currentNode <- 1:Ntip
#   while (anyNA(y)) {
#     pNode <- unique(parent[child %in% currentNode])
#     ## piping of magrittr is slower than nested function call.
#     ## pipeR is fastest, may consider to use pipeR
#     ##
#     ## child %in% currentNode %>% which %>% parent[.] %>% unique
#     ## idx <- sapply(pNode,
#     ##              function(i) all(child[parent == i] %in% currentNode))
#     idx <-
#       unlist(lapply(pNode, function(i)
#         all(child_list[[i]] %in% currentNode)))
#     newNode <- pNode[idx]
#
#     y[newNode] <- unlist(lapply(newNode, function(i) {
#       mean(y[child_list[[i]]], na.rm = TRUE)
#       ##child[parent == i] %>% y[.] %>% mean(na.rm=TRUE)
#     }))
#
#     currentNode <-
#       c(currentNode[!currentNode %in% unlist(child_list[newNode])], newNode)
#     ## currentNode <-
#     ##    c(currentNode[!currentNode %in% child[parent %in% newNode]], newNode)
#     ## parent %in% newNode %>% child[.] %>%
#     ##     `%in%`(currentNode, .) %>% `!` %>%
#     ##         currentNode[.] %>% c(., newNode)
#   }
#
#   return(y)
# }
#
# # modified for revbayes
# # from unexported treeio function read.stats_beast_internal
# .read.stats_revbayes_internal <- function(beast, tree) {
#   phylo <- ape::read.tree(text = tree) # read the tree
#   tree2 <-
#     .add_pseudo_nodelabel(phylo) # add nodelabels (if there aren't already any)
#
#   ## node name corresponding to stats
#   nn <- unlist(strsplit(unlist(strsplit(tree2, split = ",")), "\\)"))
#   nn <- gsub("\\(*", "", nn)
#   nn <- gsub("[:;].*", "", nn)
#   nn <- gsub(" ", "", nn)
#   nn <- gsub("'", "", nn)
#   nn <- gsub('"', "", nn)
#
#   # nn <- strsplit(tree2, split=",") %>% unlist %>%
#   #   strsplit(., split="\\)") %>% unlist %>%
#   #   gsub("\\(*", "", .) %>%
#   #   gsub("[:;].*", "", .) %>%
#   #   gsub(" ", "", .) %>%
#   #   gsub("'", "", .) %>%
#   #   gsub('"', "", .)
#   # get the label of each node (internal and external) in the order
#   # they appear in the newick string
#
#   phylo <- ape::read.tree(text = tree2)
#   root <- .rootNode(phylo)
#   nnode <- phylo$Nnode
#
#   tree_label <- c(phylo$tip.label, phylo$node.label)
#   ii <- match(nn, tree_label)
#
#   if (any(grepl("TRANSLATE", beast, ignore.case = TRUE))) {
#     label2 <- c(phylo$tip.label, root:treeio::getNodeNum(phylo))
#   } else {
#     label2 <- as.character(1:treeio::getNodeNum(phylo))
#   }
#   node <- label2[match(nn, tree_label)]
#
#   ## stats <- unlist(strsplit(tree, "\\["))[-1]
#   ## stats <- sub(":.+$", "", stats
#
#   ## BEAST1 edge stat fix
#   tree <- gsub("\\]:\\[&(.+?\\])", ",\\1:", tree)
#   tree <- gsub(":(\\[.+?\\])", "\\1:", tree)
#
#   if (grepl("\\]:[0-9\\.eE+\\-]*\\[", tree) ||
#       grepl("\\]\\[", tree)) {
#     ## MrBayes output
#     # stats <- strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\[") %>% unlist
#     stats <- unlist(strsplit(tree, "\\]:[0-9\\.eE+\\-]*\\["))
#     lstats <- lapply(stats, function(x) {
#       unlist(strsplit(x, split = "\\][,\\)]"))
#     })
#
#     for (i in seq_along(stats)) {
#       n <- length(lstats[[i]])
#       if (i == length(stats)) {
#         stats[i] <- lstats[[i]][n]
#       } else {
#         stats[i] <- paste0(lstats[[i]][n],
#                            sub("&", ",", lstats[[i + 1]][1]))
#       }
#     }
#     stats <- gsub("\\]\\[&", ",", stats)
#   } else {
#     ## BEAST output
#     stats <- unlist(strsplit(tree, ":"))
#   }
#
#   names(stats) <- node
#
#   stats <- stats[grep("\\[", stats)]
#   stats <- sub("[^\\[]*\\[", "", stats)
#
#   stats <- sub("^&", "", stats)
#   stats <- sub("];*$", "", stats)
#   stats <- gsub("\"", "", stats)
#
#   #this is what is breaking readTrees for the OU output
#   stats2 <- lapply(seq_along(stats), function(i) {
#     x <- stats[[i]]
#     y <- unlist(strsplit(x, ","))
#     sidx <- grep("=\\{", y)
#     eidx <- grep("\\}$", y)
#
#     flag <- FALSE
#     if (length(sidx) > 0) {
#       flag <- TRUE
#       # SETS <- lapply(seq_along(sidx), function(k) {
#       #   p <- y[sidx[k]:eidx[k]]
#       #   gsub(".*=\\{", "", p) %>% gsub("\\}$", "", .)
#       # })
#       SETS <- lapply(seq_along(sidx), function(k) {
#         p <- y[sidx[k]:eidx[k]]
#         gsub("\\}$", "", gsub(".*=\\{", "", p))
#       })
#       names(SETS) <- gsub("=.*", "", y[sidx])
#
#       kk <- unlist(lapply(seq_along(sidx), function(k) {
#         sidx[k]:eidx[k]
#       }))
#       y <- y[-kk]
#     }
#
#     if (length(y) == 0)
#       return(SETS)
#
#     name <- gsub("=.*", "", y)
#     val <- gsub(".*=", "", y)
#     val <- gsub("^\\{", "", val)
#     val <- gsub("\\}$", "", val)
#
#     # %>%
#     #   gsub("^\\{", "", .) %>%
#     #   gsub("\\}$", "", .)
#
#     if (flag) {
#       nn <- c(name, names(SETS))
#     } else {
#       nn <- name
#     }
#
#     res <- rep(NA, length(nn))
#     names(res) <- nn
#
#     for (i in seq_along(name)) {
#       # res[i] <- if(treeio:::is.numeric(val[i])) as.numeric(val[i]) else val[i]
#       res[i] <-
#         if (is.numeric(val[i]))
#           as.numeric(val[i])
#       else
#         val[i]
#     }
#     if (flag) {
#       j <- i
#       for (i in seq_along(SETS)) {
#         if (is.numeric(SETS[[i]])) {
#           res[i + j] <- list(as.numeric(SETS[[i]]))
#         } else {
#           res[i + j] <- SETS[i]
#         }
#       }
#     }
#
#     return(res)
#   })
#
#   nn <- sort(unique(unlist(lapply(stats2, names))))
#
#   # nn <- lapply(stats2, names) %>% unlist %>%
#   #   unique %>% sort
#
#   stats2 <- lapply(stats2, function(x) {
#     y <- x[nn]
#     names(y) <- nn
#     y[vapply(y, is.null, logical(1))] <- NA
#     y
#   })
#
#   stats3 <- do.call(rbind, stats2)
#   stats3 <- tibble::as_tibble(stats3)
#
#   ## no need to extract sd from prob+-sd
#   ## as the sd is stored in prob_stddev
#   ##
#   ## "prob_stddev"   "prob(percent)" "prob+-sd"
#   ##
#   ##
#   ##
#   ## idx <- grep("\\+-", colnames(stats3))
#   ## if (length(idx)) {
#   ##     for (i in idx) {
#   ##         stats3[,i] <- as.numeric(gsub("\\d+\\+-", "", stats3[,i]))
#   ##     }
#   ## }
#
#   cn <- gsub("(\\d+)%", "0.\\1", colnames(stats3))
#   cn <- gsub("\\(([^\\)]+)\\)", "_\\1", cn)
#   ## cn <- gsub("\\+-", "_", cn)
#
#   colnames(stats3) <- cn
#   stats3$node <- names(stats)
#
#   i <- vapply(stats3,
#               function(x)
#                 max(vapply(x, length, numeric(1))),
#               numeric(1))
#
#   for (j in which(i == 1)) {
#     stats3[, j] <- unlist(stats3[, j])
#   }
#   stats3$node <- as.integer(stats3$node)
#   return(stats3)
# }
#
# # unexported treeio function BEAST:
# # https://github.com/YuLab-SMU/treeio
# .beast <- function(file, treetext, stats, phylo) {
#   stats$node <- gsub("\"*'*", "", stats$node)
#
#   phylo <- .remove_quote_in_tree_label(phylo)
#
#   obj <- methods::new(
#     "treedata",
#     ## fields      = fields,
#     treetext    = treetext,
#     phylo       = phylo,
#     data        = stats,
#     file        = .filename(file)
#   )
#
#   return(obj)
# }

# # unexported treeio function filename:
# # https://github.com/YuLab-SMU/treeio
# .filename <- function(file) {
#   ## textConnection(text_string) will work just like a file
#   ## in this case, just set the filename as ""
#   file_name <- ""
#   if (is.character(file)) {
#     file_name <- file
#   }
#   return(file_name)
# }

# # unexported treeio function remove_quote_in_tree_label:
# # https://github.com/YuLab-SMU/treeio
# .remove_quote_in_tree_label <- function(phylo) {
#   if (!is.null(phylo$node.label)) {
#     phylo$node.label <- gsub("\"*'*", "", phylo$node.label)
#   }
#   if (!is.null(phylo$tip.label)) {
#     phylo$tip.label <- gsub("\"*'*", "", phylo$tip.label)
#   }
#   return(phylo)
# }

# # Unexported treeio function rootnode.phylo:
# # https://github.com/YuLab-SMU/treeio
# # works for phylo objects, not tree data
# .rootNode <- function(.data, ...) {
#   edge <- .data[["edge"]]
#   ## 1st col is parent,
#   ## 2nd col is child,
#   if (!is.null(attr(.data, "order")) &&
#       attr(.data, "order") == "postorder")
#     return(edge[nrow(edge), 1])
#
#   parent <- unique(edge[, 1])
#   child <- unique(edge[, 2])
#   ## the node that has no parent should be the root
#   root <- parent[!parent %in% child]
#   if (length(root) > 1) {
#     stop("multiple roots found...")
#   }
#   return(root)
# }

# # unexported treeio function add_pseudo_nodelabel: h
# # ttps://github.com/YuLab-SMU/treeio
# .add_pseudo_nodelabel <- function(phylo) {
#   if (is.null(phylo$node.label)) {
#     nnode <- phylo$Nnode
#     phylo$node.label <- paste("X", 1:nnode, sep = "")
#   }
#   ## if tip.label contains () which will broken node name extraction
#   phylo$tip.label <- gsub("[\\(\\)]", "_", phylo$tip.label)
#
#   treetext <- ape::write.tree(phylo)
#   return(treetext)
# }

# # Fast data.frame constructor and indexing
# # No checking, recycling etc. unless asked for
# # unexported ggplot2 function new_data_frame
# .new_data_frame <- function(x = list(), n = NULL) {
#   if (length(x) != 0 && is.null(names(x))) {
#     stop("Elements must be named")
#   }
#   lengths <- vapply(x, length, integer(1))
#   if (is.null(n)) {
#     n <- if (length(x) == 0 || min(lengths) == 0)
#       0
#     else
#       max(lengths)
#   }
#   for (i in seq_along(x)) {
#     if (lengths[i] == n)
#       next
#     if (lengths[i] != 1) {
#       stop("Elements must equal the number of rows or 1")
#     }
#     x[[i]] <- rep(x[[i]], n)
#   }
#
#   class(x) <- "data.frame"
#
#   attr(x, "row.names") <- .set_row_names(n)
#   x
# }

# # modified from
# # https://github.com/GuangchuangYu/ggimage/blob/master/R/geom_subview.R
# .geom_subview_revgadgets <-
#   function (mapping = NULL,
#             data = NULL,
#             width = 0.1,
#             height = 0.1,
#             x = NULL,
#             y = NULL,
#             subview = NULL) {
#     if (is.null(data)) {
#       data <- dplyr::tibble(x = x, y = y)
#     } else if (!inherits(data, "tbl")) {
#       data <- dplyr::as_tibble(data)
#     }
#     if (is.null(mapping)) {
#       mapping <- ggplot2::aes_(x = ~ x, y = ~ y)
#     }
#
#     mapping <- as.list(mapping)
#
#     if (is.null(mapping$x)) {
#       stop("x aesthetic mapping should be provided")
#     }
#
#     if (is.null(mapping$y)) {
#       stop("y aesthetic mapping should be provided")
#     }
#
#     if (is.null(mapping$subview) && is.null(subview)) {
#       stop("subview must be provided")
#     }
#
#     if (is.null(mapping$subview)) {
#       if (!inherits(subview, "list")) {
#         subview <- list(subview)
#       }
#       data$subview <- subview
#     } else {
#       sv_var <- rvcheck::get_aes_var(mapping, "subview")
#       data$subview <- data[[sv_var]]
#     }
#
#     xvar <- rvcheck::get_aes_var(mapping, "x")
#     yvar <- rvcheck::get_aes_var(mapping, "y")
#
#     if (is.null(mapping$width)) {
#       data$width <- width
#     } else {
#       width_var <- rvcheck::get_aes_var(mapping, "width")
#       data$width <- data[[width_var]]
#     }
#
#     if (is.null(mapping$height)) {
#       data$height <- height
#     } else {
#       height_var <- rvcheck::get_aes_var(mapping, "height")
#       data$height <- data[[height_var]]
#     }
#
#     data$xmin <- data[[xvar]] - data$width / (2 * max(data[[yvar]]))
#     data$xmax <- data[[xvar]] + data$width / (2 * max(data[[yvar]]))
#     data$ymin <- data[[yvar]] - data$width / (2 * max(data[[xvar]]))
#     data$ymax <- data[[yvar]] + data$width / (2 * max(data[[xvar]]))
#
#     # save pies as images and plot as raster grobs
#     results <- list()
#     for (i in seq_len(nrow(data))) {
#       ggplot2::ggsave(
#         ".temp.png",
#         plot = data$subview[[i]],
#         bg = "transparent",
#         width = 3,
#         height = 3,
#         units = "cm",
#         dpi = 200
#       )
#       pie <- png::readPNG(".temp.png")
#       g <- grid::rasterGrob(pie, interpolate = TRUE)
#       results[[i]] <-
#         ggplot2::annotation_custom(
#           ggplotify::as.grob(g),
#           xmin = data$xmin[i],
#           xmax = data$xmax[i],
#           ymin = data$ymin[i],
#           ymax = data$ymax[i]
#         )
#     }
#     file.remove(".temp.png")
#     return(results)
#
#     #old way of plotting pies - won't plot centered on nodes
#     #lapply(1:nrow(data), function(i) {
#     #  ggplot2::annotation_custom(
#     #    ggplotify::as.grob(data$subview[[i]]),
#     #    xmin = data$xmin[i],
#     #    xmax = data$xmax[i],
#     #    ymin = data$ymin[i],
#     #    ymax = data$ymax[i]
#     #  )
#     #})
#   }

# # modified from
# #https://github.com/YuLab-SMU/ggtree/blob/0681b23fe6afb510e5a6041a7cbf50c3b18473e8/R/inset.R
# .inset.revgadgets <-
#   function (tree_view,
#             insets,
#             width = 0.1,
#             height = 0.1,
#             hjust = 0,
#             vjust = 0,
#             x = "node",
#             pos = 0.5) {
#     df <- tree_view$data[as.numeric(names(insets)),]
#
#     # position subviews based on tree part
#     x <- match.arg(x, c("node", "branch", "edge", "parent_shoulder"))
#     if (x == "node") {
#       xx <- df$x
#     } else if (x == "parent_shoulder") {
#       xx <- df$x[match(df$parent, df$node)]
#     } else {
#       xx <- df$branch
#     }
#     yy <- df$y
#     xx <- xx - hjust # x-coordinates for nodes
#     yy <- yy - vjust # y-coordinates for nodes
#
#     if (length(width) == 1) {
#       width <- rep(width, length(insets))
#     }
#     if (length(height) == 1) {
#       height <- rep(height, length(insets))
#     }
#
#     tree_view <- tree_view +
#       .geom_subview_revgadgets(
#         subview = insets,
#         width = width,
#         height = height,
#         x = xx,
#         y = yy
#       )
#     # return treeview with subviews
#     return(tree_view)
#   }
