#' Plot tree
#'
#' Plots a single tree, such as an MCC or MAP tree.
#'
#' Plots a single tree, such as an MCC or MAP tree, with
#' optionally labeled posterior probabilities at nodes, a
#' timescale plotted on the x - axis, and 95\% CI for node ages.
#'
#'
#' @param tree (list of lists of treedata objects; no default) Name of a list of lists of
#' treedata objects, such as produced by readTrees(). This object should only contain
#' only one summary tree from one trace file. If it contains multiple trees or multiple
#' traces, only the first will be used.
#'
#' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with timescale in MYA.
#'
#' @param node_age_bars (logical; TRUE) Plot time tree with node age bars?
#'
#' @param node_age_bars_colored_by (character; NULL) Specify column to color node age bars by,
#' such as "posterior". If null, all node age bars plotted the same color, specified by
#' node_age_bars_color
#'
#' @param node_age_bars_color (character; "blue") Color for node age bars. If node_age_bars_colored_by
#' specifies a variable (not NULL), you must provide two colors, low and high values for a gradient. Colors must be either
#' R valid color names or valid hex codes.
#'
#' @param node_labels (character; NULL) Plot text labels at nodes, specified by the name of the
#' corresponding column in the tidytree object. If NULL, no text is plotted.
#'
#' @param node_labels_color (character; "black") Color to plot node_labels, either as a valid
#' R color name or a valid hex code.
#'
#' @param node_labels_size (numeric; 3) Size of node labels
#'
#' @param tip_labels (logical; TRUE) Plot tip labels?
#'
#' @param tip_labels_italics (logical; FALSE) Plot tip labels in italics?
#'
#' @param tip_labels_remove_underscore (logical; TRUE) Remove underscores in tip labels?
#'
#' @param tip_labels_color (character; "black") Color to plot tip labels, either as a valid
#' R color name or a valid hex code.
#'
#' @param tip_labels_size (numeric; 3) Size of tip labels
#'
#' @param scale_bar (logical; FALSE) Plot a scale bar in branch length units?
#'
#' @param node_pp (logical; FALSE) Plot posterior probabilities as symbols at nodes? Specify
#' symbol aesthetics with node_pp_shape, node_pp_color, and node_pp_size.
#'
#' @param node_pp_shape (integer; 1) Integer corresponding to point shape (value between 0-25).
#'
#' @param node_pp_color (character; "black") Color for node_pp symbols, either as valid R color name(s)
#' or hex code(s). Can be a single character string specifying a single color, or a vector of
#' length two specifying two colors to form a gradient. In this case, posterior probabilities
#' will be indicated by color along the specified gradient.
#'
#' @param node_pp_size (numeric or character; 1) Size for node_pp symbols. If numeric, the size will
#' be fixed at the specified value. If a character, it should specify "variable", indicating that
#' size should be scaled by the posterior value.
#'
#' @param branch_color (character; "black") A single character string specifying the color (R color
#' name or hex code) for all branches OR a vector of length 2 specifying two colors for a gradient,
#' used to color the branches according to the variable specified in color_branch_by.
#' If only 1 color is provided and you specify color_branch_by, default colorswill be chosen
#' (low = "#005ac8", high = "#fa7850").
#'
#' @param color_branch_by (character; NULL ) Optional name of one quantitative variable in the
#' treedata object to color branches, such as a rate.
#'
#' @param line_width (numeric; 1) Change line width for branches
#'
#' @return returns a single plot object.
#'
#' @examples
#' \dontrun{
#' # Example of standard tree plot
#' # Add on a scale bar using ggtree, but note that the x axis position must
#' # be negative and scales with tree height
#'
#' file <- system.file("extdata", "sub_models/primates_cytb_GTR_MAP.tre", package="RevGadgets")
#' tree <- readTrees(paths = file)
#' # Reroot tree before plotting
#' tree_rooted <- rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")
#' # Plot
#' plot <- plotTree(tree = tree_rooted, node_age_bars = FALSE, node_pp = F, node_labels = "posterior",
#'                  tip_labels_remove_underscore = T, node_labels_size = 3,
#'                  tip_labels_italics = F) + ggtree::geom_treescale(x = -0.1, y = 0)
#'
#' # Example of coloring branches by rate
#' file <- system.file("extdata", "relaxed_ou/relaxed_OU_MAP.tre", package="RevGadgets")
#' tree <- readTrees(paths = file)
#' plot <- plotTree(tree = tree, node_age_bars = FALSE, node_pp = F,
#'                  tip_labels_remove_underscore = T, tip_labels_italics = F,
#'                  color_branch_by = "branch_thetas", line_width = 1.7) +
#'        ggplot2::theme(legend.position=c(.1, .9))
#' }
#' @export


plotTree <- function(tree, timeline = FALSE, node_age_bars = TRUE, node_age_bars_color = "blue", node_age_bars_colored_by = NULL,
                     node_labels = NULL, node_labels_color = "black", node_labels_size = 3, tip_labels = TRUE,
                     tip_labels_italics = TRUE, tip_labels_remove_underscore = FALSE, tip_labels_color = "black",
                     tip_labels_size = 3, scale_bar = FALSE, node_pp = FALSE, node_pp_shape = 16, node_pp_color = "black",
                     node_pp_size = "variable", branch_color = "black", color_branch_by = NULL, line_width = 1) {
  # enforce argument matching
  if (!is.list(tree)) stop("tree should be a list of lists of treedata objects")
  if (class(tree[[1]][[1]]) != "treedata") stop("tree should be a list of lists of treedata objects")
  vars <- colnames(tree[[1]][[1]]@data)
  if (is.logical(timeline) == FALSE) stop("timeline should be TRUE or FALSE")
  if (is.logical(node_age_bars) == FALSE) stop("node_age_bars should be TRUE or FALSE")
  if (.isColor(node_age_bars_color) == FALSE) stop("node_age_bars_color should be valid color(s)")
  if (is.null(node_age_bars_colored_by) == FALSE &
      any(vars %in% node_age_bars_colored_by) == FALSE) stop("node_age_bars_colored_by should be a column in your tidytree object")
  if (is.null(node_labels) == FALSE &
      any(vars %in% node_labels) == FALSE) stop("node_labels should be NULL or a column in your tidytree object")
  if (is.null(node_labels_color) == FALSE & .isColor(node_labels_color) == FALSE) stop("node_labels_color should be NULL or a recognized color")
  if (is.logical(tip_labels) == FALSE) stop("tip_labels should be TRUE or FALSE")
  if (is.logical(tip_labels_italics) == FALSE) stop("tip_labels_italics should be TRUE or FALSE")
  if (.isColor(tip_labels_color) == FALSE) stop("tip_labels_color should be a recognized color")
  if (class(node_pp) != "logical") stop("node_pp should be TRUE or FALSE")
  if (node_pp) {
    if (length(node_pp_color) > 2) stop("node_pp_color should be of length 1 or 2")
    if (.isColor(node_pp_color) == FALSE) stop("node_pp_color should be a recognized color")
    if (node_pp_shape %in% 0:25 == FALSE) stop("node_pp_shape should be a recognized shape (value between 0 and 25)")
    if (is.numeric(node_pp_size) == FALSE & node_pp_size != "variable") stop("node_pp_size should be numeric or 'variable'")
  }
  if (length(branch_color) == 1 & !.isColor(branch_color)) stop("branch_color should be a recognized color")
  if (length(branch_color) == 2) {
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2]) == FALSE) stop("Neither values of branch_color are a recognized color")
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2])) stop("branch_color[1] is not a recognized color")
    if (.isColor(branch_color[1]) &
        .isColor(branch_color[2]) == FALSE) stop("branch_color[2] is not a recognized color")
  } else if (length(branch_color) > 2 ) stop("only 2 colors may be specified in branch_color")
  if (is.null(color_branch_by) == FALSE &
      any(vars %in% color_branch_by) == FALSE) stop("color_branch_by should be NULL or a column in your tidytree object")
  if (is.numeric(line_width) == FALSE) stop ("line_width should be numeric")

  # grab single tree from input
  phy <- tree[[1]][[1]]

    # initiate plot
  if (is.null(color_branch_by)) {
    pp <- ggtree::ggtree(phy, right = F, size = line_width, color = branch_color)
  } else if (!is.null(color_branch_by)) {
    pp <- ggtree::ggtree(phy, right = F, size = line_width)
  }


  #### paramter compatibility checks ###
  if (length(node_pp_color) == 2 & length(branch_color) == 2) stop("You may only include variable colors for either
                                                                   node_pp_label or branch_color, not for both")
  #check that if user wants node_age_bars tree, there are dated intervals in the file
  if (node_age_bars == TRUE) {
    if(!"age_0.95_HPD" %in% colnames(phy@data)) stop("You specified node_age_bars, but there is no age_0.95_HPD column in the treedata object.")
  }

  # get dimensions
  n_nodes <- treeio::Nnode(phy)
  tree_height <- max(phytools::nodeHeights(phy@phylo))
  ntips <- sum(pp$data$isTip)

  # reformat labels if necessary
  if (tip_labels_remove_underscore) { pp$data$label <- gsub("_", " ", pp$data$label)}

    # add timeline
  if (timeline == TRUE) {


    pp$data$age_0.95_HPD <- lapply(pp$data$age_0.95_HPD, function(z) {
      if (is.null(z) || is.na(z)) { return(c(NA,NA)) } else { return(as.numeric(z)) }
    })
    minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
    if (node_age_bars == FALSE) {
      max_age <- max(ape::branching.times(phy@phylo))
    } else {
      max_age <- max(minmax, na.rm =TRUE)
    }

    if (max_age > 100){
      interval <- 50
    } else {interval <- 10}

    dx <- max_age %% interval

        # set coordinates
    ### fix the xlim and ylims - if no error bars, should be a function of max age and n nodes, respectively
    ### if error bars, -x lim should be as old as the max of the error bar
    #pp <- pp + ggplot2::coord_cartesian(xlim = c(-max_age,30), ylim=c(-7, n_nodes+1.5), expand=F)
    pp <- pp + ggplot2::coord_cartesian()
    pp <- pp + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                           limits = c(-max(minmax, na.rm = T)*1.05, tree_height/2),
                                           breaks = -rev(seq(0,max_age+dx,interval)),
                                           labels = rev(seq(0,max_age+dx,interval)),
                                           )
    pp <- pp + ggtree::theme_tree2()
    #pp <- pp + ggplot2::theme(legend.position=c(.05, .955), axis.line = ggplot2::element_line(colour = "black"))
    pp <- ggtree::revts(pp)
    pp <- .add_epoch_times(pp, max_age, dy_bars=-7, dy_text=-3)
  }

    # add scale bar
  if (scale_bar == TRUE) {pp <- pp + ggtree::geom_treescale()}

  # processing for node_age_bars and tip_age_bars
  if (node_age_bars == TRUE) {
    # Encountered problems with using geom_range to plot age HPDs in ggtree. It
    # appears that geom_range incorrectly rotates the HPD relative to the height
    # of the node unnecessarily. My guess for this would be because older version
    # of ggtree primarily supported length measurements, and not height measurements
    # so the new capability to handle height might contain a "reflection" bug.
    #
    # For example, suppose a node has height 3 with HPD [2, 7]. You can think of
    # this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
    # causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
    #
    # See this excellent trick by Tauana:
    # https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
    # Adapted this code to also plot fossil tip uncertainty in red
    pp$data$age_0.95_HPD <- lapply(pp$data$age_0.95_HPD, function(z) {
      if (is.null(z) || is.na(z)) { return(c(NA,NA)) } else { return(as.numeric(z)) }
    })
    minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
    bar_df <- data.frame(node_id = as.integer(pp$data$node), isTip = pp$data$isTip, as.data.frame(minmax))
    names(bar_df) <- c("node_id", "isTip", "min", "max")
    node_df <- dplyr::filter(bar_df, isTip == FALSE)
    if (is.null(node_age_bars_colored_by) == TRUE) {
      # plot age densities
      node_df <- dplyr::left_join(node_df, pp$data, by=c("node_id"="node"))
      node_df <- dplyr::select(node_df,  node_id, min, max, y)
      pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y),
                                       data=node_df, color=node_age_bars_color, size=1.5, alpha=0.8)
     } else if (is.null(node_age_bars_colored_by) == FALSE) {
      pp$data$olena <- c(rep(NA, times = ntips),
                         as.numeric(.convertAndRound(L = unlist(pp$data[pp$data$isTip == FALSE,
                                                             node_age_bars_colored_by]))))
      node_df <- dplyr::left_join(node_df, pp$data, by=c("node_id"="node"))
      node_df <- dplyr::select(node_df,  node_id, min, max, y, olena)
      pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y, color = olena),
                                       data=node_df, size=1.5, alpha=0.8) +
                 ggplot2::scale_color_gradient(low = node_age_bars_color[1], high = node_age_bars_color[2],
                                               name = paste(.simpleCap(node_age_bars_colored_by)))
      }
  }

  # add node labels (text)
  if (is.null(node_labels) == FALSE) {

    #catch some funkiness from importing an unrooted tree
    if (node_labels == "posterior") {
      pp$data[grep("]:", unlist(pp$data[,node_labels])), node_labels] <- NA
    }
    pp$data$kula <- c(rep(NA, times = ntips),
                      .convertAndRound(L = unlist(pp$data[pp$data$isTip == FALSE,
                                                          node_labels])))
    #change any NAs that got converted to characters back to NA
    pp$data$kula[pp$data$kula == "NA"] <- NA
    pp <- pp + ggtree::geom_nodelab(ggplot2::aes(label = kula),
                                    geom = "text", color = node_labels_color,
                                    hjust = 0, size = node_labels_size)
  }

  # add tip labels (text)
  if (tip_labels == TRUE) {
    if (tip_labels_italics) {
      pp <- pp + ggtree::geom_tiplab(ggplot2::aes(label = paste0('italic(`', label, '`)')),
                                     size = tip_labels_size, offset = tree_height * 0.01,
                                     color = tip_labels_color, parse=TRUE)

    } else {
      pp <- pp + ggtree::geom_tiplab(size = tip_labels_size, offset = tree_height * 0.01,
                                     color = tip_labels_color)
      }
  }

  # add node PP (symbols)
  if (node_pp == TRUE) {
    # reformat posterior
    pp$data$posterior <- as.numeric(pp$data$posterior)

    if (length(node_pp_color) == 1 & node_pp_size == "variable") {
      pp <- pp + ggtree::geom_nodepoint(color = node_pp_color,
                                        ggplot2::aes(size = posterior),
                                        shape = node_pp_shape) +
        ggplot2::scale_size_continuous(name="Posterior")
    } else if (length(node_pp_color) == 2 & node_pp_size != "variable") {

      pp <- pp + ggtree::geom_nodepoint(size = node_pp_size,
                                        ggplot2::aes(color = posterior),
                                        shape = node_pp_shape) +
        ggplot2::scale_color_gradient(name="Posterior",
                                      low = node_pp_color[1], high = node_pp_color[2])
    }


  }

  # add branch coloration by variable
  if (is.null(color_branch_by) == FALSE) {

    #set default colors if none provided
    if (length(branch_color) != 2) {
      branch_color <- c("#005ac8", "#fa7850")
    }
     #col_num <- which(colnames(phy@data) == color_branch_by)
     #phy@data[,col_num] <- as.numeric(as.data.frame(phy@data)[,col_num]) #convert data to numeric
     #name <- .simpleCap(sub(pattern = "_", replacement = " ", color_branch_by))
     #pp <- pp +
     #  ggplot2::aes(color=I(as.data.frame(phy@data)[,col_num])) +
     #  ggplot2::scale_color_gradient(low = branch_color[1], high = branch_color[2],
     #                                name = name)

    col_num <- which(colnames(pp$data) == color_branch_by)
    pp$data[,col_num] <- as.numeric(as.data.frame(pp$data)[,col_num]) #convert data to numeric
    name <- .simpleCap(sub(pattern = "_", replacement = " ", color_branch_by))
    pp <- pp +
      ggplot2::aes(color=I(as.data.frame(pp$data)[,col_num])) +
      ggplot2::scale_color_gradient(low = branch_color[1], high = branch_color[2],
                                    name = name)

  }

  # readjust axis
  if (node_age_bars == FALSE & timeline == FALSE) {
    # add extra space on plot for tip labels
    tree_height <- max(phytools::nodeHeights(phy@phylo))
    pp <- pp + ggtree::xlim(-tree_height, tree_height/2)
    pp <- ggtree::revts(pp)
  } else if (node_age_bars == TRUE & timeline == FALSE & tip_labels == TRUE) {
    tree_height <- max(phytools::nodeHeights(phy@phylo))
    pp <- pp + ggtree::xlim(0, tree_height + tree_height/2)
  }

  # adjust legend(s)

  pp <- pp + ggplot2::theme(legend.position=c(.9, .8))
  return(pp)
}










#################### OLD VERSION ####################
# #' Plot tree
# #'
# #' Plots a single tree, such as an MCC or MAP tree.
# #'
# #' Plots a single tree, such as an MCC or MAP tree, with
# #' optionally labele posterior probabilities at nodes, a
# #' timescale plotted on the x - axis, and 95\% CI for node ages.
# #'
# #'
# #' @param tree (list of lists of treedata objects; no default) Name of a list of lists of
# #' treedata objects, such as produced by readTrees(). This object should only contain
# #' only one summary tree from one trace file. If it contains multiple trees or multiple
# #' traces, only the first will be used.
# #'
# #' @param node_age_bars (logical; FALSE) Plot time tree with node age bars.
# #'
# #' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with timescale in MYA.
# #'
# #' @param node_labels (logical; TRUE) Plot posterior probabilities as colored cirlces at
# #' nodes?
# #'
# #' @param tip_labels (logical; TRUE) Plot tip labels?
# #'
# #' @param color_branch (character; NULL ) Optional name of one quantitative variable in the
# #' treedata object to color branches, such as a rate.
# #'
# #' @param fossils (character vector; NULL) Optional vector of the tips that are fossils,
# #' to highlight age bars of fossils in red. Only works if node_age_bars is also TRUE.
# #'
# #' @param line_width (numeric; 1) Change line width for branches
# #'
# #' @return returns a single plot object.
# #'
# #' @examples
# #'
# #' @export
#
# plotTree <- function(tree, node_age_bars = FALSE, timeline = FALSE, node_labels = TRUE,
#                      tip_labels = TRUE, color_branch = NULL, fossils = NULL, line_width = 1) {
#   # enforce argument matching
#   if (!is.list(tree)) stop("tree should be a list of lists of treedata objects")
#   if (class(tree[[1]][[1]]) != "treedata") stop("tree should be a list of lists of treedata objects")
#   if (!is.logical(node_age_bars)) stop("node_age_bars should be TRUE or FALSE")
#   if (!is.logical(timeline)) stop("timeline should be TRUE or FALSE")
#   if (!is.logical(node_labels)) stop("node_labels should be TRUE or FALSE")
#   if (!is.logical(tip_labels)) stop("tip_labels should be TRUE or FALSE")
#   if (is.character(fossils) | is.null(fossils) == FALSE) stop("fossils should be NULL or character string or vector")
#
#   # grab single tree from input
#   phy <- tree[[1]][[1]]
#
#   #check that if user wants node_age_bars tree, there are dated intervals in the file
#   if (node_age_bars) {
#     if(!"age_0.95_HPD" %in% colnames(phy@data)) stop("You specified node_age_bars, but there is no age_0.95_HPD column in the treedata object.")
#   }
#
#    # initiate plot
#   pp <- ggtree::ggtree(phy, right = F, size = line_width)
#
#   # processing for node_age_bars
#   if (node_age_bars) {
#     # Encountered problems with using geom_range to plot age HPDs in ggtree. It
#     # appears that geom_range incorrectly rotates the HPD relative to the height
#     # of the node unnecessarily. My guess for this would be because older version
#     # of ggtree primarily supported length measurements, and not height measurements
#     # so the new capability to handle height might contain a "reflection" bug.
#     #
#     # For example, suppose a node has height 3 with HPD [2, 7]. You can think of
#     # this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
#     # causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
#     #
#     # See this excellent trick by Tauana:
#     # https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
#     # Adapted this code to also plot fossil tip uncertainty in red
#
#     phy@data$age_0.95_HPD <- lapply(phy@data$age_0.95_HPD, function(z) {
#       if (is.na(z)) { return(c(NA,NA)) } else { return(z) }
#     })
#
#     minmax <- t(matrix(unlist(phy@data$age_0.95_HPD), nrow = 2))
#     bar_df <- data.frame(node_id=as.integer(phy@data$node), as.data.frame(minmax))
#     names(bar_df) <- c("node_id", "min", "max")
#     if (!is.null(fossils)) {
#       fossil_df <-  dplyr::filter(bar_df, node_id %in% match(fossils,phy@phylo$tip.label))
#     }
#     node_df <- dplyr::filter(bar_df, node_id > ape::Ntip(phy@phylo))
#
#     # plot age densities
#     node_df <- dplyr::left_join(node_df, pp$data, by=c("node_id"="node"))
#     node_df <- dplyr::select(node_df,  node_id, min, max, y)
#     pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y),
#                                      data=node_df, color="blue", size=1.5, alpha=0.3)
#     if (!is.null(fossils)) {
#       fossil_df <- dplyr::left_join(fossil_df, pp$data, by=c("node_id"="node"))
#       fossil_df <- dplyr::select(fossil_df, node_id, min, max, y)
#       pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y),
#                                        data=fossil_df, color="red", size=1.5, alpha=0.4)
#     }
#   }
#
#   # add timeline
#   if (timeline) {
#     # get dimensions
#     n_nodes <- treeio::Nnode(phy)
#     max_age <- max(ape::branching.times(phy@phylo))
#     dx <- max_age %% 10
#
#     # set coordinates
#     ### fix the xlim and ylims - if no error bars, should be a function of max age and n nodes, respectively
#     ### if error bars, -x lim should be as old as the max of the error bar
#     #pp <- pp + ggplot2::coord_cartesian(xlim = c(-max_age,30), ylim=c(-7, n_nodes+1.5), expand=F)
#     pp <- pp + ggplot2::coord_cartesian()
#     pp <- pp + ggplot2::scale_x_continuous(breaks = seq(-max_age-dx,0,10), labels = rev(seq(0,max_age+dx,10)))
#     pp <- pp + ggtree::theme_tree2()
#     pp <- pp + ggplot2::labs(x="Age (Ma)")
#     pp <- pp + ggplot2::theme(legend.position=c(.05, .955), axis.line = ggplot2::element_line(colour = "black"))
#     pp <- ggtree::revts(pp)
#     pp <- .add_epoch_times(pp, max_age, dy_bars=-7, dy_text=-3)
#   }
#
#   # add node labels (PP)
#   if (node_labels) {
#
#     # format posterior data
#     phy@data$posterior[ phy@data$posterior == 1 ] <- NA
#
#     # plot clade support
#     pp$data$posterior_class = NA
#     pp$data$posterior_class[ which(pp$data$posterior>=0.99) ] = ">0.99"
#     pp$data$posterior_class[ which(pp$data$posterior<0.99&pp$data$posterior>0.95) ] = ">0.95"
#     pp$data$posterior_class[ which(pp$data$posterior<=0.95&pp$data$posterior>0.75) ] = ">0.75"
#     pp$data$posterior_class[ which(pp$data$posterior<=0.75&pp$data$posterior>0.5) ] = ">0.50"
#     pp$data$posterior_class[ which(pp$data$posterior<=0.5) ] = "<0.50"
#     pp$data$posterior_class = factor(pp$data$posterior_class, levels=c(">0.99",">0.95",">0.75", ">0.50", "<0.50"))
#     pp$data$label = sapply(pp$data$label, function(x) { gsub("_"," ",x) })
#     pp$data$label = sapply(pp$data$label, function(x) { gsub("subsp. ","",x) })
#
#     # set posterior support by color if no branch coloring
#     if(is.null(color_branch)) {
#       pp <- pp + ggtree::geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], size=2, color="black")
#       pp <- pp + ggtree::geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], ggplot2::aes(color=posterior_class), size=1.5)
#       #### add custom colors and default colors here ####
#       col_post <- c("#000000","#666666","#999999","#BBBBBB","#EEEEEE")
#       names(col_post) <- levels(pp$data$posterior_class)
#       pp <- pp + ggplot2::scale_color_manual(values=col_post, name="Posterior")
#     } else { # if branches are colored, use size for pp
#       pp <- pp +
#         ggtree::geom_nodepoint(data=pp$data[ !is.na(pp$data$posterior), ],
#                                ggplot2::aes(size=posterior_class), color="black") +
#         ggplot2::scale_size_discrete(name="Posterior")
#     }
#
#   }
#
#   # add branch coloring
#   if (!is.null(color_branch)){
#     col_num <- which(colnames(phy@data) == color_branch)
#     phy@data[,col_num] <- as.numeric(as.data.frame(phy@data)[,col_num]) #convert data to numeric
#     name <- .simpleCap(sub(pattern = "_", replacement = " ", color_branch))
#     pp <- pp +
#       ggplot2::aes(color=I(as.data.frame(phy@data)[,col_num])) +
#       ggplot2::scale_color_gradient(low = "#005ac8", high = "#fa7850",
#                                     name = name)
#   }
#
#   # add tip labels
#   if (tip_labels) {
#     # set tip & clade names
#     pp <- pp + ggtree::geom_tiplab(size=2.5, offset=0.2, color = "black")
#     # add extra space on plot for tip labels
#     tree_height <- max(phytools::nodeHeights(phy@phylo))
#     pp <- pp + ggtree::xlim(0, tree_height + tree_height/2)
#   }
#
#   return(pp)
# }
#
