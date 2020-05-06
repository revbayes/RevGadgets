#' plot MAP
#'
#' Plot character states and posterior probabilities as points on nodes (size = pp, colour = state)
#' and character state on tips (different shape from nodes, colour = state)
#'
#' @param t (treedata object; none) Output of processAncStatesDiscrete() function
#' containing tree and ancestral states.
#' @param cladogenetic (logical; FALSE) Plot shoulder states of cladogenetic analyses?
#' @param tip_labels (logical; TRUE) Label taxa at tips?
#' @param tip_labels_size (numeric; 2) Size of tip labels.
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from tree.
#' @param node_labels_as (character; "NULL") Optional plotting of text at nodes. Possible
#' values are "state" for the ancestral states , "state_posterior" for posterior probabilities
#' of the estimated ancestral state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL for not plotting any text at the nodes (default).
#' @param node_labels_size (numeric; 2) Size of node labels text. Ignored if node_labels_as = NULL.
#' @param node_labels_offset (numeric; 0.1) Horizontal offset of node labels from nodes.
#' Ignored if node_labels_as = NULL.
#' @param node_size_as (character; "state_posterior") How to vary size of node symbols. Options
#' are "state_posterior" (default) for posterior probabilities of the estimated ancestral state,
#' "node_posterior" or the posterior probability of the node on the tree, "state" for vary size by the
#' ancestral state itself in cases where there are many character states (e.g. chromosome numbers;
#' we do not recommend this option for characters with few states), or NULL for fixed symbol size.
#' @param node_color_as (character; "state") How to vary to color of node symbols. Options are
#' "state" (default) to vary by estimated ancestral states, "state_posterior" for posterior probabilities
#' of the estimated ancestral state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL to set all as one color.
#' @param node_shape_as (character; NULL) Option to vary node symbol by shape. Options are NULL
#' to keep shape constant or "state" to vary shape by ancestral state.
#' @param node_shape (integer; 19) Shape type for nodes. If node_shape_as = "state", provide a vector
#' with length of the number of states.
#' @param node_color ("character"; "default") Colors for node symbols. Defaults to default RevGadgets
#' colors. If node_color_as = "state', provide a vector of length of the character states.
#' If node_color_as = "posterior", provide a vector of length 2 to generate a color gradient.
#' @param node_size (numeric; c(2, 6)) Range of sizes, or fixed size, for node symbols.
#' If node_size_as = "state_posterior", "node_posterior", or "state", numeric vector of length two.
#' If node_size_as = NULL, numeric vector of length one.
#' @param tip_states (logical; TRUE) Plot states of taxa at tips?
#' @param tip_states_size (numeric; node_size) Size for tip symbols. Defaults to the same
#' size as node symbols.
#' @param tip_states_shape (integer; node_shape) Shape for tip symbols. Defaults to the same
#' as node symbols.
#' @param state_transparency (integer; 0.75) Alpha (transparency) of state symbols- varies from
#' 0 to 1.
#' @param show_state_legend (logical; TRUE) Plot legend for states?
#' @param show_posterior_legend (logical; TRUE) Plot legend for posterior probabilities?
#' @param tree_layout (character; "rectangular") Tree shape layout, passed to ggtree(). Options
#' are 'rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle', or 'daylight'
#'
#' @export


plotMAP <- function(t,
                    # option for plotting shoulder states
                    cladogenetic = FALSE,

                    # label taxa at tips
                    tip_labels = TRUE,
                    tip_labels_size = 2,
                    tip_labels_offset = 1,

                    # text labels at nodes
                    node_labels_as = NULL,
                    node_labels_size = 2,
                    node_labels_offset = 0.1,

                    # what to plot at nodes
                    node_size_as = "state_posterior",
                    node_color_as = "state",
                    node_shape_as = NULL,

                    # aesthetics for plotting at nodes
                    node_shape = 19,
                    node_color = "default",
                    node_size = c(2, 6),

                    # aesthetics for tip states (inherents additional aesthetics from nodes)
                    tip_states = TRUE,
                    tip_states_size = node_size,
                    tip_states_shape = node_shape,

                    state_transparency = 0.75,
                    show_state_legend = TRUE,
                    show_posterior_legend = TRUE,
                    tree_layout = "rectangular") {

  # add in parameter compatability checks!


  # get number of nodes
  tree <- attributes(t)$phylo
  n_node <- ggtree:::getNodeNum(tree)

  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)

  # set default color palette
  if (node_color == "default") {
    if (is.null(node_color_as) == TRUE) {
      colors <- .colFun(1)
    } else if (node_color_as == "state") {
      nstates <- length(unique(p$data$anc_state_1))
      colors <- .colFun(nstates)
      names(colors) <- levels(as.factor(p$data$anc_state_1))
    } else if (node_color_as == "node_posterior" |
               node_color_as == "state_posterior") {
      colors <- .colFun(2)
    }
  } else {
    colors <- node_color
  }

  # check aesthetics lengths and adjust if needed
  # shape
  if (is.null(node_shape_as) == TRUE) {
    if (length(node_shape) > 1) {
      node_shape <- node_shape[1]
    }
  }
  # color
  if (is.null(node_color_as) == TRUE) {
    if (length(colors) > 1) {
      colors <- colors[1]
    }
  }
  # size
  if (is.null(node_size_as) == TRUE) {
    if (length(node_size) > 1) {
      node_size <- node_size[1]
    }
  }

  # add tip labels
  if (tip_labels == TRUE) {
    p <- p + ggtree:::geom_tiplab(size = tip_labels_size, offset = tip_labels_offset)
  }

  # add the tip states
  if (tip_states == TRUE) {

    # unless node size should vary by state, don't allow tip sizes to vary
    if (is.null(node_size_as) == TRUE || node_size_as != "state") {
      tip_states_size <- tip_states_size[1]
    }

    # vary tip symbols by color only
    # when shape is null and size is not state
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state" &
          is.null(node_shape_as) == TRUE &
          (is.null(node_size_as) == TRUE || node_size_as != "state")) {
        p <- p + ggtree::geom_tippoint(ggtree::aes(colour = factor(anc_state_1)),
                                       size = tip_states_size, alpha = state_transparency,
                                       shape = tip_states_shape)
      }
    }

    # vary tip symbols by shape only
    # when shape is state, color is not state, size is not state
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state" &
          (is.null(node_color_as) == TRUE || node_color_as != "state") &
          (is.null(node_size_as) == TRUE || node_size_as != "state")) {
        p <- p + ggtree::geom_tippoint(ggtree::aes(shape = factor(anc_state_1)),
                                       size = tip_states_size, alpha = state_transparency,
                                       color = colors)
      }
    }

    # vary tip symbol by shape and color
    # when shape is state, color is state, and size is anything but state
    if (is.null(node_color_as) == FALSE & is.null(node_shape_as) == FALSE) {
      if (node_color_as == "state" &
          node_shape_as == "state" &
          (is.null(node_size_as) == TRUE || node_size_as != "state")) {
        p <-  p + ggtree::geom_tippoint(ggtree::aes(shape = factor(anc_state_1),
                                                    color = factor(anc_state_1)),
                                        size = tip_states_size, alpha = state_transparency)
      }
    }

    # vary tip symbol by size only
    # when size is state, color is not state, and shape is null
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state" &
          (is.null(node_color_as) == TRUE || node_color_as != "state") &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(ggtree::aes(size = as.numeric(anc_state_1)),
                                       shape = tip_states_shape, alpha = state_transparency,
                                       color = colors)
      }
    }

    # vary tip symbol by size and color
    # when size is state, color is state, and shape is null
    if (is.null(node_size_as) == FALSE & is.null(node_color_as) == FALSE) {
      if (node_size_as == "state" &
          node_color_as == "state" &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(ggtree::aes(size = anc_state_1,
                                                   color = anc_state_1),
                                       shape = tip_states_shape, alpha = state_transparency)
      }
    }

  } # end tip_states == TRUE

  # set up ancestral states dataframe if anc_states not already in phylo data object
  # Carrie note: when would this happen? this should be taken care of by process anc states....
  #if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
  #  anc_data <- data.frame(node = names(attributes(t)$data$end_state_1),
  #                        anc_state_1 = levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
  #                        anc_state_1_pp = as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
  #  `%<+%` <- ggtree::`%<+%`
  #  p <- p %<+% anc_data
  #}
  #

  # plot symbols at nodes and tips
  blank_nodes <- is.null(node_color_as) == TRUE & is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE

  if (blank_nodes == FALSE) {

    #translate to column names
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state") {p$data$node_color_as <- factor(p$data$anc_state_1)}
      if (node_color_as == "node_posterior") {p$data$node_color_as <- as.numeric(p$data$posterior)}
      if (node_color_as == "state_posterior") {p$data$node_color_as <- as.numeric(p$data$anc_state_1_pp)}
    }

    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state") {p$data$node_size_as <- p$data$anc_state_1}
      if (node_size_as == "node_posterior") {p$data$node_size_as <- as.numeric(p$data$posterior)}
      if (node_size_as == "state_posterior") {p$data$node_size_as <- as.numeric(p$data$anc_state_1_pp)}
    }

    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state") {p$data$node_shape_as <- factor(p$data$anc_state_1)}
      if (node_shape_as == "node_posterior") {p$data$node_shape_as <- as.numeric(p$data$posterior)}
      if (node_shape_as == "state_posterior") {p$data$node_shape_as <- as.numeric(p$data$anc_state_1_pp)}
    }

    ## This is annoying and will be annoying to debug- is there a way to automate this??

    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as,
                                                  size = node_size_as, shape = node_shape_as),
                                      alpha = state_transparency)
    }

    # plot if color and size vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as, size = node_size_as),
                                      shape = node_shape, alpha = state_transparency)
    }

    #plot if color and shape vary
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as, shape = node_shape_as),
                                      size = node_size, alpha = state_transparency)
    }

    #plot if size and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(shape = node_shape_as, size = node_size_as),
                                      color = colors, alpha = state_transparency)
    }

    #plot if just color varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as), size = node_size,
                                      shape = node_shape, alpha = state_transparency)
    }

    #plot if just size varies
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(size = node_size_as), color = colors,
                                      shape = node_shape, alpha = state_transparency)
    }

    #plot if just shape varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(shape = node_shape_as), color = colors,
                                      size = node_size, alpha = state_transparency)
    }

  } # end blank_nodes == FALSE

  # add cladogenetic events
  if (cladogenetic == TRUE) {
    stop("Start states not yet implemented for MAP ancestral states.")
  }

  # what does this do? - Carrie
  #pp_offset_range <- 2*(c(min(pp), max(pp)) - 0.5)
  #nd_offset_interval <- node_size_range[2] - node_size_range[1]
  #nd_offset <- node_size_range[1]
  #node_size_range <- pp_offset_range * nd_offset_interval + nd_offset

  #node_size_range[1] = node_size_range[1] * min(pp) / 0.5
  #node_size_range[2] = node_size_range[2] * max(pp)

# what does this do?
#  if (node_label_size == 0) {
#    p <- p + ggtree::geom_text(ggtree::aes(label=sprintf("%.02f", as.numeric(anc_state_1_pp))), hjust="left", nudge_x=node_label_nudge_x, size=node_pp_label_size)
#  }

  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      p <- p + ggtree::geom_nodelab(ggtree::aes(label = anc_state_1), hjust="left",
                                    nudge_x = node_labels_offset, size = node_labels_size)
    } else if (node_labels_as == "state_posterior") {
      p <- p + ggtree::geom_nodelab(ggtree::aes(label = .convertAndRound(anc_state_1_pp)), hjust="left",
                                    nudge_x = node_labels_offset, size = node_labels_size)
    }
  }

  # add custom colors, shapes, and sizes

  if (is.null(node_size_as) == FALSE) {
    p <- p + ggplot2::scale_size(range = node_size, name = node_size_as)
  }
  if (is.null(node_color_as) == FALSE) {
    p <- p + ggplot2::scale_color_manual(values = colors, name = node_color_as)
  }
  if (is.null(node_shape_as) == FALSE) {
    p <- p + ggplot2::scale_shape_manual(values = node_shape, name = node_shape_as)
  }

  # add space on x axis for tip labels

  if (tip_labels == TRUE) {
    tree_height <- max(phytools::nodeHeights(t@phylo))
    p <- p + ggtree::xlim(0, tree_height + tree_height/2)
  }
  # set up the legend
  #if (show_state_legend == TRUE) {
  #  if ()
  #  p = p + ggplot2::guides(colour=ggplot2::guide_legend("State"), order=1)
  #} else {
  #  p = p + ggplot2::guides(colour=FALSE, order=2)
  #}
  #if (show_posterior_legend == TRUE) {
  #  p = p + ggplot2::guides(size=ggplot2::guide_legend("Posterior Probability"), order=3)
  #} else {
  #  p = p + ggplot2::guides(size=FALSE, order=4)
  #}
  return(p)
}
