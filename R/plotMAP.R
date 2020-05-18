#' plot MAP
#'
#' Plot character states and posterior probabilities as points on nodes (size = pp, colour = state)
#' and character state on tips (different shape from nodes, colour = state)
#'
#' @param t (treedata object; none) Output of processAncStatesDiscrete() function
#' containing tree and ancestral states.
#' @param cladogenetic (logical; FALSE) Plot shoulder states of cladogenetic analyses?
#' @param tip_labels (logical; TRUE) Label taxa labels at tips?
#' @param tip_labels_size (numeric; 2) Size of tip labels.
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from tree.
#' @param tip_labels_italics (logical; FALSE) Italicize tip labels?
#' @param tip_labels_remove_underscore (logical; TRUE) Remove underscores from tip labels?
#' @param tip_labels_states (logical; FALSE) Optional plotting of text at tips in addition
#' to taxa labels.
#' @param tip_labels_states_size (numeric; 2) Size of state labels at tips. Ignored if
#' tip_labels_states is FALSE.
#' @param tip_labels_states_offset (numeric; 0.1) Horizontal offset of tip state labels.
#' Ignored if tip_labels_states = NULL.
#' @param node_labels_as (character; NULL) Optional plotting of text at nodes. Possible
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
                    tip_labels_italics = FALSE,
                    tip_labels_remove_underscore = TRUE,

                    # label states at tips
                    tip_labels_states = FALSE,
                    tip_labels_states_size = 2,
                    tip_labels_states_offset = 0.1,

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
  ##### parameter compatability checks! #####

  ##### calculate helper variables #####

  tree <- attributes(t)$phylo
  n_node <- ggtree:::getNodeNum(tree)

  # get names of data columns
  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE & "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "start_state_"
  } else if (cladogenetic == FALSE & "anc_state_" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }

  ##### create basic tree plot #####
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)

  ##### translate to column names #####
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {p$data$node_color_as <- factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))}
    if (node_color_as == "node_posterior") {p$data$node_color_as <- as.numeric(p$data$posterior)}
    if (node_color_as == "state_posterior") {p$data$node_color_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1", "_pp")))}
  }

  if (is.null(node_size_as) == FALSE) {
    if (node_size_as == "state") {p$data$node_size_as <- as.integer(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))}
    if (node_size_as == "node_posterior") {p$data$node_size_as <- as.numeric(p$data$posterior)}
    if (node_size_as == "state_posterior") {p$data$node_size_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1", "_pp")))}
  }

  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as == "state") {p$data$node_shape_as <- factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))}
    if (node_shape_as == "node_posterior") {p$data$node_shape_as <- as.numeric(p$data$posterior)}
    if (node_shape_as == "state_posterior") {p$data$node_shape_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1", "_pp")))}
  }

  if (cladogenetic == TRUE) {
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state") {p$data$clado_node_color_as <- factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))}
      if (node_color_as == "node_posterior") {p$data$clado_node_color_as <- 1}
      if (node_color_as == "state_posterior") {p$data$clado_node_color_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1", "_pp")))}
    }

    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state") {p$data$clado_node_size_as <- as.integer(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))}
      if (node_size_as == "node_posterior") {p$data$clado_node_size_as <- 1}
      if (node_size_as == "state_posterior") {p$data$clado_node_size_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1", "_pp")))}
    }

    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state") {p$data$clado_node_shape_as <- factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))}
      if (node_shape_as == "node_posterior") {p$data$clado_node_shape_as <- 1}
      if (node_shape_as == "state_posterior") {p$data$clado_node_shape_as <- as.numeric(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1", "_pp")))}
    }
  }

  # gather list of all character states from data
  if (cladogenetic == TRUE) {
    all_states <- unique(c(p$data$start_state_1, p$data$end_state_1))
  } else {
    all_states <- unique(factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))))
  }

  ##### color processing and checks #####
  # check if number of states exceeds default color palette options
  if (node_color == "default" & length(all_states) > 12) {
    stop(paste0(length(all_states), " states in dataset; please provide colors (default only can provide up to 12"))
  }

  # check if number of states not equal to provided colors
  if (node_color != "default" & length(node_color) < length(all_states)) {
    stop(paste0("You provided fewer colors in node_color than states in your dataset. There are ",
                length(all_states), " states and you provide ", length(node_color), " colors."))
  }
  if (node_color != "default" & length(node_color) > length(all_states)) {
    stop(paste0("You provided more colors in node_color than states in your dataset. There are ",
                length(all_states), " states and you provide ", length(node_color), " colors."))
  }

  # set default colors
  if (any(node_color == "default")) {
    if (is.null(node_color_as) == TRUE) {
      colors <- .colFun(1)
    } else if (node_color_as == "state") {
      nstates <- length(all_states)
      colors <- .colFun(nstates)
      names(colors) <- all_states
    } else if (node_color_as == "node_posterior" |
               node_color_as == "state_posterior") {
      colors <- .colFun(2)
    }
  } else {
    colors <- node_color
  }

  ##### adjust aesthetics lengths if needed #####
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

  ##### reformat labels if necessary #####
  if (tip_labels_remove_underscore) {p$data$label <- gsub("_", " ", p$data$label)}

  ##### calculate cladogenetic plotting data #####
  if (cladogenetic == TRUE) {
    x <- .getXcoord(tree)
    y <- .getYcoord(tree)
    x_anc <- numeric(n_node)
    node_index <- numeric(n_node)
    for (i in 1:n_node) {
      if (.getParent(tree, i) != 0) {
        # if not the root, get the x coordinate for the parent node
        x_anc[i] <- x[.getParent(tree, i)]
        node_index[i] <- i
      }
    }
    shoulder_data <- data.frame(node = node_index, x_anc = x_anc, y = y)
    `%<+%` <- ggtree::`%<+%`
    p <- p %<+% shoulder_data
  }
  ##### start plotting #####

  # add tip labels
  if (tip_labels == TRUE) {
    if (tip_labels_italics == TRUE) {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = paste0('italic(`', label, '`)')),
          size = tip_labels_size,
          offset = tip_labels_offset,
          parse = TRUE
        )
    } else {
      p <-
        p + ggtree::geom_tiplab(size = tip_labels_size, offset = tip_labels_offset)
    }
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
          (is.null(node_size_as) == TRUE || node_size_as != "state"))  {
        p <- p + ggtree::geom_tippoint(ggtree::aes(colour = node_color_as),
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
        p <- p + ggtree::geom_tippoint(ggtree::aes(shape = node_shape_as),
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
        p <-  p + ggtree::geom_tippoint(ggtree::aes(shape = node_shape_as,
                                                    color = node_color_as),
                                        size = tip_states_size, alpha = state_transparency)
      }
    }

    # vary tip symbol by size only
    # when size is state, color is not state, and shape is null
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state" &
          (is.null(node_color_as) == TRUE || node_color_as != "state") &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(ggtree::aes(size = node_size_as),
                                       shape = tip_states_shape, alpha = state_transparency,
                                       color = "grey")
      }
    }

      # vary tip symbol by size and color
      # when size is state, color is state, and shape is null
      if (is.null(node_size_as) == FALSE & is.null(node_color_as) == FALSE) {
        if (node_size_as == "state" &
            node_color_as == "state" &
            is.null(node_shape_as) == TRUE) {
          p <- p + ggtree::geom_tippoint(ggtree::aes(size = node_size_as,
                                                     color = node_color_as),
                                         shape = tip_states_shape, alpha = state_transparency)
        }
    }
  }

  # plot symbols at nodes and tips
  blank_nodes <- is.null(node_color_as) == TRUE & is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE
  if (blank_nodes == FALSE) {
    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as,
                                                  size = node_size_as, shape = node_shape_as),
                                      alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = clado_node_color_as,
                                                     size = clado_node_size_as,
                                                     shape = clado_node_shape_as,
                                                     x = x_anc, y = y),
                                        na.rm = TRUE, alpha = state_transparency) +
          ggtree::geom_tippoint(ggplot2::aes(colour = clado_node_color_as,
                                             size = clado_node_size_as,
                                             shape = clado_node_shape_as,
                                             x = x_anc, y = y),
                                na.rm=TRUE, alpha = state_transparency)
      }
    }


    # plot if color and size vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as, size = node_size_as),
                                      shape = node_shape, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = clado_node_color_as,
                                                     size = clado_node_size_as,
                                                     x = x_anc, y = y),
                                        shape = node_shape, na.rm = TRUE,
                                        alpha = state_transparency) +
                 ggtree::geom_tippoint(ggplot2::aes(colour = clado_node_color_as,
                                                    size = clado_node_size_as,
                                                    x = x_anc, y = y),
                                       shape = node_shape, na.rm=TRUE,
                                       alpha = state_transparency)
      }
    }

    #plot if color and shape vary
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as, shape = node_shape_as),
                                      size = node_size, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = clado_node_color_as,
                                                     shape = clado_node_shape_as,
                                                     x = x_anc, y = y),
                                        size = node_size, na.rm = TRUE,
                                        alpha = state_transparency) +
                 ggtree::geom_tippoint(ggplot2::aes(colour = clado_node_color_as,
                                                    shape = clado_node_shape_as,
                                                    x = x_anc, y = y),
                                       size = node_size, na.rm=TRUE,
                                       alpha = state_transparency)
      }
    }

    #plot if size and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(shape = node_shape_as, size = node_size_as),
                                      color = colors, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(size = clado_node_size_as,
                                                     shape = clado_node_shape_as,
                                                     x = x_anc, y = y),
                                        color = colors, na.rm = TRUE,
                                        alpha = state_transparency) +
                 ggtree::geom_tippoint(ggplot2::aes(size = clado_node_size_as,
                                                    shape = clado_node_shape_as,
                                                    x = x_anc, y = y),
                                       color = colors, na.rm=TRUE,
                                       alpha = state_transparency)
      }
    }

    #plot if just color varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour = node_color_as), size = node_size,
                                      shape = node_shape, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(color = clado_node_color_as,
                                                     x = x_anc, y = y),
                                        size = node_size, shape = node_shape, na.rm = TRUE,
                                        alpha = state_transparency) +
                 ggtree::geom_tippoint(ggplot2::aes(color = clado_node_color_as,
                                                    x = x_anc, y = y),
                                       size = node_size, shape = node_shape, na.rm = TRUE,
                                       alpha = state_transparency)
      }
    }

    #plot if just size varies
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == TRUE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(size = node_size_as), color = colors,
                                      shape = node_shape, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(size = clado_node_size_as,
                                                     x = x_anc, y = y),
                                        color = colors, shape = node_shape, na.rm = TRUE,
                                        alpha = state_transparency) +
                 ggtree::geom_tippoint(ggplot2::aes(color = clado_node_color_as,
                                                    x = x_anc, y = y),
                                       color = colors, shape = node_shape, na.rm = TRUE,
                                       alpha = state_transparency)
      }
    }

    #plot if just shape varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(shape = node_shape_as), color = colors,
                                      size = node_size, alpha = state_transparency)
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodepoint(ggplot2::aes(shape = clado_node_shape_as,
                                                     x = x_anc, y = y),
                                        color = colors, size = node_size, na.rm = TRUE,
                                        alpha = state_transparency) +
          ggtree::geom_tippoint(ggplot2::aes(shape = clado_node_shape_as,
                                             x = x_anc, y = y),
                                color = colors, size = node_size, na.rm = TRUE,
                                alpha = state_transparency)
      }
    }

  } # end blank_nodes == FALSE

  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = end_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size) +
                 ggtree::geom_text(ggplot2::aes(label = start_state_1, x = x_anc, y = y),
                                   hjust = "right", nudge_x = node_labels_offset, size = node_labels_size, na.rm = TRUE)
      } else if (cladogenetic == FALSE & state_pos_str_base == "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = anc_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      } else if (cladogenetic == FALSE & state_pos_str_base != "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = end_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = end_state_1_pp), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size) +
          ggtree::geom_text(ggplot2::aes(label = start_state_1_pp, x = x_anc, y = y),
                            hjust = "right", nudge_x = node_labels_offset, size = node_labels_size, na.rm = TRUE)
      } else if (cladogenetic == FALSE & state_pos_str_base == "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = anc_state_1_pp), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      } else if (cladogenetic == FALSE & state_pos_str_base != "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggtree::aes(label = end_state_1_pp), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      }
    }
  }

  # add tip states labels (text)
  if (tip_labels_states == TRUE) {
    if (state_pos_str_base == "anc_state_") {
      p <- p + ggtree::geom_tiplab(ggplot2::aes(label = anc_state_1), hjust="left",
                                   offset = tip_labels_states_offset, size = tip_labels_states_size)
    } else {
      p <- p + ggtree::geom_tiplab(ggplot2::aes(label = end_state_1), hjust="left",
                                   offset = tip_labels_states_offset, size = tip_labels_states_size)
    }
  }

  # add custom colors, shapes, and sizes
  if (is.null(node_size_as) == FALSE) {
    p <- p + ggplot2::scale_size(range = node_size, name = node_size_as)
  }
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      p <- p + ggplot2::scale_color_manual(values = colors, name = node_color_as)
    } else if (node_color_as == "state_posterior" | node_color_as == "node_posterior") {
      p <- p + ggplot2::scale_color_gradient(low = colors[1], high = colors[2], name = node_color_as)
    }
  }
  if (is.null(node_shape_as) == FALSE) {
    p <- p + ggplot2::scale_shape_manual(values = node_shape, name = node_shape_as)
  }

  # add space on x axis for tip labels
  if (tip_labels == TRUE) {
    tree_height <- max(phytools::nodeHeights(t@phylo))
    p <- p + ggtree::xlim(0, tree_height + tree_height/2)
  }

  return(p)
}
