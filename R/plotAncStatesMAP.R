#' plot Anc States MAP
#'
#' Plots the MAP estimates of ancestral states. Can accomodate cladogenetic reconstructions
#' by plotting on shoulders. Defaults to varying the symbols by color to indicate estimated
#' ancestral state and varying the size of the symbol to indicate the posterior probability
#' of that estimate, but symbol shape may also vary to accomodate black and white figures.
#' For more details on the aesthetics options, see parameter details below. For data with
#' many character states (such as choromosome counts), vary the size of the symbol by estimated
#' ancestral state, and vary the posterior probability of that esimate by a color gradient.
#' Text labels at nodes and tips are also available.
#'
#' @param t (treedata object; none) Output of processAncStates() function
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
#' @param tree_layout (character; "rectangular") Tree shape layout, passed to ggtree(). Options
#' are 'rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle', or 'daylight'
#'
#' @param timeline (logical; FALSE) Plot tree with labeled x-axis with timescale in MYA.
#'
#' @examples
#'
#' \dontrun{
#' # Standard ancestral state reconstruction example with various aesthetics
#'
#' # process file
#' file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
#' example <- processAncStates(file, state_labels = c("1" = "Awesome", "2" = "Beautiful", "3" = "Cool!"))
#'
#' # have states vary by color and indicate state pp with size (default)
#' plotAncStatesMAP(t = example, tip_labels_italics = T)
#'
#' # have states vary by color and indicate state pp with size , and add a timeline
#' plotAncStatesMAP(t = example, tip_labels_italics = T, timeline = T)
#'
#' # have states vary by color and symbol, label nodes with pp of states
#' plotAncStatesMAP(t = example,  node_shape_as = "state", node_size = 4, node_shape = c(15, 17,20),
#'         node_size_as = NULL, node_labels_as = "state_posterior")
#'
#' # black and white figure - state as symbols and state pp with text
#' plotAncStatesMAP(t = example, node_color_as = NULL, node_shape_as = "state", node_shape =  c(15, 17,20),
#'         node_size_as = NULL, node_size = 4, node_labels_as = "state_posterior",
#'         node_color = "grey", state_transparency = 1)
#'
#' # default with circular tree
#' plotAncStatesMAP(t = example, tree_layout = "circular")
#'
#' # Chromosome evolution example
#'
#' # process file
#' file <- system.file("extdata", "chromo/ChromEvol_simple_final.tree", package="RevGadgets")
#' chromo_example <- processAncStates(file, labels_as_numbers = TRUE)
#'
#' # plot
#' plotAncStatesMAP(t = chromo_example, node_color_as = "state_posterior", node_size_as = "state",
#'         node_color = RevGadgets:::.colFun(2), tip_labels_offset = 0.005,
#'         node_labels_as = "state", node_labels_offset = 0, tip_labels_states = TRUE,
#'         tip_labels_states_offset = 0, tip_states = FALSE)
#'
#' # DEC example (cladogenetic)
#'
#' # process file
#' file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
#' labs <- c("1" = "K", "2" = "O", "3" = "M", "4" = "H", "5" = "KO",
#' "6" = "KM", "7" = "OM", "8" = "KH", "9" = "OH", "10" = "MH", "11" = "KOM",
#' "12" = "KOH", "13" = "KMH", "14" = "OMH", "15" = "KOMH")
#' dec_example <- processAncStates(file, state_labels = labs)
#'
#' # plot
#' plotAncStatesMAP(t = dec_example, cladogenetic = TRUE, tip_labels_offset = 0.5)
#' }
#'
#' @export


plotAncStatesMAP <- function(t,
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
                    tree_layout = "rectangular",
                    timeline = FALSE) {
    ##### parameter compatability checks! #####
  if (class(t) != "treedata") stop("t should be a treedata objects")
  if (is.logical(cladogenetic) == FALSE) stop("cladogenetic should be TRUE or FALSE")
  if (is.logical(tip_labels) == FALSE) stop("tip_labels should be TRUE or FALSE")
  if (is.numeric(tip_labels_size) == FALSE) stop("tip_labels_size should be a number")
  if (is.numeric(tip_labels_offset) == FALSE) stop("tip_labels_offset should be a number")
  if (is.logical(tip_labels_italics) == FALSE) stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_remove_underscore) == FALSE) stop("tip_labels_remove_underscore should be TRUE or FALSE")
  if (is.logical(tip_labels_states) == FALSE) stop("tip_labels_states should be TRUE or FALSE")
  if (is.numeric(tip_labels_states_size) == FALSE) stop("tip_labels_states_size should be a number")
  if (is.numeric(tip_labels_states_offset) == FALSE) stop("tip_labels_states_offsetshould be a number")
  if (is.null(node_labels_as) == FALSE) {
    node_labels_as <- match.arg(node_labels_as, choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.numeric(node_labels_size) == FALSE) stop("node_labels_size should be a number")
  if (is.numeric(node_labels_offset) == FALSE) stop("node_labels_offset should be a number")
  if (is.null(node_size_as) == FALSE) {
    node_size_as <- match.arg(node_size_as, choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_color_as) == FALSE) {
    node_color_as <- match.arg(node_color_as, choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as != "state") stop("node_shape_as should be NULL or 'state'")
  }
  if (is.numeric(node_shape) == FALSE) stop("node_shape should be a number indicating symbol type")
  if (is.character(node_color) == FALSE) stop ("node_color should be 'default' or valid color(s)")
  if (node_color[1] != "default" & any(.isColor(node_color) == FALSE) ) stop("node_color should be valid color(s)")
  if (any(is.numeric(node_size) == FALSE)) stop("node_size should be a single number or a vector of two numbers")
  if (length(node_size) > 2) stop("node_size should be a single number or a vector of two numbers")
  if (is.logical(tip_states) == FALSE) stop("tip_states should be TRUE or FALSE")
  if (is.numeric(tip_states_size) == FALSE) stop("tip_states_size should be a number")
  if (is.numeric(tip_states_shape) == FALSE) stop("tip_states_shape should be a number indicating symbol type")
  if (is.numeric(state_transparency) == FALSE) stop("state_transparency should be a number between 0 - 1")
  if (state_transparency > 1 | state_transparency < 0) stop("state_transparency should be a number between 0 - 1")
  tree_layout <- match.arg(tree_layout, choices = c('rectangular', 'slanted', 'fan', 'circular', 'radial', 'equal_angle','daylight'))
  if (is.logical(timeline) == FALSE) stop("timeline should be TRUE or FALSE")
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) { stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")}
  }

  ##### calculate helper variables #####

  tree <- attributes(t)$phylo
  n_node <- ggtree:::getNodeNum(tree)

  ##### create basic tree plot #####
  # supressing warning:
  # Warning message:
  #   `tbl_df()` is deprecated as of dplyr 1.0.0.
  # Please use `tibble::as_tibble()` instead.
  p <- suppressWarnings(ggtree::ggtree(t, layout = tree_layout, ladderize = TRUE))

  ##### process column names #####

  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE & "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE & "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }

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
  if (!is.null(node_color_as) && node_color_as == "states") {
    if (node_color[1] == "default" & length(all_states) > 12) {
      stop(paste0(length(all_states), " states in dataset; please provide colors (default only can provide up to 12"))
    }

    # check if number of states not equal to provided colors
    if (node_color[1] != "default" & length(node_color) < length(all_states)) {
      stop(paste0("You provided fewer colors in node_color than states in your dataset. There are ",
                  length(all_states), " states and you provide ", length(node_color), " colors."))
    }
    if (node_color[1] != "default" & length(node_color) > length(all_states)) {
      stop(paste0("You provided more colors in node_color than states in your dataset. There are ",
                  length(all_states), " states and you provide ", length(node_color), " colors."))
    }
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
        p <- p + ggtree::geom_tippoint(ggplot2::aes(colour = node_color_as),
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
        p <- p + ggtree::geom_tippoint(ggplot2::aes(shape = node_shape_as),
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
        p <-  p + ggtree::geom_tippoint(ggplot2::aes(shape = node_shape_as,
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
        p <- p + ggtree::geom_tippoint(ggplot2::aes(size = node_size_as),
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
          p <- p + ggtree::geom_tippoint(ggplot2::aes(size = node_size_as,
                                                     color = node_color_as),
                                         shape = tip_states_shape, alpha = state_transparency)
        }
    }
  }

  # plot symbols at nodes and shoulders
  blank_nodes <- is.null(node_color_as) == TRUE & is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE
  if (blank_nodes == FALSE) {
    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = node_color_as,
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = node_color_as, size = node_size_as),
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = node_color_as, shape = node_shape_as),
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(shape = node_shape_as, size = node_size_as),
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(colour = node_color_as), size = node_size,
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(size = node_size_as), color = colors,
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
      p <- p + ggtree::geom_nodepoint(ggplot2::aes(shape = node_shape_as), color = colors,
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
  }

  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = end_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size) +
                 ggtree::geom_text(ggplot2::aes(label = start_state_1, x = x_anc, y = y),
                                   hjust = "right", nudge_x = node_labels_offset, size = node_labels_size, na.rm = TRUE)
      } else if (cladogenetic == FALSE & state_pos_str_base == "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = anc_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      } else if (cladogenetic == FALSE & state_pos_str_base != "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = end_state_1), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = .convertAndRound(end_state_1_pp) ), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size) +
          ggtree::geom_text(ggplot2::aes(label = .convertAndRound(start_state_1_pp), x = x_anc, y = y),
                            hjust = "right", nudge_x = node_labels_offset, size = node_labels_size, na.rm = TRUE)
      } else if (cladogenetic == FALSE & state_pos_str_base == "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = .convertAndRound(anc_state_1_pp) ), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      } else if (cladogenetic == FALSE & state_pos_str_base != "anc_state_") {
        p <- p + ggtree::geom_nodelab(ggplot2::aes(label = .convertAndRound(end_state_1_pp) ), hjust="left",
                                      nudge_x = node_labels_offset, size = node_labels_size)
      }
    } else if (node_labels_as == "node_posterior") {
      p <- p + ggtree::geom_nodelab(ggplot2::aes(label = .convertAndRound(posterior) ), hjust = "left",
                                    nudge_x = node_labels_offset, size = node_labels_size)
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

  # add timeline
  if (timeline == TRUE) {

    if ("age_0.95_HPD" %in% colnames(p$data)){
      p$data$age_0.95_HPD <- lapply(p$data$age_0.95_HPD, function(z) {
        if (is.null(z) || is.na(z)) { return(c(NA,NA)) } else { return(as.numeric(z)) }
      })
      minmax <- t(matrix(unlist(p$data$age_0.95_HPD), nrow = 2))
    }
      if (!"age_0.95_HPD" %in% colnames(p$data)) {
        max_age <- max(ape::branching.times(tree))
        print("Adding timeline using branch lengths. Make sure branch lengths are in MYA units!")
      } else {
        max_age <- max(minmax, na.rm =TRUE)
        print("Adding timeline using age_0.95_HPD from tree file.")
      }

    if (max_age > 100){
      interval <- 50
    } else {interval <- 10}

    dx <- max_age %% interval

    # set coordinates
    p <- p + ggplot2::coord_cartesian()
    p <- p + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                         expand = c(0, 0),
                                         limits = c(-max_age*1.05, tree_height/2),
                                         breaks = -rev(seq(0,max_age+dx,interval)),
                                         labels = rev(seq(0,max_age+dx,interval)),
    )
    p <- p + ggtree::theme_tree2()
    #pp <- p + ggplot2::theme(legend.position=c(.05, .955), axis.line = ggplot2::element_line(colour = "black"))
    p <- ggtree::revts(p)
    n_tips <- length(tree$tip.label)
    p <- p + ggplot2::scale_y_continuous(limits = c(-n_tips/20, n_tips*1.1), expand = c(0, 0))
    #p  <- p + ggplot2::annotate(geom = "segment", x = 0, xend = -max_age, y = -n_tips/20, yend = -n_tips/20)
    p <- .add_epoch_times(p, max_age, dy_bars=-n_tips/20, dy_text=-n_tips/25)
  }

  return(p)
}
