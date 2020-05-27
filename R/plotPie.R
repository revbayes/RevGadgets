#' plot Pie
#'
#' Plot character states and posterior probabilities as pies on nodes.
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
#' to taxa labels. Will plot the MAP state label.
#' @param tip_labels_states_size (numeric; 2) Size of state labels at tips. Ignored if
#' tip_labels_states is FALSE.
#' @param tip_labels_states_offset (numeric; 0.1) Horizontal offset of tip state labels.
#' Ignored if tip_labels_states = NULL.
#' @param node_labels_as (character; NULL) Optional plotting of text at nodes. Possible
#' values are "state" for the MAP ancestral states, "state_posterior" for posterior probabilities
#' of the estimated ancestral state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL for not plotting any text at the nodes (default).
#' @param node_labels_size (numeric; 2) Size of node labels text. Ignored if node_labels_as = NULL.
#' @param node_labels_offset (numeric; 0.1) Horizontal offset of node labels from nodes.
#' Ignored if node_labels_as = NULL.
#' @param node_pie_size_as (character; "posterior") How to vary size of node symbols. Options
#' are "node_posterior" or the posterior probability of the node on the tree or NULL for fixed symbol size.
#' @param node_color_as (character; "state") How to vary to color of node symbols. Options are
#' "state" (default) to vary by estimated ancestral states, "state_posterior" for posterior probabilities
#' of the estimated ancestral state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL to set all as one color.
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
#' #'
#' @export

plotPie <- function(t,
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

                    # pies aesthetics
                    pie_colors = "default",
                    node_pie_size = 0.1,
                    shoulder_pie_size = node_pie_size,
                    tip_pies = TRUE,
                    tip_pie_size = 0.05,

                    # nudges to center pies
                    node_pie_nudge_x = 0,
                    node_pie_nudge_y = 0,
                    tip_pie_nudge_x = node_pie_nudge_x,
                    tip_pie_nudge_y = node_pie_nudge_y,
                    shoulder_pie_nudge_x = node_pie_nudge_x,
                    shoulder_pie_nudge_y = node_pie_nudge_y,

                    collapse_states = FALSE,
                    collapse_states_threshold = 0.05,
                    state_transparency = 0.75) {
recover()
    ##### create basic tree plot #####
  p <- ggtree:::ggtree(t, ladderize = TRUE)

  ##### calculate helper variables #####
  tree <- attributes(t)$phylo
  n_node <- ggtree:::getNodeNum(tree)
  n_tip <- length(tree$tip.label)
  node_idx <- (n_tip+1):n_node
  tip_idx <- 1:n_tip
  all_idx <- 1:n_node

  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE & "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE & "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }

  ##### get state labels #####
  if (cladogenetic == TRUE) {
    all_states <- as.character(unique(c(p$data$start_state_1, p$data$end_state_1,
                                        p$data$start_state_2, p$data$end_state_2,
                                        p$data$start_state_3, p$data$end_state_3)))
  } else {
    all_states <- as.character(unique(c(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")),
                                        dplyr::pull(p$data, paste0(state_pos_str_base[1], "2")),
                                        dplyr::pull(p$data, paste0(state_pos_str_base[1], "3")))))
  }

  state_labels <- c(na.omit(all_states), "other")

  # add other so it gets included in legend even if states aren't collapsed
  if ("anc_state_" %in% state_pos_str_base) {p$data$anc_state_other <- "other"}
  if ("end_state_" %in% state_pos_str_base) {p$data$end_state_other <- "other"}

  ##### color processing and checks #####
  if (collapse_states == TRUE) {
    used_states <- .collect_probable_states(p, p_threshold = collapse_states_threshold)
    state_labels <- used_states
  }
  # check if number of states exceeds default color palette options
  if (pie_colors[1] == "default" & length(state_labels) > 12) {
    stop(paste0(length(state_labels), " states in dataset; please provide colors (default only can provide up to 12)"))
  }

  # check if number of states not equal to provided colors
  if (pie_colors[1] != "default" & length(pie_colors) < length(state_labels)) {
    stop(paste0("You provided fewer colors in node_color than states in your dataset. There are ",
                length(state_labels), " states and you provide ", length(pie_colors), " colors."))
  }
  if (pie_colors[1] != "default" & length(pie_colors) > length(state_labels)) {
    stop(paste0("You provided more colors in node_color than states in your dataset. There are ",
                length(state_labels), " states and you provide ", length(pie_colors), " colors."))
  }

  # set default colors
  if (any(pie_colors == "default")) {
      nstates <- length(state_labels)
      colors <- .colFun(nstates)
      names(colors) <- state_labels
  } else {
    colors <- pie_colors
  }

  ##### reformat labels if necessary #####
  if (tip_labels_remove_underscore) {p$data$label <- gsub("_", " ", p$data$label)}

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

  # collapse states
  if (collapse_states == TRUE){
    p <- p + ggplot2::scale_color_manual(values = colors,
                                         breaks = state_labels,
                                         name = "Range",
                                         limits = used_states)
  } else {
    p <- p + ggplot2::scale_color_manual(values = colors, breaks = state_labels)
  }

  # set up guides
  # gotta do this for a legend!
  if (cladogenetic == TRUE) {

    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_3), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_other), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)

  } else if (cladogenetic == FALSE & "anc_state_1" %in% colnames(t@data)) {

    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_3), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_other), size=0),na.rm=TRUE, alpha=0.0)

  } else {

    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_3), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_other), size=0),na.rm=TRUE, alpha=0.0)
  }

  p <- p + ggplot2::guides(colour = ggplot2::guide_legend("State", override.aes = list(size=4, alpha = 1.0)), order=1)
  p <- p + ggplot2::guides(size=FALSE)

  # plot pies at nodes (and shoulders)
  if (cladogenetic == TRUE) {
    # create state matrices (matrix of nodes (rows) and all possible states (columns), values are pp. )
    state_probs <- .build_state_probs(t, state_labels, include_start_states = TRUE)
    dat_state_end <- state_probs$end
    dat_state_start <- state_probs$start

    # make pies
    pies_start <- ggtree:::nodepie(dat_state_start, cols = 1:(ncol(dat_state_start) - 1),
                                  color = colors, alpha = state_transparency)
    pies_end <- ggtree:::nodepie(dat_state_end, cols = 1:(ncol(dat_state_end) - 1),
                                color = colors, alpha = state_transparency)

    # add pies to tree

    # add node pies
    p <- .inset.revgadgets(
      tree_view = p,
      insets = pies_end[node_idx],
      x = "node",
      height = node_pie_size,
      width = node_pie_size,
      hjust = node_pie_nudge_x,
      vjust = node_pie_nudge_y
    )

    # add shoulder pies
    p <- .inset.revgadgets(
      tree_view = p,
      insets = pies_start,
      x = "parent_shoulder",
      height = shoulder_pie_size,
      width = shoulder_pie_size,
      hjust = shoulder_pie_nudge_x,
      vjust = shoulder_pie_nudge_y
    )

    # add tip pies
    if (tip_pies == TRUE) {
      p <- .inset.revgadgets(
        tree_view = p,
        insets = pies_end[tip_idx],
        x = "node",
        height = tip_pie_size,
        width = tip_pie_size,
        hjust = tip_pie_nudge_x,
        vjust = tip_pie_nudge_y
      )
    }

  } else {

    # create state matrices (matrix of nodes (rows) and all possible states (columns), values are pp. )
    dat_state_anc <- .build_state_probs(t, state_labels, include_start_states = FALSE)[[1]]

    pies_anc <- ggtree::nodepie(dat_state_anc, cols=1:(ncol(dat_state_anc)-1),
                               color = colors, alpha = state_transparency)

    # add pies to tree

    p <- .inset.revgadgets(
      tree_view = p,
      insets = pies_anc[node_idx],
      x = "node",
      height = node_pie_size,
      width = node_pie_size,
      hjust = node_pie_nudge_x,
      vjust = node_pie_nudge_y
    )
    if (tip_pies == TRUE) {
      p <- .inset.revgadgets(
        tree_view = p,
        insets = pies_anc[tip_idx],
        x = "node",
        height = tip_pie_size,
        width = tip_pie_size,
        hjust = tip_pie_nudge_x,
        vjust = tip_pie_nudge_y
      )
    }
  }


  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    p <- p + ggtree::geom_nodelab(ggplot2::aes(label = .convertAndRound(posterior) ), hjust = "left",
                                    nudge_x = node_labels_offset, size = node_labels_size)
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

  # add space on x axis for tip labels
  if (tip_labels == TRUE) {
    tree_height <- max(phytools::nodeHeights(t@phylo))
    p <- p + ggtree::xlim(0, tree_height + tree_height/2)
  }

  return(p)
}


