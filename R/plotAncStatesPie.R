#' plot Anc States Pie
#'
#' Plot character states and posterior probabilities as pies on nodes.
#'
#' @param t (treedata object; none) Output of processAncStates() function
#' containing tree and ancestral states.
#' @param cladogenetic (logical; FALSE) Plot shoulder pies of cladogenetic analyses?
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
#' values are "state" for the MAP ancestral states, "node_posterior" for the posterior
#' probability of the node on the tree, or NULL for not plotting any text at the nodes (default).
#' @param node_labels_size (numeric; 2) Size of node labels text. Ignored if
#' node_labels_as = NULL.
#' @param node_labels_offset (numeric; 0.1) Horizontal offset of node labels from nodes.
#' Ignored if node_labels_as = NULL.
#' @param pie_colors ("character"; "default") Colors for states in pies. If "default", plots
#' the default RevGadgets colors. Provide a character vector of hex codes or other R-readible
#' colors the same length of the number of character states. Names of the vector should
#' correspond to state labels.
#' @param node_pie_size (numeric; 2) Size of the pies at nodes.
#' @param tip_pie_size (numeric; 1) Size of the pies at tips.
#' @param shoulder_pie_size (numeric; node_pie_size) Size of the pies at shoulders for
#' cladogenetic plots.
#' @param tip_pies (logical; TRUE) Plot pies tips?
#' @param node_pie_nudge_x (numeric; 0) If pies aren't centered, ajust by nudging
#' @param node_pie_nudge_y (numeric; 0) If pies aren't centered, ajust by nudging
#' @param tip_pie_nudge_x (numeric; node_pie_nudge_x) If pies aren't centered, ajust by nudging
#' @param tip_pie_nudge_y (numeric; node_pie_nudge_y) If pies aren't centered, ajust by nudging
#' @param shoulder_pie_nudge_x (numeric; node_pie_nudge_x) If pies aren't centered, ajust by nudging
#' @param shoulder_pie_nudge_y (numeric; node_pie_nudge_y) If pies aren't centered, ajust by nudging
#' @param collapse_states (logical; FALSE) Collapse low-probability states into "other" category?
#' @param collapse_states_threshold (numeric; 0.05) Probability threshold for collapsing states.
#' Varies from 0 to 1.
#' @param state_transparency (integer; 0.75) Alpha (transparency) of state symbols- varies from
#' 0 to 1.
#'
#' @examples
#'
#' \dontrun{
#'
#' # Standard ancestral state reconstruction example
#'
#' # process file and assign state labels
#' file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
#` example <- processAncStates(file,
#'                                     state_labels = c("1" = "Awesome",
#'                                                      "2" = "Beautiful",
#'                                                      "3" = "Cool!"))
#' # plot
#' plotAncStatesPie(t = example)
#'
#' # DEC Biogeographic range evolution example
#'
#' # process file
#' file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
#'
#' # labels that correspond to each region/ possible combination of regions
#' labs <- c("1" = "K", "2" = "O", "3" = "M", "4" = "H", "5" = "KO",
#'           "6" = "KM", "7" = "OM", "8" = "KH", "9" = "OH", "10" = "MH",
#'           "11" = "KOM", "12" = "KOH", "13" = "KMH", "14" = "OMH", "15" = "KOMH")
#' # Use the state_labels in the returned tidytree object to define color palette
#' # These state_labels may be a subset of the labels you provided
#' # (not all possible regions may be sampled in the dataset)
#' colors <- colorRampPalette(RevGadgets:::.colFun(12))(length(dec_example@state_labels))
#' names(colors) <- dec_example@state_labels
#'
#' # plot
#' p <- plotAncStatesPie(t = dec_example, pie_colors = colors, tip_labels_size = 3,
#'              cladogenetic = TRUE, tip_labels_offset = 0.25) +
#'              ggplot2::theme(legend.position = c(0.1, 0.75))
#' }
#'
#' @export

plotAncStatesPie <- function(t,
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
                    node_pie_size = 1,
                    shoulder_pie_size = node_pie_size,
                    tip_pies = TRUE,
                    tip_pie_size = 0.5,

                    # nudges to center pies
                    node_pie_nudge_x = 0,
                    node_pie_nudge_y = 0,
                    tip_pie_nudge_x = node_pie_nudge_x,
                    tip_pie_nudge_y = node_pie_nudge_y,
                    shoulder_pie_nudge_x = node_pie_nudge_x,
                    shoulder_pie_nudge_y = node_pie_nudge_y,

                    state_transparency = 0.75) {

  ##### parameter compatibility checks #####
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
  if (is.character(pie_colors) == FALSE) stop ("pie_colors should be 'default' or valid color(s)")
  if (pie_colors[1] != "default" & any(.isColor(pie_colors) == FALSE) ) stop("pie_colors should be valid color(s)")
  if (any(is.numeric(node_pie_size) == FALSE)) stop("node_pie_size should be a single number")
  if (length(node_pie_size) > 1) stop("node_pie_size should be a single number")
  if (any(is.numeric(shoulder_pie_size) == FALSE)) stop("shoulder_pie_size should be a single number")
  if (length(shoulder_pie_size) > 1) stop("shoulder_pie_size should be a single number")
  if (any(is.numeric(tip_pie_size) == FALSE)) stop("tip_pie_size should be a single number")
  if (length(tip_pie_size) > 1) stop("tip_pie_size should be a single number")
  if (is.numeric(node_pie_nudge_x) == FALSE) stop("node_pie_nudge_x should be a single number")
  if (is.numeric(node_pie_nudge_y) == FALSE) stop("node_pie_nudge_y should be a single number")
  if (is.numeric(tip_pie_nudge_x) == FALSE) stop("tip_pie_nudge_x should be a single number")
  if (is.numeric(tip_pie_nudge_y) == FALSE) stop("tip_pie_nudge_y should be a single number")
  if (is.numeric(shoulder_pie_nudge_x) == FALSE) stop("shoulder_pie_nudge_x should be a single number")
  if (is.numeric(shoulder_pie_nudge_y) == FALSE) stop("shoulder_pie_nudge_y should be a single number")
  if (length(node_pie_nudge_x) != 1) stop("node_pie_nudge_x should be a single number")
  if (length(node_pie_nudge_y) != 1) stop("node_pie_nudge_y should be a single number")
  if (length(tip_pie_nudge_x) != 1) stop("tip_pie_nudge_x should be a single number")
  if (length(tip_pie_nudge_y) != 1) stop("tip_pie_nudge_y should be a single number")
  if (length(shoulder_pie_nudge_x) != 1) stop("shoulder_pie_nudge_x should be a single number")
  if (length(shoulder_pie_nudge_y) != 1) stop("shoulder_pie_nudge_y should be a single number")
  if (is.numeric(state_transparency) == FALSE) stop("state_transparency should be a number between 0 - 1")
  if (state_transparency < 0 || state_transparency > 1) stop("state_transparency should be a number between 0 - 1")

    ##### create basic tree plot #####
  p <- ggtree:::ggtree(t, ladderize = TRUE)
recover()
  ##### calculate helper variables #####
  tree <- attributes(t)$phylo
  tree_height <- max(phytools::nodeHeights(t@phylo))
  n_node <- ggtree:::getNodeNum(tree)
  n_tip <- length(tree$tip.label)
  node_idx <- (n_tip+1):n_node
  tip_idx <- 1:n_tip
  all_idx <- 1:n_node

  ##### calculate pie sizes #####
  # multiply by 1.05? uniform padding of 5% on each side of plot
  node_pie_size <-  ((n_tip * tree_height) / 15 ) * node_pie_size
  shoulder_pie_size <- ((n_tip * tree_height) / 15 ) * shoulder_pie_size
  tip_pie_size <- ((n_tip * tree_height) / 15 ) * tip_pie_size

  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE & "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE & "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }

  ##### color and label processing #####

  # check if number of states exceeds default color palette options
  if (pie_colors[1] == "default" &
      length(t@state_labels) > 12) {
    stop(paste0(length(t@state_labels),
                " states in dataset; please provide colors
                (default only can provide up to 12)"))
  }

  # check if number of states not equal to provided colors
  if (pie_colors[1] != "default" &
      length(pie_colors) < length(t@state_labels)) {
    stop(paste0("You provided fewer colors in node_color
                than states in your dataset. There are ",
                length(t@state_labels), " states and you provide ",
                length(pie_colors), " colors."))
  }

  # set colors, add "other" if necessary
  otherpp <- as.numeric(dplyr::pull(p$data,
                                    var = paste0(state_pos_str_base[1],
                                                 "other_pp")))
  if (sum(otherpp, na.rm = TRUE) == 0) {
    state_labels <- t@state_labels
    # set default colors
    if (any(pie_colors == "default")) {
      nstates <- length(state_labels)
      colors <- c(.colFun(nstates))
      names(colors) <- state_labels
    } else {
      colors <- pie_colors
    }

  } else if (sum(otherpp, na.rm = TRUE) != 0) {

    state_labels <- c(t@state_labels, "other")

    if ("anc_state_" %in% state_pos_str_base) {
      p$data$anc_state_other <- "other"
    }
    if ("end_state_" %in% state_pos_str_base) {
      p$data$end_state_other <- "other"
    }

    # add other to user-set colors
    if (pie_colors[1] != "default") {
      pc_names <- names(pie_colors)
      pie_colors <- c(pie_colors, "grey50")
      names(pie_colors) <- c(pc_names, "other")
    }

    # set default colors
    if (any(pie_colors == "default")) {
      nstates <- length(state_labels) - 1
      colors <- c(.colFun(nstates), "grey50")
      names(colors) <- state_labels
    } else {
      colors <- pie_colors
    }
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

  # set up guides
  # gotta do this for a legend!
  if (cladogenetic == TRUE) {
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_3), size=0),na.rm=TRUE, alpha=0.0)
    if ("other" %in% state_labels) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_other), size=0),na.rm=TRUE, alpha=0.0)
    }
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)

  } else if (cladogenetic == FALSE & "anc_state_1" %in% colnames(t@data)) {

    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_3), size=0),na.rm=TRUE, alpha=0.0)
    if ("other" %in% state_labels) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_other), size=0),na.rm=TRUE, alpha=0.0)
    }

  } else {
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_3), size=0),na.rm=TRUE, alpha=0.0)
    if ("other" %in% state_labels) {
      p <- p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_other), size=0),na.rm=TRUE, alpha=0.0)
    }
  }

  p <- p + ggplot2::scale_color_manual(values = colors, breaks = state_labels)
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
    if (sum(otherpp, na.rm = TRUE) == 0) {dat_state_anc$other <- NULL}
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
    p <- p + ggtree::xlim(0, tree_height + tree_height/2)
  }

  return(p)
}


