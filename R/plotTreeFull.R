#' Plot Full tree
#'
#' Plots a tree, such as an MCC or MAP tree
#'
#' Plots a single tree, such as an MCC or MAP tree, with
#' the full set of functionality to be called by plotTree()
#' and plotFBDTree()
#'
#' @param tree (list of lists of treedata objects; no default) Name of a list
#' of lists of treedata objects, such as produced by readTrees(). This object
#' should only contain only one summary tree from one trace file. If it
#' contains multiple trees or multiple traces, only the first will be used.
#'
#' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with
#' timescale in MYA.
#'#'
#' @param geo (logical; timeline) Add a geological timeline? Defaults to the
#' same as timeline.
#'
#' @param time_bars (logical; timeline) Add vertical gray bars to indicate
#' geological timeline units if geo == TRUE or regular time intervals (in MYA)
#' if geo == FALSE.
#'
#' @param geo_units (list; list("epochs", "periods")) Which geological units to
#' include in the geo timescale. May be "periods", "epochs", "stages", "eons", 
#' "eras", or a list of two of those units.
#'
#' @param node_age_bars (logical; TRUE) Plot time tree with node age bars?
#'
#' @param age_bars_colored_by (character; NULL) Specify column to color
#' node/tip age bars by, such as "posterior". If null, all age bars plotted the
#' same color, specified by age_bars_color
#'
#' @param age_bars_color (character; "blue") Color for node/tip age bars.
#' If age_bars_colored_by specifies a variable (not NULL), you must provide
#' two colors, low and high values for a gradient. Colors must be either R
#' valid color names or valid hex codes.
#' 
#' @param age_bars_width (numeric; 1) Change line width for age bars
#'
#' @param node_labels (character; NULL) Plot text labels at nodes, specified by
#' the name of the corresponding column in the tidytree object. If NULL, no
#' text is plotted.
#'
#' @param node_labels_color (character; "black") Color to plot node_labels,
#' either as a valid R color name or a valid hex code.
#'
#' @param node_labels_size (numeric; 3) Size of node labels
#'
#' @param tip_labels (logical; TRUE) Plot tip labels?
#'
#' @param tip_labels_italics (logical; FALSE) Plot tip labels in italics?
#'
#' @param tip_labels_formatted (logical; FALSE) Do the tip labels contain 
#' manually added formatting information? Will set parse = TRUE in geom_text()
#' and associated functions to interpret formatting. See ?plotmath for more.
#' Cannot be TRUE if tip_labels_italics = TRUE.  
#'
#' @param tip_labels_remove_underscore (logical; FALSE) Should underscores be
#' replaced by spaces in tip labels?
#'
#' @param tip_labels_color (character; "black") Color to plot tip labels, either
#' as a valid
#' R color name or a valid hex code.
#'
#' @param tip_labels_size (numeric; 3) Size of tip labels
#'
#' @param node_pp (logical; FALSE) Plot posterior probabilities as symbols at
#' nodes? Specify symbol aesthetics with node_pp_shape, node_pp_color, and
#' node_pp_size.
#'
#' @param node_pp_shape (integer; 1) Integer corresponding to point shape
#' (value between 0-25). See ggplot2 documentation for details:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}
#'
#' @param node_pp_color (character; "black") Color for node_pp symbols,
#' either as valid R color name(s) or hex code(s). Can be a single character
#' string specifying a single color, or a vector of length two specifying two
#' colors to form a gradient. In this case, posterior probabilities will be
#' indicated by color along the specified gradient.
#'
#' @param node_pp_size (numeric or character; 1) Size for node_pp symbols.
#' If numeric, the size will be fixed at the specified value. If a character,
#' it should specify "variable", indicating that size should be scaled by the
#' posterior value. Size regulates the area of the shape, following ggplot2
#' best practices:
#' \url{https://ggplot2.tidyverse.org/reference/scale_size.html})
#'
#' @param tip_age_bars (logical; FALSE) Plot node age bars for the tips as
#' well? Useful for plotting serial sampled analyses or fossilized birth-death
#' analyses, or any cases where some tip ages are estimated.
#'
#' @param branch_color (character; "black") A single character string
#' specifying the color (R color name or hex code) for all branches OR a vector
#' of length 2 specifying two colors for a gradient, used to color the branches
#' according to the variable specified in color_branch_by. If only 1 color is
#' provided and you specify color_branch_by, default colors will be chosen
#' (low = "#005ac8", high = "#fa7850").
#'
#' @param color_branch_by (character; NULL ) Optional name of one quantitative
#' variable in the treedata object to color branches, such as a rate.
#'
#' @param line_width (numeric; 1) Change line width for branches
#'
#' @param label_sampled_ancs (logical; FALSE) Label any sampled ancestors? Will
#' inherent tip labels aesthetics for size and color.
#'
#' # added in from plotTree
#' @param tree_layout (character; "rectangular") Tree shape layout, passed to
#' ggtree(). Options are 'rectangular', 'cladogram', 'slanted', 'ellipse',
#' 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle',
#' 'daylight' or 'ape'.
#'
#' @param node_labels_offset (numeric; 0) Horizontal offset of node labels
#' from nodes.
#'
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from
#' tree.
#'
#' @param ... (various) Additional arguments passed to ggtree::ggtree().
#'
#' @seealso called by \link[RevGadgets]{plotTree} and
#' \link[RevGadgets]{plotFBDTree}
#'
#' @return returns a single plot object.



plotTreeFull <- function(tree,

                         timeline,
                         geo,
                         geo_units,
                         time_bars,

                         node_age_bars,
                         tip_age_bars,
                         age_bars_color,
                         age_bars_colored_by,
                         age_bars_width,

                         node_labels,
                         node_labels_color,
                         node_labels_size,
                         node_labels_offset,

                         tip_labels,
                         tip_labels_italics,
                         tip_labels_formatted,
                         tip_labels_remove_underscore,
                         tip_labels_color,
                         tip_labels_size,
                         tip_labels_offset,

                         label_sampled_ancs,

                         node_pp,
                         node_pp_shape,
                         node_pp_color,
                         node_pp_size,

                         branch_color,
                         color_branch_by,
                         line_width,

                         tree_layout,
                         ...) {
  # enforce argument matching
  if (!is.list(tree))
    stop("tree should be a list of lists of treedata objects")
  if (!methods::is(tree[[1]][[1]], "treedata"))
    stop("tree should be a list of lists of treedata objects")
  vars <- colnames(tree[[1]][[1]]@data)
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (is.logical(node_age_bars) == FALSE)
    stop("node_age_bars should be TRUE or FALSE")
  if (any(.isColor(age_bars_color) == FALSE))
    stop("age_bars_color should be valid color(s)")
  if (is.null(age_bars_colored_by) == FALSE &
      any(vars %in% age_bars_colored_by) == FALSE)
    stop("age_bars_colored_by should be a column in your tidytree object")
  if (is.numeric(age_bars_width) == FALSE)
    stop ("age_bars_width should be numeric")
  if (is.null(node_labels) == FALSE &
      any(vars %in% node_labels) == FALSE)
    stop("node_labels should be NULL or a column in your tidytree object")
  if (is.null(node_labels_color) == FALSE &
      .isColor(node_labels_color) == FALSE)
    stop("node_labels_color should be NULL or a recognized color")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (.isColor(tip_labels_color) == FALSE)
    stop("tip_labels_color should be a recognized color")
  if (!methods::is(node_pp,"logical"))
    stop("node_pp should be TRUE or FALSE")
  if (node_pp) {
    if (length(node_pp_color) > 2)
      stop("node_pp_color should be of length 1 or 2")
    if (.isColor(node_pp_color) == FALSE)
      stop("node_pp_color should be a recognized color")
    if (node_pp_shape %in% 0:25 == FALSE)
      stop("node_pp_shape should be a recognized shape
           (value between 0 and 25)")
    if (is.numeric(node_pp_size) == FALSE &
        node_pp_size != "variable")
      stop("node_pp_size should be numeric or 'variable'")
  }
  if (is.logical(tip_age_bars) == FALSE)
    stop("tip_age_bars should be TRUE or FALSE")
  if (length(branch_color) == 1 &
      !.isColor(branch_color))
    stop("branch_color should be a recognized color")
  if (length(branch_color) == 2) {
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2]) == FALSE)
      stop("Neither values of branch_color are a recognized color")
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2]))
      stop("branch_color[1] is not a recognized color")
    if (.isColor(branch_color[1]) &
        .isColor(branch_color[2]) == FALSE)
      stop("branch_color[2] is not a recognized color")
  } else if (length(branch_color) > 2)
    stop("only 2 colors may be specified in branch_color")
  if (is.null(color_branch_by) == FALSE &
      any(vars %in% color_branch_by) == FALSE)
    stop("color_branch_by should be NULL or a column in your tidytree object")
  if (is.numeric(line_width) == FALSE)
    stop ("line_width should be numeric")
  if (is.logical(label_sampled_ancs) == FALSE)
    stop("label_sampled_ancs should be TRUE or FALSE")
  if (is.list(geo_units)) {
    if (length(geo_units) != 2)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[1]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[2]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  } else {
    if (geo_units %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras') == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  }
  if (is.numeric(tip_labels_offset) == FALSE)
    stop ("tip_labels_offset should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop ("node_labels_offset should be a number")
  tree_layout <-
    match.arg(
      tree_layout,
      choices = c(
        'rectangular',
        'slanted',
        'ellipse',
        'cladogram',
        'roundrect',
        'fan',
        'circular',
        'inward_circular',
        'radial',
        'equal_angle',
        'daylight',
        'ape'
      )
    )
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) {
      stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")
    }
  }
  # grab single tree from input
  phy <- tree[[1]][[1]]

  ### fix for trees with sampled ancestors ###
  phylo    <- phy@phylo
  node_map <- .matchNodesTreeData(phy, phylo)
  if ("sampled_ancestor" %in% colnames(phy@data)) {
    phy@data$node <-
      as.character(node_map[match(as.numeric(phy@data$index),
                                  node_map$Rev), ]$R)
  }

  ### set up tree layout ###

  if (tree_layout == "cladogram") {
    tree_layout <- "rectangular"
    BL <- "none"
  } else {
    BL <- "branch.length"
  }

  # initiate plot
  if (is.null(color_branch_by)) {
    pp <- ggtree::ggtree(
      phy,
      right = FALSE,
      size = line_width,
      color = branch_color,
      branch.length = BL,
      layout = tree_layout,
      ...
    )
  } else if (!is.null(color_branch_by)) {
    pp <- ggtree::ggtree(
      phy,
      right = FALSE,
      size = line_width,
      branch.length = BL,
      layout = tree_layout,
      ...
    )
  }

  #### parameter compatibility checks ###
  if (length(node_pp_color) == 2 &
      length(branch_color) == 2)
    stop(
      "You may only include variable colors for either node_pp_label or
      branch_color, not for both"
    )

  #check that if user wants node_age_bars, there are dated intervals in  file
  if (node_age_bars == TRUE) {
    if (!"age_0.95_HPD" %in% colnames(phy@data))
      stop(
        "You specified node_age_bars, but there is no age_0.95_HPD column
        in the treedata object."
      )
  }

  # get dimensions
  n_node <- treeio::Nnode(phy)
  tree_height <- max(phytools::nodeHeights(phy@phylo))
  ntips <- sum(pp$data$isTip)

  # reformat labels if necessary
  if (tip_labels_remove_underscore) {
    pp$data$label <- gsub("_", " ", pp$data$label)
  }

  # check that if user wants to label sampled ancs,
  # there are sampled ancs in the files
  if (label_sampled_ancs == TRUE &
      "sampled_ancestor" %in% colnames(pp$data)) {
    sampled_ancs <-
      pp$data[!pp$data$isTip & !is.na(pp$data$sampled_ancestor),]
    if (nrow(sampled_ancs) < 1) {
      label_sampled_acs <- FALSE
    }
   }

  # add timeline
  if (timeline == TRUE) {
    warning("Plotting with default axis label (Age (Ma))")
    if (node_age_bars == FALSE) {
      minmax <- phytools::nodeHeights(phy@phylo)
      max_age <- tree_height
    } else {
      pp$data$age_0.95_HPD <- lapply(pp$data$age_0.95_HPD, function(z) {
        if (any(is.null(z)) ||
            any(is.na(z))) {
          return(c(NA, NA))
        } else {
          return(as.numeric(z))
        }
      })
      minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
      max_age <- max(minmax, na.rm = TRUE)
    }

    interval <- max_age / 5
    dx <- max_age %% interval

    # add geo
    tick_height <- ntips / 100
    if (geo == TRUE) {
      #determine whether to include quaternary
      if (tree_height > 50) {
        skipit <- c("Quaternary", "Holocene", "Late Pleistocene")
      } else {
        skipit <- c("Holocene", "Late Pleistocene")
      }
      # add deep timescale
      if (length(geo_units) == 1) {
        pp <- pp + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-max(minmax, na.rm = TRUE), tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          height = grid::unit(4, "line"),
          skip = skipit,
          abbrv = FALSE,
          rot = 90,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      } else if (length(geo_units) == 2) {
        pp <- pp + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-max(minmax, na.rm = TRUE), tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          skip = skipit,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      }
    }
    #add axis title
    pp <- pp + ggplot2::scale_x_continuous(
      name = "Age (Ma)",
      expand = c(0, 0),
      limits = c(-max_age, tree_height /
                   2),
      breaks = -rev(seq(0, max_age +
                          dx, interval)),
      labels = rev(seq(0, max_age +
                         dx, interval))
    )
    pp <- ggtree::revts(pp)

    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )

    pp <-
      pp + ggplot2::geom_segment(ggplot2::aes(
        x = 0,
        y = -tick_height * 5,
        xend = -max_age,
        yend = -tick_height * 5
      )) +
      ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      )) +
      ggplot2::annotate(
        "text",
        x = -rev(xline),
        y = -tick_height * 5 + tick_height * 2,
        label = rev(xline),
        size = tip_labels_size
      )

    # add vertical gray bars
    if (time_bars) {
      if (geo) {
        if ("epochs" %in% geo_units) {
          x_pos <- -rev(c(0, deeptime::get_scale_data("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::get_scale_data("periods")$max_age))
        }
      } else if (!geo) {
        x_pos <- -rev(xline)
      }
      for (k in 2:(length(x_pos))) {
        box_col <- "gray92"
        if (k %% 2 == 1)
          box_col <- "white"
        box <-
          ggplot2::geom_rect(
            xmin = x_pos[k - 1],
            xmax = x_pos[k],
            ymin = -tick_height * 5,
            ymax = ntips,
            fill = box_col
          )
        pp <- gginnards::append_layers(pp, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      pp <-
        pp +
        ggplot2::theme(axis.title.x =
                         ggplot2::element_text(hjust = max_age / (2 * tot)))
    }
  }

  # processing for node_age_bars and tip_age_bars
  if (node_age_bars == TRUE) {
    # Encountered problems with using geom_range to plot age HPDs in ggtree. It
    # appears that geom_range incorrectly rotates the HPD relative to the height
    # of the node unnecessarily. My guess for this would be because older
    # version of ggtree primarily supported length measurements, and not
    # height measurements so the new capability to handle height might contain
    # a "reflection" bug.
    # For example, suppose a node has height 3 with HPD [2, 7]. You can think of
    # this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
    # causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
    #
    # See this excellent trick by Tauana:
    # https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
    # Adapted this code to also plot fossil tip uncertainty in red
    pp$data$age_0.95_HPD <-
      lapply(pp$data$age_0.95_HPD, function(z) {
        if (any(is.null(z)) ||
            any(is.na(z))) {
          return(c(NA, NA))
        } else {
          return(as.numeric(z))
        }
      })

    minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
    bar_df <-
      data.frame(
        node_id = as.integer(pp$data$node),
        isTip = pp$data$isTip,
        as.data.frame(minmax)
      )
    names(bar_df) <- c("node_id", "isTip", "min", "max")
    if (tip_age_bars == TRUE) {
      tip_df <- dplyr::filter(bar_df, isTip == TRUE & !is.na(min))
      tip_df <-
        dplyr::left_join(tip_df, pp$data, by = c("node_id" = "node"))
      tip_df <- dplyr::select(tip_df, node_id, min, max, y)
    }
    if (is.null(age_bars_colored_by) == TRUE) {
      # plot age densities
      bar_df <-
        dplyr::left_join(bar_df, pp$data, by = c("node_id" = "node"))
      bar_df <- dplyr::select(bar_df,  node_id, min, max, y)
      pp <-
        pp + ggplot2::geom_segment(
          ggplot2::aes(
            x = -min,
            y = y,
            xend = -max,
            yend = y
          ),
          data = bar_df,
          color = age_bars_color,
          linewidth = age_bars_width,
          alpha = 0.8
        )
    } else if (is.null(age_bars_colored_by) == FALSE) {
      if (length(age_bars_color) == 1) {
        age_bars_color <- colFun(2)[2:1]
      }

      if ("sampled_ancestor" %in% colnames(pp$data) == TRUE) {
        sampled_tip_probs <-
          1 - as.numeric(pp$data$sampled_ancestor[pp$data$isTip == TRUE])
        sampled_tip_probs[is.na(sampled_tip_probs)] <- 0
      } else {
        sampled_tip_probs <- rep(1, sum(pp$data$isTip))
      }

      pp$data$olena <-
        c(sampled_tip_probs,
          as.numeric(.convertAndRound(L =
                                        unlist(pp$data[pp$data$isTip == FALSE,
                                                       age_bars_colored_by]))))

      bar_df <-
        dplyr::left_join(bar_df, pp$data, by = c("node_id" = "node"))
      bar_df <-
        dplyr::select(bar_df,  node_id, min, max, y, olena, isTip = isTip.x)
      if (tip_age_bars == FALSE) {
        bar_df <- dplyr::filter(bar_df, isTip == FALSE)
      }
      pp <-
        pp + ggplot2::geom_segment(
          ggplot2::aes(
            x = -min,
            y = y,
            xend = -max,
            yend = y,
            color = olena
          ),
          data = bar_df,
          linewidth = age_bars_width,
          alpha = 0.8
        ) +
        ggplot2::scale_color_gradient(
          low = age_bars_color[1],
          high = age_bars_color[2],
          name = .titleFormat(age_bars_colored_by),
          breaks = pretty(pp$data$olena)
        )
    }
  }

  # label sampled ancestors
  if (label_sampled_ancs == TRUE &
      "sampled_ancestor" %in% colnames(pp$data)) {
    sampled_ancs <-
      pp$data[!pp$data$isTip & !is.na(pp$data$sampled_ancestor),]
    space_labels <- ntips / 30
    if (tip_labels_italics == TRUE) {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": italic(`",
                         sampled_ancs$label,
                         "`)"),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0,
          parse = TRUE
        )
    } else if (tip_labels_italics == TRUE) {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": ",
                         sampled_ancs$label),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0,
          parse = TRUE
        )
    } else {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": ",
                         sampled_ancs$label),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0
        )
    }

    t_height <- ntips / 200
    df <- data.frame(
      x = sampled_ancs$x,
      vx = sampled_ancs$x,
      y = sampled_ancs$y + t_height,
      vy = sampled_ancs$y - t_height
    )
    pp <-
      pp + ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      ))

  }

  # add node labels (text)
  if (is.null(node_labels) == FALSE) {
    # catch some funkiness from importing an unrooted tree
    if (node_labels == "posterior") {
      pp$data[grep("]:", unlist(pp$data[, node_labels])), node_labels] <-
        NA
    }
    pp$data$kula <-
      c(rep(NA, times = ntips),
        .convertAndRound(L =
                           unlist(pp$data[pp$data$isTip == FALSE,
                                          node_labels])))

    # change any NAs that got converted to characters back to NA
    pp$data$kula[pp$data$kula == "NA"] <- NA
    pp <- pp + ggtree::geom_text(
      ggplot2::aes(label = kula),
      color = node_labels_color,
      nudge_x = node_labels_offset,
      hjust = 0,
      size = node_labels_size
    )
  }

  # add tip labels (text)
  if (tip_labels == TRUE) {
    if (tip_age_bars == TRUE) {
      pp$data$extant <- !pp$data$node %in% tip_df$node_id
    } else {
      pp$data$extant <- TRUE
    }
    if (tip_labels_italics == TRUE) {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(
          subset = extant & isTip,
          label = paste0('italic(`', label, '`)')
        ),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color,
        parse = TRUE
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = paste0('italic(`', new_tip_df$label, '`)'),
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size,
            parse = TRUE
          )
      }
    } else if (tip_labels_formatted == TRUE ) {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(subset = extant & isTip,
                     label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color,
        parse = TRUE
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = new_tip_df$label,
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size,
            parse = TRUE 
          )
      }
    } else {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(subset = extant & isTip,
                     label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = new_tip_df$label,
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size
          )
      }
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
        ggplot2::scale_size_continuous(name = "Posterior")
    } else if (length(node_pp_color) == 2 &
               node_pp_size != "variable") {
      pp <- pp + ggtree::geom_nodepoint(size = node_pp_size,
                                        ggplot2::aes(color = posterior),
                                        shape = node_pp_shape) +
        ggplot2::scale_color_gradient(
          name = "Posterior",
          low = node_pp_color[1],
          high = node_pp_color[2],
          breaks = pretty(pp$data$posterior)
        )
    }
  }

  # add branch coloration by variable
  if (is.null(color_branch_by) == FALSE) {
    #set default colors if none provided
    if (length(branch_color) != 2) {
      branch_color <- c("#005ac8", "#fa7850")
    }
    col_num <- which(colnames(pp$data) == color_branch_by)
    pp$data[, col_num] <-
      as.numeric(as.data.frame(pp$data)[, col_num])
    name <- .titleFormat(color_branch_by)
    pp <- pp + 
      ggplot2::aes(color = as.data.frame(pp$data)[, col_num]) +
      ggplot2::scale_color_gradient(
        low = branch_color[1],
        high = branch_color[2],
        breaks = pretty(as.data.frame(pp$data)[, col_num]),
        name = name
      )
  }

  # readjust axis for non-timeline plots
  if (timeline == FALSE & BL != "none") {
    if (node_age_bars == FALSE) {
      xlim_min <- -tree_height
    } else {
      xlim_min <- -max(t(matrix(
        unlist(pp$data$age_0.95_HPD),
        nrow = 2
      )), na.rm = TRUE)
    }

    if (tip_labels == TRUE) {
      xlim_max <- tree_height / 2
    } else {
      xlim_max <- 0
    }

    pp <- pp + ggtree::xlim(xlim_min, xlim_max)
    pp <- ggtree::revts(pp)

  }

  # readjust axis for cladograms
  if (timeline == FALSE & BL == "none") {
    xlim_min <- range(pp$data$x)[1]

    if (tip_labels == TRUE) {
      xlim_max <- range(pp$data$x)[2] * 1.5
    } else {
      xlim_max <- range(pp$data$x)[2]
    }

    pp <- pp + ggtree::xlim(xlim_min, xlim_max)

  }

  return(pp)
}
