#' plot Ancestral States MAP
#'
#' Plots the MAP estimates of ancestral states. Can accommodate cladogenetic
#' reconstructions by plotting on shoulders. Defaults to varying the symbols by
#' color to indicate estimated ancestral state and varying the size of the
#' symbol to indicate the posterior probability of that estimate, but symbol
#' shape may also vary to accommodate black and white figures. For more details
#' on the aesthetics options, see parameter details below. For data with many
#' character states (such as chromosome counts), vary the size of the symbol
#' by estimated ancestral state, and vary the posterior probability of that
#' estimate by a color gradient. Text labels at nodes and tips are also
#' available.
#'
#' @param t (treedata object; none) Output of processAncStates() function
#' containing tree and ancestral states.
#' @param cladogenetic (logical; FALSE) Plot shoulder states of cladogenetic
#' analyses?
#' @param tip_labels (logical; TRUE) Label taxa labels at tips?
#' @param tip_labels_size (numeric; 2) Size of tip labels.
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from
#' tree.
#' @param tip_labels_italics (logical; FALSE) Italicize tip labels?
#' @param tip_labels_remove_underscore (logical; TRUE) Remove underscores from
#' tip labels?
#' @param tip_labels_formatted (logical; FALSE) Do the tip labels contain 
#' manually added formatting information? Will set parse = TRUE in geom_text()
#' and associated functions to interpret formatting. See ?plotmath for more.
#' Cannot be TRUE if tip_labels_italics = TRUE.  
#' @param tip_labels_states (logical; FALSE) Optional plotting of text at tips
#' in addition
#' to taxa labels.
#' @param tip_labels_states_size (numeric; 2) Size of state labels at tips.
#' Ignored if
#' tip_labels_states is FALSE.
#' @param tip_labels_states_offset (numeric; 0.1) Horizontal offset of tip
#' state labels.
#' Ignored if tip_labels_states = NULL.
#' @param node_labels_as (character; NULL) Optional plotting of text at nodes.
#' Possible values are "state" for the ancestral states , "state_posterior"
#' for posterior probabilities of the estimated ancestral state,
#' "node_posterior" or the posterior probability of the node on the tree,
#' or NULL for not plotting any text at the nodes (default).
#' @param node_labels_size (numeric; 2) Size of node labels text. Ignored if
#' node_labels_as = NULL.
#' @param node_labels_offset (numeric; 0.1) Horizontal offset of node labels
#' from nodes. Ignored if node_labels_as = NULL.
#' @param node_labels_centered (logical; FALSE) Should node labels be centered
#' over the nodes? Defaults to FALSE: adjusting node labels to the right of
#' nodes and left of shoulders.
#' @param node_size_as (character; "state_posterior") How to vary size of
#' node symbols. Options are "state_posterior" (default) for posterior
#' probabilities of the estimated ancestral state, "node_posterior" or the
#' posterior probability of the node on the tree, "state" for vary size by the
#' ancestral state itself in cases where there are many character states
#' (e.g. chromosome numbers; we do not recommend this option for characters
#' with few states), or NULL for fixed symbol size.
#' @param node_color_as (character; "state") How to vary to color of node
#' symbols. Options are "state" (default) to vary by estimated ancestral states,
#' "state_posterior" for posterior probabilities of the estimated ancestral
#' state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL to set all as one color.
#' @param node_shape_as (character; NULL) Option to vary node symbol by shape.
#' Options are NULL to keep shape constant or "state" to vary shape by
#' ancestral state.
#' @param node_shape (integer; 19) Shape type for nodes. If node_shape_as =
#' "state", provide a vector with length of the number of states. See ggplot2
#' documentation for details:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}
#' @param node_color ("character"; "default") Colors for node symbols. Defaults
#' to default RevGadgets colors. If node_color_as = "state', provide a vector of
#' length of the character states. If your color vector is labeled with state
#' labels, the legend will be displayed in the order of the labels. If
#' node_color_as = "posterior", provide a vector of length 2 to generate a
#' color gradient.
#' @param node_size (numeric; c(2, 6)) Range of sizes, or fixed size, for node
#' symbols. If node_size_as = "state_posterior", "node_posterior", or "state",
#' numeric vector of length two. If node_size_as = NULL, numeric vector of
#' length one. Size regulates the area of the symbol, following ggplot2 best
#' practices: \url{https://ggplot2.tidyverse.org/reference/scale_size.html})
#' @param tip_states (logical; TRUE) Plot states of taxa at tips?
#' @param tip_states_size (numeric; node_size) Size for tip symbols. Defaults
#' to the same size as node symbols.
#' @param tip_states_shape (integer; node_shape) Shape for tip symbols.
#' Defaults to the same as node symbols.
#' @param state_transparency (integer; 0.75) Alpha (transparency) of state
#' symbols- varies from 0 to 1.
#' @param tree_layout (character; "rectangular") Tree shape layout, passed to
#' ggtree(). Options are 'rectangular', 'slanted', 'ellipse', 'roundrect',
#' 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight'
#' or 'ape'. When cladogenetic = TRUE, only "rectangular" and 'circular' are
#' available.
#' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with
#' timescale in MYA.
#' @param geo (logical; timeline) Add a geological timeline? Defaults to the
#' same as timeline.
#' @param time_bars (logical; timeline) Add vertical gray bars to indicate
#' geological timeline units if geo == TRUE or regular time intervals (in MYA)
#' if geo == FALSE.
#' @param geo_units (list; list("epochs", "periods")) Which geological units to
#' include in the geo timescale. May be "periods", "epochs", "stages", "eons", 
#' "eras", or a list of two of those units.
#' @param ... (various) Additional arguments passed to ggtree::ggtree().
#'
#' @return A ggplot object
#'
#' @examples
#'
#' \donttest{
#' # Standard ancestral state reconstruction example with various aesthetics
#'
#' # process file
#' file <- system.file("extdata",
#'                     "comp_method_disc/ase_freeK.tree",
#'                     package="RevGadgets")
#' example <- processAncStates(file,
#'                             state_labels = c("1" = "Awesome",
#'                                              "2" = "Beautiful",
#'                                              "3" = "Cool!"))
#'
#' # have states vary by color and indicate state pp with size (default)
#' plotAncStatesMAP(t = example)
#'
#' # have states vary by color and indicate state pp with size ,
#' # and add a timeline
#' plotAncStatesMAP(t = example, timeline = TRUE)
#'
#' # have states vary by color and symbol, label nodes with pp of states
#' plotAncStatesMAP(t = example,  node_shape_as = "state",
#'                  node_size = 4, node_shape = c(15, 17,20),
#'                  node_size_as = NULL, node_labels_as = "state_posterior")
#'
#' # black and white figure - state as symbols and state pp with text
#' plotAncStatesMAP(t = example, node_color_as = NULL,
#'                  node_shape_as = "state", node_shape =  c(15, 17,20),
#'                  node_size_as = NULL, node_size = 4,
#'                  node_labels_as = "state_posterior",
#'                  node_color = "grey", state_transparency = 1)
#'
#' # default with circular tree
#' plotAncStatesMAP(t = example, tree_layout = "circular")
#'
#'
#' # Chromosome evolution example
#'
#' # process file
#' file <- system.file("extdata",
#'                     "chromo/ChromEvol_simple_final.tree",
#'                     package="RevGadgets")
#' chromo_example <- processAncStates(file, labels_as_numbers = TRUE)
#'
#' # plot
#' plotAncStatesMAP(t = chromo_example, node_color_as = "state_posterior",
#'                  node_size_as = "state", node_color = colFun(2),
#'                  tip_labels_offset = 0.005, node_labels_as = "state",
#'                  node_labels_offset = 0, tip_labels_states = TRUE,
#'                  tip_labels_states_offset = 0, tip_states = FALSE)
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
#' plotAncStatesMAP(t = dec_example,
#'                  cladogenetic = TRUE,
#'                  tip_labels_offset = 0.5)
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
                             tip_labels_formatted = FALSE,
                             tip_labels_remove_underscore = TRUE,

                             # label states at tips
                             tip_labels_states = FALSE,
                             tip_labels_states_size = 2,
                             tip_labels_states_offset = 0.1,

                             # text labels at nodes
                             node_labels_as = NULL,
                             node_labels_size = 2,
                             node_labels_offset = 0.1,
                             node_labels_centered = FALSE,

                             # what to plot at nodes
                             node_size_as = "state_posterior",
                             node_color_as = "state",
                             node_shape_as = NULL,

                             # aesthetics for plotting at nodes
                             node_shape = 19,
                             node_color = "default",
                             node_size = c(2, 6),

                             # aesthetics for tip states (inherents additional
                             # aesthetics from nodes)
                             tip_states = TRUE,
                             tip_states_size = node_size,
                             tip_states_shape = node_shape,

                             state_transparency = 0.75,
                             tree_layout = "rectangular",

                             timeline = FALSE,
                             geo = timeline,
                             geo_units = list("epochs", "periods"),
                             time_bars = timeline,

                             ...) {
  ##### parameter compatability checks! #####
  if (!methods::is(t, "treedata"))
    stop("t should be a treedata objects")
  if (is.logical(cladogenetic) == FALSE)
    stop("cladogenetic should be TRUE or FALSE")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.numeric(tip_labels_size) == FALSE)
    stop("tip_labels_size should be a number")
  if (is.numeric(tip_labels_offset) == FALSE)
    stop("tip_labels_offset should be a number")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (is.logical(tip_labels_remove_underscore) == FALSE)
    stop("tip_labels_remove_underscore should be TRUE or FALSE")
  if (is.logical(tip_labels_states) == FALSE)
    stop("tip_labels_states should be TRUE or FALSE")
  if (is.numeric(tip_labels_states_size) == FALSE)
    stop("tip_labels_states_size should be a number")
  if (is.numeric(tip_labels_states_offset) == FALSE)
    stop("tip_labels_states_offsetshould be a number")
  if (is.null(node_labels_as) == FALSE) {
    node_labels_as <-
      match.arg(node_labels_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.numeric(node_labels_size) == FALSE)
    stop("node_labels_size should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop("node_labels_offset should be a number")
  if (is.logical(node_labels_centered) == FALSE)
    stop("node_labels_centered should be TRUE or FALSE")
  if (is.null(node_size_as) == FALSE) {
    node_size_as <-
      match.arg(node_size_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_color_as) == FALSE) {
    node_color_as <-
      match.arg(node_color_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as != "state")
      stop("node_shape_as should be NULL or 'state'")
  }
  if (is.numeric(node_shape) == FALSE)
    stop("node_shape should be a number indicating symbol type")
  if (is.character(node_color) == FALSE)
    stop ("node_color should be 'default' or valid color(s)")
  if (node_color[1] != "default" &
      any(.isColor(node_color) == FALSE))
    stop("node_color should be valid color(s)")
  if (any(is.numeric(node_size) == FALSE))
    stop("node_size should be a single number or a vector of two numbers")
  if (length(node_size) > 2)
    stop("node_size should be a single number or a vector of two numbers")
  if (is.logical(tip_states) == FALSE)
    stop("tip_states should be TRUE or FALSE")
  if (is.numeric(tip_states_size) == FALSE)
    stop("tip_states_size should be a number")
  if (is.numeric(tip_states_shape) == FALSE)
    stop("tip_states_shape should be a number indicating symbol type")
  if (is.numeric(state_transparency) == FALSE)
    stop("state_transparency should be a number between 0 - 1")
  if (state_transparency > 1 |
      state_transparency < 0)
    stop("state_transparency should be a number between 0 - 1")
  if (cladogenetic == FALSE) {
    tree_layout <-
      match.arg(
        tree_layout,
        choices = c(
          'rectangular',
          'slanted',
          'ellipse',
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
  } else if (cladogenetic == TRUE) {
    tree_layout <-
      match.arg(tree_layout, choices = c('rectangular', 'circular'))
  }
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) {
      stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")
    }
  }
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
  ##### calculate helper variables #####
  tree <- attributes(t)$phylo

  ##### create basic tree plot #####
  p <- ggtree::ggtree(t, layout = tree_layout, ...)

  # get dimensions
  n_node <- ape::Nnode(tree, internal.only = FALSE)
  tree_height <- max(phytools::nodeHeights(t@phylo))
  ntips <- sum(p$data$isTip)

  ##### process column names #####
  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE &
             "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE &
             "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }

  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (!is.factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))) {
        p$data$node_color_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
      } else {
        p$data$node_color_as <-
          dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      }
      # double check levels are alphabetical (if not numeric)
      if (suppressWarnings(any(is.na(as.integer(levels(p$data$node_color_as)))))) {
        levels(p$data$node_color_as) <-
          sort(levels(p$data$node_color_as))
      }

    }
    if (node_color_as == "node_posterior") {
      p$data$node_color_as <- as.numeric(p$data$posterior)
    }
    if (node_color_as == "state_posterior") {
      p$data$node_color_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }

  if (is.null(node_size_as) == FALSE) {
    if (node_size_as == "state") {
      size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      if (is.factor(size_tmp)) {
        p$data$node_size_as <- as.integer(levels(size_tmp))[size_tmp]
      } else {
        p$data$node_size_as <- as.integer(size_tmp)
      }
    }
    if (node_size_as == "node_posterior") {
      p$data$node_size_as <- as.numeric(p$data$posterior)
    }
    if (node_size_as == "state_posterior") {
      p$data$node_size_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }

  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as == "state") {

      p$data$node_shape_as <-
        factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
    }
    if (node_shape_as == "node_posterior") {
      p$data$node_shape_as <- as.numeric(p$data$posterior)
    }
    if (node_shape_as == "state_posterior") {
      p$data$node_shape_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }

  if (cladogenetic == TRUE) {
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state") {
        if (!is.factor(dplyr::pull(p$data,
                                   paste0(state_pos_str_base[2], "1")))) {
          p$data$clado_node_color_as <-
            factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
        } else {
          p$data$clado_node_color_as <-
            dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        }

        if (suppressWarnings(any(is.na(as.integer(levels(p$data$clado_node_color_as)))))) {
          levels(p$data$clado_node_color_as) <-
            sort(levels(p$data$clado_node_color_as))
        }
      }
      if (node_color_as == "node_posterior") {
        p$data$clado_node_color_as <- 1
      }
      if (node_color_as == "state_posterior") {
        p$data$clado_node_color_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }

    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state") {
        clado_size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        if (is.factor(clado_size_tmp)) {
          p$data$clado_node_size_as <- as.integer(levels(clado_size_tmp))[clado_size_tmp]
        } else {
          p$data$clado_node_size_as <- as.integer(clado_size_tmp)
        }
      }
      if (node_size_as == "node_posterior") {
        p$data$clado_node_size_as <- 1
      }
      if (node_size_as == "state_posterior") {
        p$data$clado_node_size_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }

    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state") {
        p$data$clado_node_shape_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
      }
      if (node_shape_as == "node_posterior") {
        p$data$clado_node_shape_as <- 1
      }
      if (node_shape_as == "state_posterior") {
        p$data$clado_node_shape_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
  }

  # gather list of all character states from data
  if (cladogenetic == TRUE) {
    all_states <- unique(c(p$data$start_state_1, p$data$end_state_1))
  } else {
    all_states <-
      na.omit(unique(factor(dplyr::pull(
        p$data, paste0(state_pos_str_base[1], "1")
      ))))
  }
  all_states <- sort(all_states)

  ##### color processing and checks #####
  # check if number of states exceeds default color palette options
  if (!is.null(node_color_as) && node_color_as == "states") {
    if (node_color[1] == "default") {
      nstates <- length(all_states)
      if (nstates <= 12) {
        node_color <- colFun(nstates)
      } else {
        node_color <- grDevices::colorRampPalette(colFun(12))(nstates)
      }
    }

    # check if number of states not equal to provided colors
    if (node_color[1] != "default" &
        length(node_color) < length(all_states)) {
      stop(
        paste0(
          "You provided fewer colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
    if (node_color[1] != "default" &
        length(node_color) > length(all_states)) {
      stop(
        paste0(
          "You provided more colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
  }

  # set default colors
  if (any(node_color == "default")) {
    if (is.null(node_color_as) == TRUE) {
      colors <- colFun(1)
    } else if (node_color_as == "state") {
      nstates <- length(all_states)
      colors <- colFun(nstates)
      # name colors if unnamed
      names(colors) <- sort(all_states)
    } else if (node_color_as == "node_posterior" |
               node_color_as == "state_posterior") {
      colors <- colFun(2)
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
  if (tip_labels_remove_underscore) {
    p$data$label <- gsub("_", " ", p$data$label)
  }

  ##### get hjust values #####
  if (node_labels_centered) {
    hjust_node <- 0.5
    hjust_shoulder <- 0.5
  }
  if (!node_labels_centered) {
    hjust_node <- 0
    hjust_shoulder <- 1
  }
  ##### calculate cladogenetic plotting data #####
  if (cladogenetic == TRUE) {
    x <- ggtree::fortify(tree)$x
    y <- ggtree::fortify(tree)$y
    x_anc <- numeric(n_node)
    node_index <- numeric(n_node)
    for (i in 1:n_node) {
      if (.getParent(tree, i) != 0) {
        # if not the root, get the x coordinate for the parent node
        x_anc[i] <- x[.getParent(tree, i)]
        node_index[i] <- i
      }
    }
    shoulder_data <-
      data.frame(node = node_index,
                 x_anc = x_anc,
                 y = y)
    if (timeline == TRUE) {
      shoulder_data$x_anc <- shoulder_data$x_anc - tree_height
    }
    `%<+%` <- ggtree::`%<+%`
    p <- p %<+% shoulder_data
  }

  ##### start plotting #####

  # add timeline
  if (timeline == TRUE) {
    max_age <- tree_height

    if (max_age > 100) {
      interval <- 50
    } else {
      interval <- 10
    }
    dx <- max_age %% interval
    # set coordinates
    ### Fix the xlim and ylims - if no error bars, should be a function of
    ### max age and n nodes, respectively.
    ### If error bars, -x lim should be as old as the max of the error bar
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
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.1, tree_height /
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
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.05, tree_height /
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
    p <- p + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                         limits = c(-tree_height, tree_height /
                                                      2))
    p <- ggtree::revts(p)
    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )

    p <-
      p + ggplot2::geom_segment(ggplot2::aes(
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
          x_pos <- -rev(c(0, deeptime::getScaleData("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::getScaleData("periods")$max_age))
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
        p <- gginnards::append_layers(p, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      p <-
        p + ggplot2::theme(axis.title.x =
                             ggplot2::element_text(hjust =
                                                     max_age /  (2 * tot)))
    }
  }

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
    }
    else if (tip_labels_formatted == TRUE ) {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        parse = TRUE
      )
    } else {
      p <-
        p + ggtree::geom_tiplab(size = tip_labels_size,
                                offset = tip_labels_offset)
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
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state"))  {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(colour = node_color_as),
          size = tip_states_size,
          alpha = state_transparency,
          shape = tip_states_shape
        )
      }
    }

    # vary tip symbols by shape only
    # when shape is state, color is not state, size is not state
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as),
          size = tip_states_size,
          alpha = state_transparency,
          color = colors
        )
      }
    }

    # vary tip symbol by shape and color
    # when shape is state, color is state, and size is anything but state
    if (is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      if (node_color_as == "state" &
          node_shape_as == "state" &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <-  p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as,
                       color = node_color_as),
          size = tip_states_size,
          alpha = state_transparency
        )
      }
    }

    # vary tip symbol by size only
    # when size is state, color is not state, and shape is null
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as),
          shape = tip_states_shape,
          alpha = state_transparency,
          color = "grey"
        )
      }
    }

    # vary tip symbol by size and color
    # when size is state, color is state or PP, and shape is null
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE) {
      if (node_size_as == "state" &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as,
                       color = node_color_as),
          shape = tip_states_shape,
          alpha = state_transparency
        )
      }
    }
  }
  
  # plot symbols at nodes and shoulders
  blank_nodes <-
    is.null(node_color_as) == TRUE &
    is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE
  if (blank_nodes == FALSE) {
    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(
        ggplot2::aes(
          colour = node_color_as,
          size = node_size_as,
          shape = node_shape_as
        ),
        alpha = state_transparency
      )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }


    # plot if color and size vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, size = node_size_as),
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }

    #plot if color and shape vary
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, shape = node_shape_as),
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }

    #plot if size and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as, size = node_size_as),
          color = colors,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }

    #plot if just color varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as),
          size = node_size,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }

    #plot if just size varies
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(size = node_size_as),
          color = colors,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }

    #plot if just shape varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as),
          color = colors,
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
  }
  
  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = start_state_1,
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = anc_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = .convertAndRound(start_state_1_pp),
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(anc_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "node_posterior") {
      p <-
        p + ggtree::geom_nodelab(
          ggplot2::aes(label = .convertAndRound(posterior)),
          hjust = hjust_node,
          nudge_x = node_labels_offset,
          size = node_labels_size
        )
    }
  }

  # add tip states labels (text)
  if (tip_labels_states == TRUE) {
    if (state_pos_str_base[1] == "anc_state_") {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = anc_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size
        )
    } else {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = end_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size
        )
    }
  }

  # add custom colors, shapes, and sizes
  if (is.null(node_size_as) == FALSE) {
    p <-
      p + ggplot2::scale_size(range = node_size,
                              name = .titleFormat(node_size_as))
  }
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (is.null(names(colors))) {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as)
        )
      } else {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as),
          breaks = names(colors)
        )
      }

    } else if (node_color_as == "state_posterior" |
               node_color_as == "node_posterior") {
      if (cladogenetic) {
        prettify <- c(p$data$node_color_as, p$data$clado_node_color_as)
      } else {
        prettify <- p$data$node_color_as
      }
      p <- p + ggplot2::scale_color_gradient(
        low = colors[1],
        high = colors[2],
        breaks = pretty(prettify),
        name = .titleFormat(node_color_as)
      )
    }
  }
  if (is.null(node_shape_as) == FALSE) {
    p <-
      p + ggplot2::scale_shape_manual(values = node_shape,
                                      name = .titleFormat(node_shape_as))
  }

  # add space on x axis for tip labels
  if (tip_labels == TRUE) {
    if (timeline == FALSE) {
      p <- p + ggtree::xlim(0, tree_height + tree_height / 2)
    }
  }

  return(p)
}
