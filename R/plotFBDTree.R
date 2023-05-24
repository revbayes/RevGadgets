#' Plot FBD tree
#'
#' Plots a single FBD tree, such as an MCC or MAP tree.
#'
#' Plots a single tree, such as an MCC or MAP tree, with
#' special features for plotting the output of fossilized
#' birth-death analyses.
#'
#' @param tree (list of lists of treedata objects; no default) Name of a list
#' of lists of treedata objects, such as produced by readTrees(). This
#' object should only contain only one summary tree from one trace file. If it
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
#' @param age_bars_width (numeric; 1.5) Change line width for age bars
#'
#' @param node_labels (character; NULL) Plot text labels at nodes, specified
#' by the name of the corresponding column in the tidytree object. If NULL,
#' no text is plotted.
#'
#' @param node_labels_color (character; "black") Color to plot node_labels,
#' either as a valid R color name or a valid hex code.
#'
#' @param node_labels_size (numeric; 3) Size of node labels
#'
#' @param node_labels_offset (numeric; 0) Horizontal offset of node labels
#' from nodes.
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
#' @param tip_labels_color (character; "black") Color to plot tip labels,
#' either as a valid R color name or a valid hex code.
#'
#' @param tip_labels_size (numeric; 3) Size of tip labels
#'
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from
#' tree.
#'
#' @param node_pp (logical; FALSE) Plot posterior probabilities as symbols at
#' nodes? Specify symbol aesthetics with node_pp_shape, node_pp_color, and
#' node_pp_size.
#'
#' @param node_pp_shape (integer; 1) Integer corresponding to point shape
#' (value between 0-25). See ggplot2 documentation for details:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}
#'
#' @param node_pp_color (character; "black") Color for node_pp symbols, either
#' as valid R color name(s) or hex code(s). Can be a single character string
#' specifying a single color, or a vector of length two specifying two colors
#' to form a gradient. In this case, posterior probabilities will be indicated
#' by color along the specified gradient.
#'
#' @param node_pp_size (numeric or character; 1) Size for node_pp symbols.
#' If numeric, the size will be fixed at the specified value. If a character,
#' it should specify "variable", indicating that size should be scaled by the
#' posterior value. Size regulates the area of the shape,
#' following ggplot2 best practices:
#' \url{https://ggplot2.tidyverse.org/reference/scale_size.html})
#'
#' @param tip_age_bars (logical; FALSE) Plot node age bars for the tips as
#' well? Useful for plotting serial sampled analyses or fossilized birth-death
#' analyses, or any cases where some tip ages are estimated.
#'
#' @param branch_color (character; "black") A single character string
#' specifying the color (R color name or hex code) for all branches OR a
#' vector of length 2 specifying two colors for a gradient, used to color the
#' branches according to the variable specified in color_branch_by. If only 1
#' color is provided and you specify color_branch_by, default colors will be
#' chosen (low = "#005ac8", high = "#fa7850").
#'
#' @param color_branch_by (character; NULL ) Optional name of one quantitative
#' variable in the treedata object to color branches, such as a rate.
#'
#' @param line_width (numeric; 1) Change line width for branches
#'
#' @param label_sampled_ancs (logical; FALSE) Label any sampled ancestors?
#' Will inherent tip labels aesthetics for size and color.
#'
#' @param ... (various) Additional arguments passed to ggtree::ggtree().
#'
#' @return returns a single plot object.
#'
#' @examples
#' \donttest{
#' file <- system.file("extdata", "fbd/bears.mcc.tre", package="RevGadgets")
#' tree <- readTrees(paths = file)
#' plotFBDTree(tree = tree, timeline = TRUE, tip_labels_italics = FALSE,
#'             tip_labels_remove_underscore = TRUE,
#'             node_age_bars = TRUE, age_bars_colored_by = "posterior",
#'             age_bars_color = rev(colFun(2))) +
#'   ggplot2::theme(legend.position=c(.25, .85))
#' }
#' @export


plotFBDTree <- function(
                        tree,
                        timeline = FALSE,
                        geo = timeline,
                        geo_units = list("epochs", "periods"),
                        time_bars = timeline,

                        node_age_bars = TRUE,
                        tip_age_bars = TRUE,
                        age_bars_color = "blue",
                        age_bars_colored_by = NULL,
                        age_bars_width = 1.5,

                        node_labels = NULL,
                        node_labels_color = "black",
                        node_labels_size = 3,
                        node_labels_offset = 0,

                        tip_labels = TRUE,
                        tip_labels_italics = FALSE,
                        tip_labels_formatted = FALSE,
                        tip_labels_remove_underscore = TRUE,
                        tip_labels_color = "black",
                        tip_labels_size = 3,
                        tip_labels_offset = 0,

                        node_pp = FALSE,
                        node_pp_shape = 16,
                        node_pp_color = "black",
                        node_pp_size = "variable",

                        branch_color = "black",
                        color_branch_by = NULL,
                        line_width = 1,

                        label_sampled_ancs = FALSE,
                        ...
                        ) {

  plotTreeFull(tree = tree,

               timeline = timeline,
               geo_units = geo_units,
               geo = geo,
               time_bars = time_bars,

               node_age_bars = node_age_bars,
               tip_age_bars = tip_age_bars,
               age_bars_color = age_bars_color,
               age_bars_colored_by = age_bars_colored_by,
               age_bars_width = age_bars_width,

               node_labels = node_labels,
               node_labels_color = node_labels_color,
               node_labels_offset = node_labels_offset,
               node_labels_size = node_labels_size,

               tip_labels = tip_labels,
               tip_labels_italics = tip_labels_italics,
               tip_labels_formatted = tip_labels_formatted,
               tip_labels_remove_underscore = tip_labels_remove_underscore,
               tip_labels_color = tip_labels_color,
               tip_labels_size = tip_labels_size,
               tip_labels_offset = tip_labels_offset,

               node_pp = node_pp,
               node_pp_shape = node_pp_shape,
               node_pp_color = node_pp_color,
               node_pp_size = node_pp_size,

               branch_color = branch_color,
               color_branch_by = color_branch_by,
               line_width = line_width,

               label_sampled_ancs = label_sampled_ancs,

               # Turn off plotTree-specific aspects
               tree_layout = "rectangular",
               ...
  )
}

