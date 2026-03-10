#' plotStochMaps
#'
#' @param tree (treedata object; none) Output of readTrees() function
#' containing tree.
#' 
#' @param maps (dataframe; no default) Dataframe with processed maps,
#' as in the output of processStochMaps()
#' 
#' @param colors (named character vector; no default) Named character vector
#' where items are colors and names are the corresponding states. The order the 
#' vector will determine the order of the legend.
#' 
#' @param color_by (character string; "prob") How to color the branches. 
#' Options are "MAP" for assigning color by the MAP state, or "prob" for 
#' assigning color as the weighted average of states based on the posterior
#' probabilities.
#'
#' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with
#' timescale in MYA.
#'
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
#' @param tip_labels (logical; TRUE) Plot tip labels?
#'
#' @param tip_labels_italics (logical; FALSE) Plot tip labels in italics?
#'
#' @param tip_labels_formatted (logical; FALSE) Do the tip labels contain
#' manually added formatting information? Will set parse = TRUE in geom_text()
#' and associated functions to interpret formatting. See ?plotmath for more.
#' Cannot be TRUE if tip_labels_italics = TRUE.
#'
#' @param tip_labels_remove_underscore (logical; TRUE) Remove underscores in tip
#' labels?
#'
#' @param tip_labels_color (character; "black") Color to plot tip labels, either
#' as a valid R color name or a valid hex code.
#'
#' @param tip_labels_size (numeric; 3) Size of tip labels
#'
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from
#' tree.
#'
#' @param line_width (numeric; 1) Change line width for branches
#'
#' @param tree_layout (character; "rectangular") Tree shape layout, passed
#' to ggtree(). Options are 'rectangular', 'fan', 'circular', or 'inward_circular'.
#' 
#' @param label_sampled_ancs (logical; FALSE) Label any sampled ancestors?
#' Will inherent tip labels aesthetics for size and color.
#'
#' @param ... (various) Additional arguments passed to ggtree::ggtree().
#'
#' @return returns a single plot object.
#'
#' @examples
#'
#' \donttest{
#'
#' # download the example dataset to working directory
#' 
#' tree_url <- 
#' "https://revbayes.github.io/tutorials/morph_ase/data/solitariness_ase_hrm.tree"
#' tree_dest_path <- "solitariness_ase_hrm.tree"
#' download.file(tree_url, tree_dest_path)
#' 
#' maps_url <- 
#' "https://revbayes.github.io/tutorials/morph_ase/data/solitariness_hrm_stoch_char_map.log"
#' maps_dest_path <- "solitariness_hrm_stoch_char_map.log"
#' download.file(maps_url, maps_dest_path)
#' 
#' # to run on your own data, change these to the paths to your data files
#' tree_file <- tree_dest_path
#' maps_file <- maps_dest_path
#' 
#' # read in tree
#' tree <- readTrees(tree_file)
#'
#' # default, don't rename states
#' stoch_map_df <- processStochMaps(tree,
#'                                  maps_file, 
#'                                  state_labels = as.character(c(0:3)), 
#'                                  burnin = 0.1)
#' # rename states 
#' stoch_map_df_named <- processStochMaps(tree,
#'                                        maps_file, 
#'                                        state_labels = c("no - slow" = "0", "yes - slow" = "1",
#'                                                         "no - fast" = "2", "yes - fast" = "3"), 
#'                                        burnin = 0.1)
#' 
#' # rename and combine states
#' stoch_map_df_combined <- processStochMaps(tree, 
#'                                           maps_file, 
#'                                           state_labels = c("no" = "0", "yes" = "1",
#'                                                            "no" = "2", "yes" = "3"), 
#'                                           burnin = 0.1)
#' 
#' # plot by MAP with default colors
#' plotStochMaps(tree = tree,
#'               maps = stoch_map_df,
#'               color_by = "MAP",
#'               colors = "default",
#'               tip_labels = FALSE) 
#' 
#' # plot by map but with custom colors and labels, order matters
#' 
#' clrs <- c("no - slow" = "#a6cee3",
#'           "no - fast" = "#1f78b4",
#'           "yes - slow" = "#b2df8a",
#'           "yes - fast" = "#33a02c")
#' 
#' plotStochMaps(tree = tree,
#'               maps = stoch_map_df_named,
#'               color_by = "MAP",
#'               colors = clrs,
#'               tip_labels = FALSE) 
#' 
#' # plot by probability, only two states shown
#' 
#' # default
#' plotStochMaps(tree = tree,
#'               maps = stoch_map_df_combined,
#'               color_by = "prob",
#'               colors = "default",
#'               tip_labels = FALSE)  
#'
#' # custom colors               
#' clrs <- clrs <- c("no" = "#a6cee3",
#'                   "yes" = "#33a02c")
#' 
#' plotStochMaps(tree = tree,
#'               maps = stoch_map_df_combined,
#'               color_by = "prob",
#'               colors = clrs,
#'               tip_labels = FALSE) 
#'
#' # remove files
#' # WARNING: only run for example dataset!
#' # otherwise you might delete your data!
#' file.remove(tree_dest_path, maps_dest_path)
#'
#' }
#'
#' @export

plotStochMaps <- function(tree,
                          maps,
                          colors = "default",
                          color_by = "prob",
                          tree_layout = "rectangular",
                          line_width = 1,
                          tip_labels = TRUE,
                          tip_labels_italics = FALSE,
                          tip_labels_formatted = FALSE,
                          tip_labels_remove_underscore = TRUE,
                          tip_labels_color = "black",
                          tip_labels_size = 3,
                          tip_labels_offset = 0,
                          timeline = FALSE,
                          geo_units = list("epochs", "periods"),
                          geo = timeline,
                          time_bars = timeline,
                          label_sampled_ancs = FALSE,
                          ...) {
  # pull tree from list object if necessary
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {
      stop("tree should contain only one tree object")
    }
  }
  
  # we do this twice because sometimes the tree is doubly nested (inside
  # of two lists)    
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {
      stop("tree should contain only one tree object")
    }
  }
  
  # plot the base tree
  p <-  plotTreeFull(
    tree = list(list(tree)),
    tree_layout = tree_layout,
    line_width = line_width,
    
    tip_labels = tip_labels,
    tip_labels_italics = tip_labels_italics,
    tip_labels_formatted = tip_labels_formatted,
    tip_labels_remove_underscore = tip_labels_remove_underscore,
    tip_labels_color = tip_labels_color,
    tip_labels_size = tip_labels_size,
    tip_labels_offset = tip_labels_offset,
    
    timeline = timeline,
    geo_units = geo_units,
    geo = timeline,
    time_bars = timeline,
    
    label_sampled_ancs = label_sampled_ancs,
    
    node_age_bars = FALSE,
    age_bars_color = "blue",
    age_bars_colored_by = NULL,
    age_bars_width = 1,
    
    node_labels = NULL,
    node_labels_color = "black",
    node_labels_size = 3,
    node_labels_offset = 0,
    
    node_pp = FALSE,
    node_pp_shape = 16,
    node_pp_color = "black",
    node_pp_size = "variable",
    
    branch_color = "black",
    color_branch_by = NULL,
    
    tip_age_bars = FALSE,
    lineend = "square",
    ...
  )
  
  # plot the colors on branches
  if (color_by == "MAP") {
    p <- .plotStochMapsMAP(tree,
                           p,
                           maps,
                           colors,
                           tree_layout,
                           line_width,
                           tip_labels,
                           tip_labels_italics,
                           tip_labels_formatted,
                           tip_labels_remove_underscore,
                           tip_labels_color,
                           tip_labels_size,
                           tip_labels_offset,
                           timeline,
                           geo_units,
                           geo,
                           time_bars,
                           label_sampled_ancs,
                           ...)
  } else if (color_by == "prob") {
    p <- .plotStochMapsProbs(tree,
                             p,
                             maps,
                             colors,
                             tree_layout,
                             line_width,
                             tip_labels,
                             tip_labels_italics,
                             tip_labels_formatted,
                             tip_labels_remove_underscore,
                             tip_labels_color,
                             tip_labels_size,
                             tip_labels_offset,
                             timeline,
                             geo_units,
                             geo,
                             time_bars,
                             label_sampled_ancs,
                             ...)
  } else {
    stop("color_by must be either 'MAP' or 'prob'")
  }
  
  return(p)
    
}



