#' Plot discrete ancestral states on tree
#'
#' Function to plot ancestral states and the associated uncertainty
#' for continuous and discrete characters. Based on ggtree plotting
#' functionality
#'
#' For discrete characters 'summary_statistic="MAP"' should be used,
#' and for continuous characters 'summary_statistic="mean"'. If the
#' tree and tip labels do not fit in the screen, adjust the visible
#' area using the 'xlim_visible' argument.
#'
#' If 'summary_statistic="MAP"', the maximum a posteriori ancestral
#' state will be plotted on each node. The color corresponds to the
#' character state, and the size of the circle represents the posterior
#' probability of that state. Cladogenetic models that estimate
#' ancestral states for both the beginning and end of each branch
#' are plotted by setting "include_start_states=TRUE".
#'
#' Maximum a posteriori ancestral chromosome numbers can be plotted
#' with 'summary_statistic="MAPChromosome"'. For chromosomes the
#' color represents the posterior probability and the size of the
#' circle represents the chromosome number.
#'
#' For 'summary_statistic="mean"' the color represents the size of
#' the 95% confidence interval, and the size of the cirlce represents
#' represents the mean character state.
#'
#' @param t (treedata; no default) The processed output from processAncStatesDiscrete
#' @param summary_statistic (character; no default) The
#' @param tree_layout (character; "rectangular) How tree is represented (see ggtree doc).
#' @param include_start_states (logical; FALSE) Whether start states for each anagenetic branch is to be plotted on the shoulders. Only applicable for cladogenetic models where state changes may occur at nodes
#' @param xlim_visible (numeric; NULL) Plotting limits of the x-axis in cartesian space. If specified, needs to be a vector of length 2.
#' @param ylim_visible (numeric; NULL) Plotting limits of the y-axis in cartesian space. If specified, needs to be a vector of length 2.
#' @param tip_label_size (numeric; 4) Font size of tip labels
#' @param tip_label_offset (numeric; 5) Horizontal offset of tip labels from end of tips
#' @param tip_label_italics (logical; FALSE) Whether or not tip labels should be italicized
#' @param tip_node_size (numeric; 2) Size of tip points; parsed to geom_tippoint()
#' @param tip_node_shape (numeric; 15) Shape of tip points; parsed to geom_tippoint()
#' @param node_label_size (numeric; 4) Font size of node labels
#' @param node_pp_label_size (numeric; 0) Font size of node labels for posterior probabilities
#' @param node_label_nudge_x (numeric; 0.1) Degree of horizontal offset of node labels.
#' @param node_pp_label_nudge_x (numeric; 0.1) Degree of horizontal offet of node posterior probability labels
#' @param shoulder_label_size (numeric; 3) Font size of labels at node shoulders
#' @param shoulder_label_nudge_x (numeric; 3) Degree of horizontal offset of labels at node shoulders
#' @param node_pie_diameter (numeric; 1.10) Size of pies at nodes (shoulders pies are 0.9x in size)
#' @param tip_pie_diameter (numeric; 1.08) Size of pie at tips
#' @param pie_nudge_x (numeric; 0.0) Degree of horizontal offset of pie charts
#' @param pie_nudge_y (numeric; 0.0) Degree of vertical offset of pie charts
#' @param alpha (numeric; 0.5) Degree of transparency in pies and points
#' @param node_size_range (numeric; c(6,15)) Size range of points (used in plotMAP)
#' @param color_low (character;"#D55E00")
#' @param color_mid (character; "#F0E442")
#' @param color_high (character; "#009E73")
#' @param show_state_legend (logical, TRUE),
#' @param show_posterior_legend (logical, TRUE)
#' @param show_tree_scale (logical, TRUE)
#' @param title (character, NULL) Title of plot
#' @param state_colors (character, NULL) Vector of colours to be used for states. Must be of equal length to the number of states.
#' @param ... (various)
#' @return 
#' @examples
#' @export

# libraries
#require(colorspace)
#require(RColorBrewer)
#require(ggplot2)
#require(ggtree)

plotAncStatesDiscrete = function(t,
                                 summary_statistic,
                                 tree_layout="rectangular",
                                 include_start_states=FALSE,
                                 xlim_visible=c(0, 40),
                                 ylim_visible=NULL,
                                 tip_label_size=4,
                                 tip_label_offset=5,
                                 tip_label_italics=FALSE,
                                 tip_node_size=2,
                                 tip_node_shape=15,
                                 node_label_size=4,
                                 node_pp_label_size=0,
                                 node_label_nudge_x=0.1,
                                 node_pp_label_nudge_x=0.1,
                                 shoulder_label_size=3,
                                 shoulder_label_nudge_x=-0.1,
                                 node_pie_diameter=1.10,
                                 tip_pie_diameter=1.08,
                                 pie_nudge_x=0.0,
                                 pie_nudge_y=0.0,
                                 alpha=0.5,
                                 node_size_range=c(6, 15),
                                 color_low="#D55E00",
                                 color_mid="#F0E442",
                                 color_high="#009E73",
                                 show_state_legend=TRUE,
                                 show_posterior_legend=TRUE,
                                 show_tree_scale=TRUE,
                                 state_labels=NULL,
                                 state_colors=NULL,
                                 title="",
                                 fig_height=7,
                                 fig_width=7,
                                 ...) {

  # Argument checking
  if ( !(summary_statistic %in% c("MAP", "mean", "MAPChromosome", "MAPRange", "PieRange", "PieState")) ) {
    stop("Invalid summary statistic specified.")
  }

  # Some checks (needs to be done for all scenarios?)
  if(summary_statistic == "PieRange"){
    # Checking to see if PieRange data has start and end states
    if(sum(grepl(pattern = "start|end", names(attributes(t)$data))) == 0){
      stop("PieRange summary statistic chosen, but input data does not have start or end states. If you are doing ancestral state reconstruction, perhaps you would like to try PieState instead")
    }
  }
  if(summary_statistic == "PieState"){
    if(sum(grepl(pattern = "anc", names(attributes(t)$data))) == 0){
      stop("PieState summary statistic chosen, but input data does not contain ancestral states.")
    }
  }

  if ( summary_statistic %in% c("MAPRange", "PieRange", "PieState") ) {

    # State labels (should have been either provided to or generated by processAncStatesDiscrete)
    if (is.null(state_labels)) {
      state_labels <- attributes(t)$state_labels
    }

    # State colours
    if( is.null(state_colors) ) {
      print("State colors not provided by user. Defaults will be used")
      state_colors <- RColorBrewer::brewer.pal(n = length(state_labels), name = "Set3")
      names(state_colors) <- state_labels
    }
  }

  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)

  # General plotting of tree with tip labels ===========
  p = ggtree::ggtree(t, layout=tree_layout, ladderize=TRUE)

  # MJL: need to re-enable tip_label_italics if desired
  # p = p + ggtree::geom_tiplab(size=tip_label_size, offset=tip_label_offset, parse=tip_label_italics)
  p = p + ggtree::geom_tiplab(size=tip_label_size, offset=tip_label_offset)

  if (summary_statistic == "MAPChromosome") {

    p <- plotMAPchromosome(p, t, include_start_states, shoulder_label_nudge_x, shoulder_label_size, node_label_nudge_x, node_label_size, alpha, show_state_legend, show_posterior_legend, color_low, color_mid, color_high, n_node)

  }

  else if (summary_statistic == "MAPRange") {
    p <- plotMAPrange(p, t, include_start_states, tip_node_size, alpha, show_state_legend, show_posterior_legend)

  }

  else if (summary_statistic == "MAP") {
    p <- plotMAP(p, t,include_start_states, node_label_nudge_x, node_label_size, node_size_range, alpha, show_state_legend, show_posterior_legend, node_pp_label_size, tip_node_size, tip_node_shape)

  }

  else if (summary_statistic == "mean") {

    p <- plotMean(p, t, include_start_states, node_label_nudge_x, node_label_size, color_low, color_mid, color_high, alpha, show_state_legend, show_posterior_legend)

  }

  else if (summary_statistic == "PieState"){
    
    p <- plotPieState(p, t, include_start_states, show_state_legend, state_colors, state_labels, alpha, tip_pie_diameter, node_pie_diameter, pie_nudge_y, pie_nudge_x)

  }

  else if (summary_statistic == "PieRange") {

    p <- plotPieRange(p, t, show_state_legend, state_colors, state_labels, tip_pie_diameter, node_pie_diameter, pie_nudge_x, pie_nudge_y, include_start_states)
  }

  # if (use_state_colors) {
  #   #print(state_colors)
  #   #print(state_labels)
  #   p = p + ggplot2::scale_color_manual(values=state_colors, breaks=as.vector(state_labels))
  # }

  p = p + ggplot2::scale_radius(range = node_size_range)
  p = p + ggplot2::theme(legend.position="left")

  # show title
  p = p + ggplot2::ggtitle(title)

  # set visible area
  #p = p + ggplot2::coord_cartesian(xlim = xlim_visible, ylim=ylim_visible, expand=TRUE)

  return(p)
}


# modified from inset
inset.revgadgets = function (tree_view, insets, width = 0.1, height = 0.1, hjust = 0,
                             vjust = 0, x = "node", pos = 0.5){
  df <- tree_view$data[as.numeric(names(insets)), ]

  # position subviews based on tree part
  x <- match.arg(x, c("node", "branch", "edge", "parent_shoulder"))
  if (x == "node") {
    xx <- df$x
  }
  else if (x == "parent_shoulder") {
    xx <- df$x[ match(df$parent, df$node) ]
  }
  else {
    xx <- df$branch
  }
  yy <- df$y
  xx <- xx - hjust
  yy <- yy - vjust
  print("There are the x-coordinates (for nodes) :")
  print(xx)
  print("There are the y-coordinates (for nodes) :")
  print(yy)
  if (length(width)==1) width = rep(width, length(insets))
  if (length(height)==1) height = rep(height, length(insets))

  # add subviews
  tree_view = tree_view + geom_subview_revgadgets(subview = insets, width = width,
                                       height = height, x = xx, y = yy)
    #ggimage:::geom_subview(subview = insets, width = width,
                           #                      height = height, x = xx, y = yy)

  # return treeview with subviews
  return(tree_view)
}
require("rvcheck")
require("ggimage")
geom_subview_revgadgets <- function (mapping = NULL, data = NULL, width = 0.1, height = 0.1, 
                                     x = NULL, y = NULL, subview = NULL) 
  # This is basically just a copy of ggimage:::geom_subview with print statements
{
  if (is.null(data)) {
    data <- data_frame(x = x, y = y)
  }
  else if (!inherits(data, "tbl")) {
    data <- as_data_frame(data)
  }
  if (is.null(mapping)) {
    mapping <- aes_(x = ~x, y = ~y)
  }
  mapping <- as.list(mapping)
  if (is.null(mapping$x)) {
    stop("x aesthetic mapping should be provided")
  }
  if (is.null(mapping$y)) {
    stop("y aesthetic mapping should be provided")
  }
  if (is.null(mapping$subview) && is.null(subview)) {
    stop("subview must be provided")
  }
  if (is.null(mapping$subview)) {
    if (!inherits(subview, "list")) {
      subview <- list(subview)
    }
    data$subview <- subview
  }
  else {
    sv_var <- get_aes_var(mapping, "subview")
    data$subview <- data[[sv_var]]
  }
  xvar <- get_aes_var(mapping, "x")
  yvar <- get_aes_var(mapping, "y")
  if (is.null(mapping$width)) {
    data$width <- width
  }
  else {
    width_var <- get_aes_var(mapping, "width")
    data$width <- data[[width_var]]
  }
  if (is.null(mapping$height)) {
    data$height <- height
  }
  else {
    height_var <- get_aes_var(mapping, "height")
    data$height <- data[[height_var]]
  }
  data$xmin <- data[[xvar]] - data$width/2
  data$xmax <- data[[xvar]] + data$width/2
  data$ymin <- data[[yvar]] - data$height/2
  data$ymax <- data[[yvar]] + data$height/2
  print("These are the x-min coordinates:")
  print(data$xmin)
  print("These are the x-max coordinates:")
  print(data$xmax)
  print("These are the y-min coordinates:")
  print(data$ymin)
  print("These are the y-max coordinates:")
  print(data$ymax)
  lapply(1:nrow(data), function(i) {
    annotation_custom(as.grob(data$subview[[i]]), xmin = data$xmin[i], 
                      xmax = data$xmax[i], ymin = data$ymin[i], ymax = data$ymax[i])
  })
}

plotMAPchromosome <- function(p, t, include_start_states, shoulder_label_nudge_x, shoulder_label_size, node_label_nudge_x, node_label_size, alpha, show_state_legend, show_posterior_legend, color_low, color_mid, color_high, n_node){
  if (include_start_states) {

    if (!("start_state_1" %in% colnames(attributes(t)$data))) {
      print("Start states not found in input tree.")
      return()
    }

    # set the root's start state to NA
    attributes(t)$data$start_state_1[n_node] = NA

    # add clado daughter lineage start states on "shoulders" of tree
    # get x, y coordinates of all nodes
    tree = t@phylo
    x = getXcoord(tree)
    y = getYcoord(tree)
    x_anc = numeric(n_node)
    node_index = numeric(n_node)
    for (i in 1:n_node) {
      if (getParent(tree, i) != 0) {
        # if not the root, get the x coordinate for the parent node
        x_anc[i] = x[getParent(tree, i)]
        node_index[i] = i
      }
    }
    shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
    `%<+%` = ggtree::`%<+%`
    p = p %<+% shoulder_data

    # plot the states on the "shoulders"
    p = p + ggtree::geom_text(ggtree::aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)

    # add ancestral states as node labels
    p = p + ggtree::geom_text(ggtree::aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

    # show ancestral states as size / posteriors as color
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=as.numeric(end_state_1_pp), size=as.numeric(end_state_1)), alpha=alpha)

  } else {

    # add ancestral states as node labels
    p = p + ggtree::geom_text(ggtree::aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

    # show ancestral states as size / posteriors as color
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=as.numeric(end_state_1_pp), size=as.numeric(end_state_1)), alpha=alpha)

  }

  min_low = 0.0
  max_up = 1.0
  p = p + ggplot2::scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=0.5)
  if (show_state_legend) {
    p = p + ggplot2::guides(size=ggplot2::guide_legend("Chromosome Number"))
  } else {
    p = p + ggplot2::guides(size=FALSE)
  }
  if (show_posterior_legend) {
    p = p + ggplot2::guides(colour=ggplot2::guide_legend("Posterior Probability", override.aes = list(size=8)))
  } else {
    p = p + ggplot2::guides(colour=FALSE)
  }
}


plotMAPrange <- function(p, t, include_start_states, tip_node_size, alpha, show_state_legend, show_posterior_legend){
  # Plots maximum a posteriori states
  if (!include_start_states) {
    warning("Ignoring that include_start_states is set to FALSE")
  }
  if (!("start_state_1" %in% colnames(attributes(t)$data))) {
    stop("Start states not found in input tree.")
  }

  # add ancestral states as node labels
  #p = p + ggtree::geom_text(ggtree::aes(label=end_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

  # set the root's start state to NA
  attributes(t)$data$start_state_1[n_node] = NA

  # add clado daughter lineage start states on "shoulders" of tree
  # get x, y coordinates of all nodes
  x = getXcoord(tree)
  y = getYcoord(tree)
  x_anc = numeric(n_node)
  node_index = numeric(n_node)
  for (i in 1:n_node) {
    if (getParent(tree, i) != 0) {
      # if not the root, get the x coordinate for the parent node
      x_anc[i] = x[getParent(tree, i)]
      node_index[i] = i
    }
  }
  shoulder_data = data.frame(node=node_index, x_anc=x_anc, y=y)
  `%<+%` = ggtree::`%<+%`
  p = p %<+% shoulder_data

  # plot the states on the "shoulders"
  p = p + geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
  p = p + geom_nodepoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)
  p = p + geom_tippoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)

  # show tip states as color
  #print(shoulder_data)
  #print(x_anc)
  #print(c(attributes(t)$data$start_state_1,attributes(t)$data$end_state_1))

  p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha)

  # show ancestral states as color / posteriors as size
  p = p + geom_nodepoint(aes(colour=factor(end_state_1), size=as.numeric(end_state_1_pp)), alpha=alpha)

  if (show_state_legend) {
    p = p + guides(colour=guide_legend("Range", override.aes = list(size=8), order=1))
  } else {
    p = p + guides(colour=FALSE)
  }

  if (show_posterior_legend) {
    p = p + guides(size=guide_legend("Posterior probability", order=2))
  } else {
    p = p + guides(size=FALSE)
  }
  return(p)
}


plotMAP <- function(p, t, include_start_states, node_label_nudge_x, node_label_size, node_size_range, alpha, show_state_legend, show_posterior_legend, node_pp_label_size, tip_node_size, tip_node_shape) {
  if (include_start_states) {
    stop("Start states not yet implemented for MAP ancestral states.")
  }
  if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
    anc_data = data.frame(node=names(attributes(t)$data$end_state_1),
                          anc_state_1=levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
                          anc_state_1_pp=as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
    `%<+%` = ggtree::`%<+%`
    p = p %<+% anc_data
  }

  # add ancestral states as node labels
  p = p + ggtree::geom_text(ggtree::aes(label=anc_state_1), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

  # show ancestral states as color / posteriors as size
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_1), size=as.numeric(anc_state_1_pp)), alpha=alpha)

  pp = as.numeric( as.vector( attributes(t)$data$anc_state_1_pp) )

  if (!F) {
    pp_offset_range = 2*(c(min(pp), max(pp)) - 0.5)
    nd_offset_interval = node_size_range[2] - node_size_range[1]
    nd_offset = node_size_range[1]
    node_size_range = pp_offset_range * nd_offset_interval + nd_offset
    #node_size_range[1] = node_size_range[1] * min(pp) / 0.5
    #node_size_range[2] = node_size_range[2] * max(pp)
  }

  if (node_label_size == 0) {
    p = p + ggtree::geom_text(ggtree::aes(label=sprintf("%.02f", as.numeric(anc_state_1_pp))), hjust="left", nudge_x=node_label_nudge_x, size=node_pp_label_size)
  }
  #p = p = ggplot2::scale_fill_continuous(breaks=c(0.6, 0.7, 0.8, 0.9, 1.0))

  # show the tip values
  p = p + ggtree::geom_tippoint(ggtree::aes(colour=factor(anc_state_1)), size=tip_node_size, alpha=alpha, shape=tip_node_shape)

  # set up the legend
  if (show_state_legend) {
    p = p + ggplot2::guides(colour=ggplot2::guide_legend("State"), order=1)
  } else {
    p = p + guides(colour=FALSE, order=2)
  }
  if (show_posterior_legend) {
    p = p + ggplot2::guides(size=ggplot2::guide_legend("Posterior Probability"), order=3)
  } else {
    p = p + ggplot2::guides(size=FALSE, order=4)
  }
  return(p)
}


plotMean <- function(p, t, include_start_states, node_label_nudge_x, node_label_size, color_low, color_mid, color_high, alpha, show_state_legend, show_posterior_legend){
  if (include_start_states) {
    print("Start states not implemented for mean ancestral states.")
    return()
  }

  # add ancestral states as node labels
  p = p + geom_text(aes(label=round(mean, 2)), hjust="left", nudge_x=node_label_nudge_x, size=node_label_size)

  # show the size of the 95% CI as color
  lowers = as.numeric(levels(attributes(t)$data$lower_0.95_CI))[attributes(t)$data$lower_0.95_CI]
  uppers = as.numeric(levels(attributes(t)$data$upper_0.95_CI))[attributes(t)$data$upper_0.95_CI]
  diffs = uppers - lowers
  diffs_df = data.frame(node=names(attributes(t)$data$lower_0.95_CI), diff_vals=diffs)
  `%<+%` = ggtree::`%<+%`
  p = p %<+% diffs_df

  min_low = min(diffs, na.rm=TRUE)
  max_up = max(diffs, na.rm=TRUE)
  mid_val = min_low + (max_up - min_low) / 2.0
  p = p + scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high, limits=c(min_low, max_up), midpoint=mid_val)
  p = p + geom_nodepoint(aes(size=mean, colour=diff_vals), alpha=alpha)

  # show the tip values
  p = p + geom_tippoint(aes(size=mean), color="grey", alpha=alpha)

  # set up the legend
  if (show_state_legend) {
    legend_text = "Mean State"
    p = p + guides(size=guide_legend(legend_text))
  } else {
    p = p + guides(size=FALSE)
  }
  if (show_posterior_legend) {
    p = p + guides(colour=guide_legend("95% CI Width", override.aes=list(size=4)))
  } else {
    p = p + guides(colour=FALSE)
  }
  return(p)
}

plotPieState <- function(p, t, include_start_states, show_state_legend, state_colors, state_labels, alpha, tip_pie_diameter, node_pie_diameter, pie_nudge_y, pie_nudge_x){
  
  if (include_start_states) {
    print("Start states not yet implemented for PieState ancestral states.")
    return()
  }
  
  if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
    anc_data = data.frame(node=names(attributes(t)$data$end_state_1),
                          anc_state_1=levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
                          anc_state_1_pp=as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
    `%<+%` = ggtree::`%<+%`
    p = p %<+% anc_data
  }
  
  # print tips
  #p = p + ggtree::geom_tippoint(ggtree::aes(colour=factor(anc_state_1)), size=1, alpha = 0.0)
  
  # plot invisible node states (for legend)
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_1), size=0),na.rm=TRUE, alpha=0.0)
  #p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_2), size=0),na.rm=TRUE, alpha=0.0)
  #p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_3), size=0),na.rm=TRUE, alpha=0.0)
  
  # set up the legend
  if (show_state_legend) {
    p = p + ggplot2::guides(colour=ggplot2::guide_legend("State", override.aes = list(size=length(state_labels), alpha = 1.0)), order=1)
  } else {
    p = p + ggplot2::guides(colour=FALSE, order=2)
  }
  p = p + ggplot2::guides(size=FALSE)

  #p = p + ggplot2::scale_color_manual(values=state_colors, breaks=state_labels)

  # position legend
  p = p + ggplot2::theme(legend.position="left")
  
  # get anc state matrices (for pie/bar charts)
  dat_state_anc = build_state_probs(t, state_labels, include_start_states)$anc
  
  # make pie objects
  pies_anc = ggtree::nodepie(dat_state_anc, cols=1:(ncol(dat_state_anc)-1), alpha=alpha)

  p = p + ggtree::geom_inset(pies_anc, height=node_pie_diameter)
  
  return(p)
}


plotPieRange <- function(p, t, show_state_legend, state_colors, state_labels, tip_pie_diameter, node_pie_diameter, pie_nudge_x, pie_nudge_y, include_start_states){
  if (!("start_state_1" %in% colnames(attributes(t)$data))) {
    stop("Start states not found in input tree.")
  }

  # set the root's start state to NA
  #attributes(t)$data$start_state_1[n_node] = NA

  # print tips
  p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=1e-2, alpha=alpha)

  # plot invisible node states (for legend)
  p = p + geom_nodepoint(aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
  p = p + geom_nodepoint(aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
  p = p + geom_nodepoint(aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)

  # set up the legend
  if (show_state_legend) {
    p = p + guides(colour=guide_legend("State"), order=1)
  } else {
    p = p + guides(colour=FALSE, order=2)
  }
  p = p + guides(size=FALSE)
  p = p + guides(colour = guide_legend(override.aes = list(size=5)))

  # if (use_state_colors) {
    used_states = collect_probable_states(p)
    p = p + scale_color_manual(values=state_colors, breaks=state_labels,  name="Range", limits = used_states)
  # }
  p = p + theme(legend.position="left")

  # # MJL: to remove later
  # break_legend = F
  # if (break_legend) {
  #     p$data$x = p$data$x + (15 - max(p$data$x))
  #     x_breaks = 0:15
  #     x_labels = rep("", 16)
  #     x_labels[ c(0,5,10,15)+1 ] = c(0,5,10,15)
  #     p = p + scale_x_continuous(breaks = x_breaks, labels = rev(x_labels), sec.axis = sec_axis(~ ., breaks = 15-c(6.15, 4.15, 2.55, 1.2), labels=c("+K","+O","+M","+H") ))
  #     p = p + theme_tree2()
  #     p = p + coord_cartesian(xlim = c(0,20), expand=TRUE)
  #     p = p + labs(x="Age (Ma)")
  #     p = add_island_times(p)
  #     p = p + theme(legend.position="left", axis.line = element_line(colour = "black"))
  #     p = p + guides(colour = guide_legend(override.aes = list(size=5), nrow=6))
  # }

  # get anc state matrices (for pie/bar charts)
  #print(t)
  state_probs <- build_state_probs(t, state_labels, include_start_states)
  dat_state_end = state_probs$end
  dat_state_start = state_probs$start

  # make pie objects
  n_tips = length(tree$tip.label)
  n_nodes = 2 * n_tips - 1
  node_idx = (n_tips+1):n_nodes
  tip_idx = 1:n_tips
  all_idx = 1:n_nodes

  pies_end = nodepie(dat_state_end,cols=1:(ncol(dat_state_end)-1),color=state_colors,alpha=alpha)
  pies_start = nodepie(dat_state_start,cols=1:(ncol(dat_state_start)-1),color=state_colors,alpha=alpha)

  pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )

  p_node = inset.revgadgets(tree_view=p,
                            insets=pies_end[all_idx],
                            x="node",
                            height=pd,
                            width=pd,
                            hjust=pie_nudge_x,
                            vjust=pie_nudge_y)


  p = inset.revgadgets(tree_view=p_node,
                       insets=pies_start,
                       x="parent_shoulder",
                       height=node_pie_diameter*0.9,
                       width=node_pie_diameter*0.9,
                       hjust=pie_nudge_x,
                       vjust=pie_nudge_y)
  return(p)
}
