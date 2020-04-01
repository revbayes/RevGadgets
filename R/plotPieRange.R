## plotPieRange
# Description:
#     Plots character states as pie charts on nodes, shoulders and tips
#     Area of pie for each state is proportional to the posterior probability of that state
#    Tips are only a single colour as character states for tips are observed.

## TODO
# - Arguments list
# - Unit test

require("rvcheck")
require("ggimage")
require("tibble")
require("ggplot2")
require("ggtree")

plotPieRange <- function(t,
                         show_state_legend = T,
                         state_colors,
                         state_labels,
                         tip_label_size = 5,
                         tip_label_offset = 0,
                         tip_pie_diameter = 3,
                         node_pie_diameter = 3,
                         shoulder_pie_diameter = 3,
                         pie_nudge_x = 0,
                         pie_nudge_y = 0,
                         include_start_states = T,
                         alpha = 0.5,
                         tree_layout = "rectangular"){
  
  # Error checking
  if (!("start_state_1" %in% colnames(attributes(t)$data))) {
    stop("Start states not found in input tree.")
  }
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
  # set the root's start state to NA
  #attributes(t)$data$start_state_1[n_node] = NA
  
  p = p + ggtree::geom_tippoint(ggtree::aes(colour=factor(end_state_1)), size=1e-2, alpha=0.5)
  
  # plot invisible node states (for legend)
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)
  
  # set up the legend
  if (show_state_legend) {
    p = p + ggplot2::guides(colour=guide_legend("State"), order=1)
  } else {
    p = p + ggplot2::guides(colour=FALSE, order=2)
  }
  p = p + ggplot2::guides(size=FALSE)
  p = p + ggplot2::guides(colour = guide_legend(override.aes = list(size=5)))
  
  used_states = collect_probable_states(p)
  p = p + ggplot2::scale_color_manual(values=state_colors, breaks=state_labels,  name="Range", limits = used_states)
  p = p + theme(legend.position="left")
  
  # get anc state matrices (for pie/bar charts)
  state_probs <- build_state_probs(t, state_labels, include_start_states)
  dat_state_end = state_probs$end
  dat_state_start = state_probs$start
  
  # make pie objects
  n_tips = length(tree$tip.label)
  n_nodes = 2 * n_tips - 1
  node_idx = (n_tips+1):n_nodes
  tip_idx = 1:n_tips
  all_idx = 1:n_nodes
  
  pies_end = ggtree:::nodepie(dat_state_end,cols=1:(ncol(dat_state_end)-1),
                              color=state_colors,alpha=alpha)
  pies_start = ggtree:::nodepie(dat_state_start,cols=1:(ncol(dat_state_start)-1),
                                color=state_colors,alpha=alpha)
  
  pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )
  
  # plot
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
                       height=shoulder_pie_diameter,
                       width=shoulder_pie_diameter,
                       hjust=pie_nudge_x,
                       vjust=pie_nudge_y)
  
  
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
  # print("There are the x-coordinates (for nodes) :")
  # print(xx)
  # print("There are the y-coordinates (for nodes) :")
  # print(yy)
  if (length(width)==1) width = rep(width, length(insets))
  if (length(height)==1) height = rep(height, length(insets))
  
  # add subviews
  # old way
  # tree_view = tree_view + geom_subview_revgadgets(subview = insets, width = width,
  #                                      height = height, x = xx, y = yy)
  
  # new way using geom_inset
  tree_view = tree_view + ggtree::geom_inset(insets = insets,
                                             width = width, 
                                             height = height,
                                             x)
  
  
  # return treeview with subviews
  return(tree_view)
}

geom_subview_revgadgets <- function (mapping = NULL, data = NULL, width = 0.1, height = 0.1,
                                     x = NULL, y = NULL, subview = NULL)
  # This is basically just a copy of ggimage:::geom_subview with print statements
{
  if (is.null(data)) {
    data <- tibble(x = x, y = y)
  }
  else if (!inherits(data, "tbl")) {
    data <- as_tibble(data)
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
  # print("These are the x-min coordinates:")
  # print(data$xmin)
  # print("These are the x-max coordinates:")
  # print(data$xmax)
  # print("These are the y-min coordinates:")
  # print(data$ymin)
  # print("These are the y-max coordinates:")
  # print(data$ymax)
  lapply(1:nrow(data), function(i) {
    annotation_custom(as.grob(data$subview[[i]]), xmin = data$xmin[i],
                      xmax = data$xmax[i], ymin = data$ymin[i], ymax = data$ymax[i])
  })
}


# For testing
# t <- treeio::read.beast("simple.ase.tre")
# state_info <- read.csv("simple.state_info.txt")
# 
# PieRange_plot <- plotPieRange(t,
#                               state_labels = state_info$state,
#                               state_colors = c(state_info$color, "grey20"),
#                               pie_nudge_x = 0.01,
#                               pie_nudge_y= 0.1)
# ggsave(PieRange_plot, filename = "PieRange_plot.pdf", width = 7, height = 7)
