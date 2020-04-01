# plotMAPrange
# Description:
#     Plot character state and posterior probability as points on tips, nodes, shoulders (colour = character state, size = pp)

plotMAPrange <- function(t,
                         include_start_states = T,
                         tip_label_size = 4,
                         tip_label_offset = 1,
                         tip_node_size = 2,
                         shoulder_label_size = 3,
                         shoulder_label_nudge_x = -0.1,
                         alpha = 0.5,
                         show_state_legend = T,
                         show_posterior_legend = T,
                         tree_layout = "rectangular"){
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
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
  p = p + ggtree::geom_text(aes(label=start_state_1, x=x_anc, y=y), hjust="right", nudge_x=shoulder_label_nudge_x, size=shoulder_label_size, na.rm=TRUE)
  p = p + ggtree::geom_nodepoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)
  p = p + ggtree::geom_tippoint(aes(colour=factor(start_state_1), x=x_anc, y=y, size=as.numeric(end_state_1_pp)),na.rm=TRUE, alpha=alpha)
  
  # show tip states as color
  #print(shoulder_data)
  #print(x_anc)
  #print(c(attributes(t)$data$start_state_1,attributes(t)$data$end_state_1))
  
  p = p + ggtree::geom_tippoint(ggtree::aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha)
  
  # show ancestral states as color / posteriors as size
  p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(end_state_1), size=as.numeric(end_state_1_pp)), alpha=alpha)
  
  if (show_state_legend) {
    p = p + ggplot2::guides(colour=guide_legend("Range", override.aes = list(size=8), order=1))
  } else {
    p = p + ggplot2::guides(colour=FALSE)
  }
  
  if (show_posterior_legend) {
    p = p + ggplot2::guides(size=guide_legend("Posterior probability", order=2))
  } else {
    p = p + ggplot2::guides(size=FALSE)
  }
  return(p)
}


# Testing
# t <- treeio::read.beast("chromosomes_ancestral_states.tree")
# maprange_plot <- plotMAPrange(t)
# ggsave("MAPrange_plot.pdf", maprange_plot, width = 12, height = 6)