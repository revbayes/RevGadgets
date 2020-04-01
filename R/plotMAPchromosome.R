
plotMAPchromosome <- function(t,
                              include_start_states = FALSE,
                              show_state_legend = TRUE,
                              show_posterior_legend = TRUE,
                              tip_label_size = 4,
                              tip_label_offset  = 1,
                              node_label_nudge_x = 0.1,
                              node_label_size = 4,
                              shoulder_label_nudge_x = 0.1,
                              shoulder_label_size = 3,
                              alpha = 0.5,
                              color_low = "#D55E00",
                              color_mid = "#F0E442",
                              color_high = "#009E73",
                              tree_layout = "rectangular"){
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
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
  
  p = p + ggplot2::scale_colour_gradient2(low=color_low, mid=color_mid, high=color_high,
                                          limits=c(0.0, 1.0), midpoint=0.5) 
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
  return(p)
}


# For testing
# t <- treeio::read.beast("chromosomes_ancestral_states.tree")
# mapchromo_plot <- plotMAPchromosome(t = t, include_start_states = T)
# mapchromo_plot <- mapchromo_plot + ggplot2::coord_cartesian(xlim = c(0,15))
# ggsave("MAPchromosome_plot.pdf", mapchromo_plot, width = 12, height = 6)