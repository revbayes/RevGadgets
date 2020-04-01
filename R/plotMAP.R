# plotMAP
# Description:
#     Plot character states and posterior probabilities as points on nodes (size = pp, colour = state)
#     and character state on tips (different shape from nodes, colour = state)

plotMAP <- function(t,
                    include_start_states= F,
                    tip_label_size = 4,
                    tip_label_offset = 0.1,
                    tip_node_size = 2,
                    tip_node_shape = 15,
                    node_label_size = 4,
                    node_label_nudge_x = 0.1,
                    node_pp_label_size = 0.1,
                    node_size_range = c(6, 15),
                    alpha = 0.5,
                    show_state_legend = T,
                    show_posterior_legend = T,
                    tree_layout = "rectangular") {
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
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
    p = p + ggplot2::guides(colour=FALSE, order=2)
  }
  if (show_posterior_legend) {
    p = p + ggplot2::guides(size=ggplot2::guide_legend("Posterior Probability"), order=3)
  } else {
    p = p + ggplot2::guides(size=FALSE, order=4)
  }
  return(p)
}

# Testing
# t <- treeio::read.beast("ancestral_states.tree")
# map_plot <- plotMAP(t, tip_label_size = 2, tip_label_offset = 1, node_label_size = 2, node_pp_label_size = 2,
#                     tree_layout = "rectangular") + coord_cartesian(xlim = c(0, 20))
# map_plot <- plotMAP(t, tip_label_size = 2, tip_label_offset = 1, node_label_size = 0, node_pp_label_size = 2,
#                     tree_layout = "rectangular") + coord_cartesian(xlim = c(0, 20)) # posterior probabilities only plotted if node_label_size is 0.
# ggsave(map_plot, filename = "MAP_plot.pdf", width = 8, height = 6)
