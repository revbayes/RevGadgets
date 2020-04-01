# plotMean
# Description:
#    I have no idea what this does, I don't have any example data

plotMean <- function(t,
                     include_start_states = F,
                     tip_label_size = 4,
                     tip_label_offset = 1,
                     node_label_size = 4,
                     node_label_nudge_x = 0.1,
                     color_low = "#D55E00",
                     color_mid = "#F0E442",
                     color_high = "#009E73",
                     alpha = 0.5,
                     show_state_legend = T,
                     show_posterior_legend = T,
                     tree_layout = "rectangular"){
  if (include_start_states) {
    print("Start states not implemented for mean ancestral states.")
    return()
  }
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
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
    p = p + ggplot2::guides(size=guide_legend(legend_text))
  } else {
    p = p + ggplot2::guides(size=FALSE)
  }
  if (show_posterior_legend) {
    p = p + ggplot2::guides(colour=guide_legend("95% CI Width", override.aes=list(size=4)))
  } else {
    p = p + ggplot2::guides(colour=FALSE)
  }
  return(p)
}

# Testing
# t <- treeio::read.beast("ancestral_states.tree")
# map_plot <- plotMean(t) #+ coord_cartesian(xlim = c(0, 20))
# map_plot <- plotMAP(t, tip_label_size = 2, tip_label_offset = 1, node_label_size = 0, node_pp_label_size = 2,
#                     tree_layout = "rectangular") + coord_cartesian(xlim = c(0, 20)) # posterior probabilities only plotted if node_label_size is 0.
# ggsave(map_plot, filename = "MAP_plot.pdf", width = 8, height = 6)
