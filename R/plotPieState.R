# plotPieState
# Description:
#    Plots character states as pie charts on nodes and tips
#    Area of pies for each state (on nodes) is proportional to the posterior probability
#    Tips are only a single colour as character states for tips are observed.

plotPieState <- function(t,
                         include_start_states = F,
                         show_state_legend = T,
                         tip_label_size = 5,
                         tip_label_offset = 1,
                         alpha = 0.5,
                         tip_pie_diameter = 2,
                         node_pie_diameter = 2,
                         pie_nudge_y = 0,
                         pie_nudge_x = 0,
                         state_colors = NULL,
                         state_labels = NULL, 
                         tree_layout = "rectangular"){
  
  # get number of nodes
  tree = attributes(t)$phylo
  n_node = ggtree:::getNodeNum(tree)
  
  # create basic tree plot
  p <- ggtree:::ggtree(t, layout = tree_layout, ladderize = TRUE)
  p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  
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
  
  p = p + ggplot2::scale_color_manual(values=state_colors, breaks=state_labels)
  
  # position legend
  p = p + ggplot2::theme(legend.position="left")
  
  # get anc state matrices (for pie/bar charts)
  dat_state_anc = build_state_probs(t, state_labels, include_start_states)$anc
  
  # make pie objects
  pies_anc = ggtree::nodepie(dat_state_anc, cols=1:(ncol(dat_state_anc)-1),
                             color=state_colors, alpha=alpha)
  
  #height and width are the size of the insets relative to the axes, so we should calculate
  #p = p + ggtree::geom_inset(pies_anc, height=node_pie_diameter, hjust=pie_nudge_x, vjust=pie_nudge_y)
  p = p + ggtree::geom_inset(pies_anc, height=0.1, hjust=pie_nudge_x, vjust=pie_nudge_y)
  #return(p)
}

# Testing
# t <- treeio::read.beast("ancestral_states_pie.tree")
# PieState_plot <- plotPieState(t, state_labels = c(0:2),
#                               state_colors = c("red", "yellow", "purple", "grey"),
#                               tip_pie_diameter = 1,
#                               node_pie_diameter = 2)
# ggsave("PieState_plot.pdf", PieState_plot, width = 10, height = 10)

