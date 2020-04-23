#' Plot ancestral states (discrete, pie)
#'
#' [ function tags to be written ]
#'
#' @export
plotAncStateDiscretePie <- function(t,
                                    cladogenetic = T,
                                    # show_tip_pie = T,
                                    # show_node_pie = T,
                                    # show_shoulder_pie = T,
                                    tip_pie_diameter = 1, # change to tip_diameter? we shouldn't be plotting tip states as pies, no?
                                    node_pie_diameter = 1,
                                    shoulder_pie_diameter = 1,
                                    tip_labels_remove_underscore = T,
                                    tip_labels_italics = F,
                                    pie_nudge_x = 0,
                                    pie_nudge_y = 0,
                                    show_tip_labels = F,
                                    tip_label_size = 5,
                                    tip_label_offset = 0,
                                    alpha = 0.5,
                                    state_labels,
                                    state_colors,
                                    collapse_states = T,
                                    show_legend = T) {



  # check data
  if(cladogenetic == T & any(grepl( names(attributes(t)$data), pattern = "start\\_state")) == F){
    stop("Your analysis does not include a cladogenetic component. Please specify `cladogenetic = FALSE`")
  }
  if(cladogenetic == F & any(grepl( names(attributes(t)$data), pattern = "start\\_state")) == T){
    stop("Your analysis included a cladogenetic component. Please specify `cladogenetic = TRUE`.
         If you do not wish to plot character states after cladogenesis, you may additionally specify `show_shoulder_pie = F`")
  }

  # reformat some data

  # reformat tip labels if necessary
  if (tip_labels_remove_underscore & !tip_labels_italics) {
    t@phylo$tip.label <- gsub("_", " ", t@phylo$tip.label)
  } else if (tip_labels_remove_underscore & tip_labels_italics) {
    stop("removing underscores and italicizing tip labels is not currently supported")
  }


  # get some basic tree information
  tree = attributes(t)$phylo
  n_tips = length(tree$tip.label)
  n_nodes = 2 * n_tips - 1
  node_idx = (n_tips+1):n_nodes
  tip_idx = 1:n_tips
  all_idx = 1:n_nodes

  # create basic tree plot
  p <- ggtree:::ggtree(t, ladderize = TRUE)

  # plot tip labels
  if(show_tip_labels == T){
    p <- p + ggtree:::geom_tiplab(size = tip_label_size, offset = tip_label_offset)
  }

  # old code that would coerce code from a cladogenetic model by taking the end state
  # if (!("anc_state_1" %in% colnames(attributes(t)$data))) {
  #   anc_data = data.frame(node=names(attributes(t)$data$end_state_1),
  #                         anc_state_1=levels(attributes(t)$data$end_state_1)[attributes(t)$data$end_state_1],
  #                         anc_state_1_pp=as.numeric(levels(attributes(t)$data$end_state_1_pp))[attributes(t)$data$end_state_1_pp])
  #   `%<+%` = ggtree::`%<+%`
  #   p = p %<+% anc_data
  # }

  # set up the legend
  # plot invisible node states (for legend)
  #p = p + ggtree::geom_tippoint(ggtree::aes(colour=factor(anc_state_1)), size=1, alpha = 0.0)
  if(cladogenetic == TRUE){
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)
  } else {
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_1), size=0),na.rm=TRUE, alpha=0.0)
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_2), size=0),na.rm=TRUE, alpha=0.0)
    p = p + ggtree::geom_nodepoint(ggtree::aes(colour=factor(anc_state_3), size=0),na.rm=TRUE, alpha=0.0)
  }

  if (show_legend == TRUE) {
    p = p + ggplot2::guides(colour=ggplot2::guide_legend("State", override.aes = list(size=4, alpha = 1.0)), order=1)
  } else {
    p = p + ggplot2::guides(colour=FALSE, order=2)
  }
  p = p + ggplot2::guides(size=FALSE)

  # collapse states
  if(collapse_states == TRUE){
    print(RevGadgets:::.colFun(2))
    used_states = RevGadgets:::.collect_probable_states(p)
    p = p + ggplot2::scale_color_manual(values=state_colors, breaks=state_labels,  name="Range", limits = used_states)
  } else {
    p = p + ggplot2::scale_color_manual(values=state_colors, breaks=state_labels)
  }

  p = p + ggplot2::theme(legend.position="left")


  # plot pies
  if(cladogenetic == TRUE){
    # create state matrices (matrix of nodes (rows) and all possible states (columns), values are pp. )
    state_probs <- .build_state_probs(t, state_labels, include_start_states = TRUE)
    dat_state_end = state_probs$end
    dat_state_start = state_probs$start

    # make pies
    pies_start = ggtree:::nodepie(dat_state_start,cols=1:(ncol(dat_state_start)-1),
                                  color=state_colors,alpha=alpha)
    pies_end = ggtree:::nodepie(dat_state_end,cols=1:(ncol(dat_state_end)-1),
                                color=state_colors,alpha=alpha)

    # define pie size
    pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )

    # add pies to tree
    p = inset.revgadgets(tree_view=p,
                              insets=pies_end[all_idx],
                              x="node",
                              height=pd,
                              width=pd,
                              hjust=pie_nudge_x,
                              vjust=pie_nudge_y)

    if(show_shoulder_pie == T){
      p = inset.revgadgets(tree_view=p,
                           insets=pies_start,
                           x="parent_shoulder",
                           height=shoulder_pie_diameter,
                           width=shoulder_pie_diameter,
                           hjust=pie_nudge_x,
                           vjust=pie_nudge_y)
    }

  } else {
    # create state matrices (matrix of nodes (rows) and all possible states (columns), values are pp. )
    dat_state_anc = .build_state_probs(t, state_labels, include_start_states = FALSE)$anc

    # make pies
    pies_anc = ggtree::nodepie(dat_state_anc, cols=1:(ncol(dat_state_anc)-1),
                               color=state_colors, alpha=alpha)

    # define pie size
    pd = c( rep(tip_pie_diameter, n_tips), rep(node_pie_diameter, n_nodes-n_tips) )

    # add pies to tree
    #p = p + ggtree::geom_inset(pies_anc, height=pd, hjust=pie_nudge_x, vjust=pie_nudge_y)
    p = inset.revgadgets(tree_view=p,
                         insets=pies_anc[all_idx],
                         x="node",
                         height=pd,
                         width=pd,
                         hjust=pie_nudge_x,
                         vjust=pie_nudge_y)
  }

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
   xx <- xx - hjust # x-coordinates for nodes
   yy <- yy - vjust # y-coordinates for nodes

   if (length(width)==1) width = rep(width, length(insets))
   if (length(height)==1) height = rep(height, length(insets))

   # add subviews
   # old way
   tree_view = tree_view + geom_subview_revgadgets(subview = insets, width = width,
                                        height = height, x = xx, y = yy)

   # new way using geom_inset
   # tree_view = tree_view + ggtree::geom_inset(insets = insets,
   #                                            width = width,
   #                                            height = height)
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

   lapply(1:nrow(data), function(i) {
     annotation_custom(as.grob(data$subview[[i]]), xmin = data$xmin[i],
                       xmax = data$xmax[i], ymin = data$ymin[i], ymax = data$ymax[i])
   })
 }


# # For testing
# t_clado <- treeio::read.beast("inst/extdata/dec/simple.ase.tre")
# state_info_range <- read.csv("inst/extdata/dec/simple.state_info.txt")
#
# PieRange_plot <- plotAncStateDiscretePie(t_clado,
#                                          cladogenetic = T,
#                                          tip_pie_diameter = 2,
#                                          node_pie_diameter = 2,
#                                          shoulder_pie_diameter = 2,
#                                          state_labels = state_info$state,
#                                          state_colors = c(as.vector(state_info$color), "grey20"),
#                                          pie_nudge_x = 0.01,
#                                          pie_nudge_y= 0.1,
#                                          show_tip_labels = F)
# ggsave("PieRange_plot.pdf", PieRange_plot, width = 7, height = 7)
#
#
# t_nonclado <- treeio::read.beast("ancestral_states_pie.tree")
#
# PieState_plot <- plotAncStateDiscretePie(t_nonclado,
#                                          cladogenetic = F,
#                                          state_labels = c(0:2),
#                                          state_colors = c("red", "yellow", "purple", "grey"),
#                                          tip_pie_diameter = 2,
#                                          node_pie_diameter = 2,
#                                          pie_nudge_x = 0.01,
#                                          pie_nudge_y= 0.1,
#                                          show_tip_labels = F)
# ggsave("PieState_plot.pdf", PieState_plot, width = 7, height = 7)



## TO DO
# collapsing states in the legend no longer seems to work in this function (this needs to be done using processAncStateDiscrete)

#require("rvcheck")
#require("ggimage")
#require("tibble")
#require("ggplot2")
#require("ggtree")
#require("tibble")

