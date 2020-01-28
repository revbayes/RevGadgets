#' Plot tree
#'
#' Plots a single tree, such as an MCC or MAP tree.
#'
#' Plots a single tree, such as an MCC or MAP tree, with
#' optionally labele posterior probabilities at nodes, a
#' timescale plotted on the x - axis, and 95\% CI for node ages.
#'
#'
#' @param tree (list of lists of tidytree objects; no default) Name of a list of lists of
#' tidytree objects, such as produced by readTrees(). This object should only contain
#' only one summary tree from one trace file. If it contains multiple trees or multiple
#' traces, only the first will be used.
#'
#' @param chrono (logical; FALSE) Plot as a chronogram with node age bars and time-labeled
#' x-axis?
#'
#' @param node_labels (logical; TRUE) Plot posterior probabilities as colored cirlces at
#' nodes?
#'
#' @param tip_labels (logical: TRUE) Plot tip labels?
#'
#' @param fossils (character vector; NULL) Optional vector of the tips that are fossils,
#' to be sure that these tips are highlighted in red.
#'
#' @return returns a single plot object.
#'
#' @examples
#'
#' @export

plotTree <- function(tree, chrono = FALSE, node_labels = TRUE,
                     tip_labels = TRUE, fossils = NULL) {
  recover()
  # grab single tree from input
  phy <- tree[[1]][[1]]

  # format posterior data
  phy@data$posterior[ phy@data$posterior == 1 ] <- NA

  # Encountered problems with using geom_range to plot age HPDs in ggtree. It
  # appears that geom_range incorrectly rotates the HPD relative to the height
  # of the node unnecessarily. My guess for this would be because older version
  # of ggtree primarily supported length measurements, and not height measurements
  # so the new capability to handle height might contain a "reflection" bug.
  #
  # For example, suppose a node has height 3 with HPD [2, 7]. You can think of
  # this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
  # causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
  #
  # See this excellent trick by Tauana:
  # https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
  # Adapted this code to also plot fossil tip uncertainty in red

  phy@data$age_0.95_HPD <- lapply(phy@data$age_0.95_HPD, function(z) {
    if (is.na(z)) { return(c(NA,NA)) } else { return(z) }
  })

  minmax <- t(matrix(unlist(phy@data$age_0.95_HPD), nrow=2))
  bar_df <- data.frame(node_id=as.integer(phy@data$node), as.data.frame(minmax))
  names(bar_df) <- c("node_id", "min", "max")
  fossil_df <-  dplyr::filter(bar_df, node_id %in% match(fossils,phy@phylo$tip.label))
  node_df <- dplyr::filter(bar_df, node_id > ape::Ntip(phy@phylo))

  # get dimensions
  n_nodes <- treeio::Nnode(phy)
  max_age <- 90
  dx <- max_age %% 10

  # plot
  pp <- ggtree::ggtree(phy, right=F )

  # plot age densities
  node_df <- dplyr::left_join(node_df, pp$data, by=c("node_id"="node"))
  node_df <- dplyr::select(node_df,  node_id, min, max, y)
  fossil_df <- dplyr::left_join(fossil_df, pp$data, by=c("node_id"="node"))
  fossil_df <- dplyr::select(fossil_df, node_id, min, max, y)
  pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y),
                                   data=node_df, color="blue", size=1.5, alpha=0.3)
  pp <- pp + ggplot2::geom_segment(ggplot2::aes(x=-min, y=y, xend=-max, yend=y),
                                   data=fossil_df, color="red", size=1.5, alpha=0.4)

  # set coordinates
  pp <- pp + ggplot2::coord_cartesian(xlim = c(-max_age,30), ylim=c(-7, n_nodes+1.5), expand=F)
  pp <- pp + ggplot2::scale_x_continuous(breaks = seq(-max_age-dx,0,10), labels = rev(seq(0,max_age+dx,10)))
  pp <- pp + ggtree::theme_tree2()
  pp <- pp + ggplot2::labs(x="Age (Ma)")
  pp <- pp + ggplot2::theme(legend.position=c(.05, .955), axis.line = ggplot2::element_line(colour = "black"))
  pp <- ggtree::revts(pp)
  pp <- .add_epoch_times(pp, max_age, dy_bars=-7, dy_text=-3)

  # plot clade support
  pp$data$posterior_class = NA
  pp$data$posterior_class[ which(pp$data$posterior>=0.99) ] = ">0.99"
  pp$data$posterior_class[ which(pp$data$posterior<0.99&pp$data$posterior>0.95) ] = ">0.95"
  pp$data$posterior_class[ which(pp$data$posterior<=0.95&pp$data$posterior>0.75) ] = ">0.75"
  pp$data$posterior_class[ which(pp$data$posterior<=0.75&pp$data$posterior>0.5) ] = ">0.50"
  pp$data$posterior_class[ which(pp$data$posterior<=0.5) ] = "<0.50"
  pp$data$posterior_class = factor(pp$data$posterior_class, levels=c(">0.99",">0.95",">0.75", ">0.50", "<0.50"))
  pp$data$label = sapply(pp$data$label, function(x) { gsub("_"," ",x) })
  pp$data$label = sapply(pp$data$label, function(x) { gsub("subsp. ","",x) })

  # set posterior support
  pp <- pp + ggtree::geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], size=2, color="black")
  pp <- pp + ggtree::geom_nodepoint( data=pp$data[ !is.na(pp$data$posterior), ], ggplot2::aes(color=posterior_class), size=1.5)
  col_post <- c("#000000","#666666","#999999","#BBBBBB","#EEEEEE")
  names(col_post) <- levels(pp$data$posterior_class)
  pp <- pp + ggplot2::scale_color_manual(values=col_post, name="Posterior")

  # set tip & clade names
  pp <- pp + ggtree::geom_tiplab(size=2.5, offset=0.2)

  return(pp)
}
