# run "stoch_map_test_tmp/simmaps.Rev" to generate output for this test

# read a tree
tree <- readTrees("stoch_map_test_tmp/tree.nexus")[[1]][[1]]

# process samples
stoch_map_df <- RevGadgets:::processStochMaps(tree, "stoch_map_test_tmp/maps.log", states = c("0","1"), burnin = 0.1)
maps <- stoch_map_df
#plot samples:

colors <- c("0" = "red", "1" = "blue")
plotStochMaps <- function(tree,
                          maps,
                          colors,
                          layout) {
  states <- names(colors)
  
  p <- plotTree(tree = list(list(tree)),
                lineend = "square")
  
  dat <- dplyr::left_join(maps, p$data, by = "node")

  max <- apply(dat[,states], MARGIN = 1, which.max) 
  names(max) <- NULL
  dat$map_state <- states [unlist(max) ]

  # horizontal segments
  
  seg_horiz <- data.frame(x = dat$x - dat$x0,
                          xend = dat$x -dat$x1,
                          y = dat$y,
                          yend = dat$y,
                          col = dat$map_state)
  
  seg_vert <- data.frame(x = dat$xend,
                         xend = dat$x -dat$x1,
                         y = dat$y,
                         yend = dat$y,
                         col = dat$map_state)
  
  p + ggplot2::geom_segment(data = seg_horiz, ggplot2::aes(x = x, y = y,
                                                    xend = xend, yend = yend,
                                                    color = col),
                            lineend = "square") +
    ggplot2::scale_color_manual(values = colors)
  
  
  
}