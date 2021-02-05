context("tests the plotFBDTree function")

test_that("plots FBD tree", {

  # get files
  tree_file <- system.file("extdata", "fbd/bears.mcc.tre", package="RevGadgets")
  plot_file <- system.file("extdata", "graphs/plotFBDTree.rds", package="RevGadgets")

  # make a new plot
  example <- readTrees(paths = tree_file)
  plot_new <- plotFBDTree(tree = example, timeline = T, tip_labels_italics = F,
                          tip_labels_remove_underscore = T, tip_age_bars = T,
                          node_age_bars = T, node_age_bars_colored_by = "posterior",
                          node_age_bars_color = rev(RevGadgets:::.colFun(2))) + ggplot2::theme(legend.position=c(.25, .85))
  print(plot_new)
  # read original plot object
  plot_orig <- readRDS(plot_file)

  # compare plot objects
  for (i in 1:length(plot_new)) {
    if (names(plot_new[i]) != "plot_env") {
      expect_equal(plot_new[[i]], plot_orig[[i]])
    }
  }

})
