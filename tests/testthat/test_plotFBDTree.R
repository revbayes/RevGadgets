context("tests the plotFBDTree function")

test_that("plots FBD tree", {
  # get files
  tree_file <-
    system.file("extdata", "fbd/bears.mcc.tre", package = "RevGadgets")
  plot_file <-
    system.file("extdata", "graphs/plotFBDTree_df.rds", package = "RevGadgets")

  # make a new plot
  example <- readTrees(paths = tree_file)
  plot_new <-
    plotFBDTree(
      tree = example,
      timeline = T,
      tip_labels_italics = F,
      tip_labels_remove_underscore = T,
      tip_age_bars = T,
      node_age_bars = T,
      age_bars_colored_by = "posterior",
      age_bars_color = rev(colFun(2))
    ) + ggplot2::theme(legend.position = c(.25, .85))

  # read original plot data object
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))
  # test for errors in plot_new
  expect_error(print(plot_new), NA)
  dev.off()

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

})
