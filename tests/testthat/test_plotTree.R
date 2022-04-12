context("tests the plotTree function")

test_that("plot basic, not-yet-rooted phylogeny", {
  # load in the trace file
  file_1 <-  system.file("extdata",
                         "sub_models/primates_cytb_GTR_MAP.tre",
                         package = "RevGadgets")
  tree <- readTrees(paths = file_1)
  # produce the plot pi parameters object
  plot_new <- plotTree(tree = tree, node_labels = "posterior")
  #print(plot_new)
  # load the saved plot for comparison
  file_2 <- system.file("extdata",
                        "graphs/plotTree_df.rds",
                        package = "RevGadgets")
  plot_orig <- readRDS(file_2) # loads an object called 'plot'

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))

  # test for errors in plot_new
  expect_error(print(plot_new), NA)

  #  compare plot data objects
  expect_equal(plot_new$data, plot_orig)

  dev.off()
})
