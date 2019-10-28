context("tests the plotAncStatesDiscrete function")

test_that("plots discrete ancestral states", {

  # get file names
  ase_tree_file    <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  state_label_file <- system.file("extdata", "dec/range_colors.n4.txt", package="RevGadgets")
  original_file    <- system.file("extdata", "dec/simple.plot.rds", package="RevGadgets")

  # read in the label/color file
  state_info <- read.csv(state_label_file, header=T, stringsAsFactors=F, sep=",")

  # read in the tree file
  tree_obj <- processAncStatesDiscrete(t=ase_tree_file, state_labels=state_info$range)

  # create the plotting object
  plot_obj <- plotAncStatesDiscrete(tree=tree_obj, state_colors=state_info$color)

  # read original process object
  original_obj <- readRDS(original_file)

  # strict match of test and original objects
  expect_equal(identical(plot_obj, original_obj), TRUE)

})
