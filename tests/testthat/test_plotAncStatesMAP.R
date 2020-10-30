context("tests the plotAncStatesMAP function")

test_that("plots MAP of ancestral states", {

  # get files
  tree_file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
  plot_file <- system.file("extdata", "graphs/plotAncStatesMAP.RData", package="RevGadgets")

  # make a new plot
  example <- processAncStates(tree_file, state_labels = c("1" = "Awesome", "2" = "Beautiful", "3" = "Cool!"))
  plot_new <- plotAncStatesMAP(t = example)

  # read original plot object
  load(plot_file)

  # compare plot objects
  expect_equal(plot[[1]], plot_new[[1]])

})
