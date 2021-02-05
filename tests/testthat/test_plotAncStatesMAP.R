context("tests the plotAncStatesMAP function")

test_that("plots MAP of ancestral states", {
  # get files
  tree_file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
  plot_file <- system.file("extdata", "graphs/plotAncStatesMAP.rds", package="RevGadgets")

  # make a new plot
  example <- processAncStates(tree_file, state_labels = c("1" = "Awesome", "2" = "Beautiful", "3" = "Cool!"))
  plot_new <- plotAncStatesMAP(t = example)
  print(plot_new)

  # read original plot object
  plot_orig <- readRDS(plot_file)

  # compare plot objects

  for (i in 1:length(plot_new)){
    if (names(plot_new[i]) != "plot_env"){
      expect_equal(plot_new[[i]], plot_orig[[i]])
    }
  }
})
