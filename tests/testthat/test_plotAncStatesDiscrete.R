context("tests the plotAncStatesDiscrete function")

test_that("plots discrete ancestral states", {

  # get file names
  test_file     <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  original_file <- system.file("extdata", "dec/simple.plot.rds", package="RevGadgets")
  label_file    <- system.file("extdata", "dec/simple.state_info.txt", package="RevGadgets")

  # collect state info
  state_info <- read.csv(label_file, header=TRUE, stringsAsFactors=FALSE, sep=",")
  state_labels <- state_info$label; names(state_labels) <- state_info$state
  state_colors <- state_info$color; names(state_colors) <- state_labels

  # process data for new test object
  process_obj <- processAncStatesDiscrete(tree_file = test_file, state_labels=state_labels)

  # create the plotting object
  plot_obj <- plotAncStatesDiscrete(t=process_obj, state_colors=state_colors, summary_statistic="PieRange", include_start_states=TRUE)

  # read original process object (temporarily disabled)
  # original_obj <- readRDS(original_file)

  # strict match of test and original objects (temporarily disabled)
  # expect_equal(identical(plot_obj, original_obj), TRUE)

})
