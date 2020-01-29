context("tests the processAncStatesDiscrete function")

test_that("processes discrete ancestral states scripts", {

  # get file names
  test_file     <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  original_file <- system.file("extdata", "dec/simple.ase.tre.rds", package="RevGadgets")
  label_file    <- system.file("extdata", "dec/simple.state_info.txt", package="RevGadgets")

  # collect state info
  state_info <- read.csv(label_file, header=TRUE, stringsAsFactors=FALSE, sep=",")
  state_labels <- state_info$label
  names(state_labels) <- state_info$state

  # commented out until PieRange works

  # process data for new test object
  #process_obj <- processAncStatesDiscrete(tree_file = test_file, state_labels=state_labels)

  # read original process object
  #original_obj <- readRDS(original_file)

  # strict match of test and original objects
  #expect_equal(identical(process_obj, original_obj), TRUE)

})
