context("tests the processAncStatesDiscrete function")

test_that("processes discrete ancestral states scripts", {

  # get file names
  test_file     <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  original_file <- system.file("extdata", "dec/simple.ase.tre.rds", package="RevGadgets")

  # process data for new test object
  process_obj <- processAncStatesDiscrete(tree_file = test_file)

  # read original process object
  original_obj <- readRDS(original_file)

  # strict match of test and original objects
  expect_equal(identical(process_obj, original_obj), TRUE)

})
