context("tests the processAncStatesDiscrete function")

# relies heavily on readTrace(), so we only test for elements not
# included in the readTrace() testing.

test_that("processes discrete ancestral states scripts", {
  test_file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  original_file <- system.file("extdata", "dec/simple.ase.tre.rds", package="RevGadgets")

  # process data for new test object
  test_obj <- processAncStatesDiscrete(tree_file = test_file)

  # read original process object
  original_obj <- readRDS(original_file)

  # strict match of test and original objects
  expect_equal(identical(test_obj, original_obj), TRUE)

})
