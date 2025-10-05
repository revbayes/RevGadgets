context("tests the Validation output processing function")

test_that("processValidation produces correct output", {
  # path of the test directory
  validation_test_path <- system.file("extdata", "output_validation_test", package = "RevGadgets")

  # run the processing function on the test data
  validation_results <- processValidation(path = validation_test_path, n_reps = 5, verbose = FALSE)

  # find the path to our pre-computed "correct" results
  expected_file <- system.file("extdata", "validation_results_expected.rds", package = "RevGadgets")
  expected_results <- readRDS(expected_file)

  # the output of the function should be identical to the saved, correct output
  expect_equal(validation_results, expected_results)

  # check the structure of the output
  expect_true(is.list(validation_results))
  expect_equal(sort(names(validation_results)), c("alpha", "sigma"))
})
