context("tests the SBC processing function")

test_that("processSBC produces correct output", {
  # path of the test directory
  sbc_test_path <- system.file("extdata", "output_sbc_test", package = "RevGadgets")

  # run the processing function on the test data
  sbc_results <- processSBC(path = sbc_test_path, n_reps = 5, verbose = FALSE)

  # find the path to our pre-computed "correct" results
  expected_file <- system.file("extdata", "sbc_results_expected.rds", package = "RevGadgets")
  expected_results <- readRDS(expected_file)

  # the output of the function should be identical to the saved, correct output
  expect_equal(sbc_results, expected_results)

  # check the structure of the output
  expect_true(is.list(sbc_results))
  expect_equal(sort(names(sbc_results)), c("alpha", "sigma"))
})
