context("tests the Validation results plotting function")

test_that("plotValidation generates valid ggplot objects", {
  # load the results to use as input
  expected_file <- system.file("extdata", "validation_results_expected.rds", package = "RevGadgets")
  validation_results <- readRDS(expected_file)

  # generate a plot for a single parameter
  p_alpha <- plotValidation(validation_results, parameter = "alpha")

  # generate a list of plots for all parameters
  p_list <- plotValidation(validation_results)

  # the function should return a ggplot object
  expect_s3_class(p_alpha, "ggplot")

  # the lists of plots should contain ggplot objects
  expect_true(is.list(p_list))
  expect_s3_class(p_list$sigma, "ggplot")

  # the plots should be able to be printed without error
  pdf(NULL)
  expect_error(print(p_alpha), NA)
  expect_error(print(p_list$alpha), NA)
  dev.off()
})
