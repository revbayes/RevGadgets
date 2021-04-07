context("tests the plotMuSSE function")

test_that("plot MuSSE", {

  # read in and process file
  bisse_file <- system.file("extdata", "sse/primates_BiSSE_activity_period.p", package="RevGadgets")
  pdata <- processSSE(bisse_file)
  plot_new <- plotMuSSE(pdata)

  plot_file <- system.file("extdata", "graphs/plotMuSSE.rds", package="RevGadgets")
  plot_orig <- readRDS(plot_file)

  # test for errors in plot_new
  expect_error(print(plot_new), NA)

  #  compare plot data objects
  expect_equal(plot_new$data, plot_orig$data)

})
