context("tests the plotHiSSE function")

test_that("plot HiSSE", {

  # read in and process file
  hisse_file <- system.file("extdata", "sse/primates_HiSSE_2.p", package="RevGadgets")
  pdata <- processSSE(hisse_file)
  plot_new <- plotHiSSE(pdata)

  plot_file <- system.file("extdata", "graphs/plotHiSSE.rds", package="RevGadgets")
  plot_orig <- readRDS(plot_file)

  # test for errors in plot_new
  expect_error(print(plot_new), NA)

  #  compare plot data objects
  expect_equal(plot_new$data, plot_orig$data)

})
