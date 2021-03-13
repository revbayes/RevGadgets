context("tests the readTrace function")

test_that("reads one trace", {
  file <- system.file("extdata", "sub_models/primates_cytb_GTR.p", package="RevGadgets")
  trace_single <- readTrace(path = file, burnin = 0)
  expect_equal(length(trace_single), 1)
  expect_equal(class(trace_single), "list")
  expect_equal(class(trace_single[[1]]), "data.frame")
  expect_equal(nrow(trace_single[[1]]), 1001)
  expect_equal(ncol(trace_single[[1]]), 59)
})

test_that("reads multiple traces", {
  file_1 <- system.file("extdata", "sub_models/primates_cytb_GTR_run_1.p", package="RevGadgets")
  file_2 <- system.file("extdata", "sub_models/primates_cytb_GTR_run_2.p", package="RevGadgets")
  trace_multi <- readTrace(path = c(file_1, file_2), burnin = 0 )
  expect_equal(length(trace_multi), 2)
  expect_equal(class(trace_multi), "list")
  expect_equal(class(trace_multi[[1]]), "data.frame")
  expect_equal(nrow(trace_multi[[1]]), 501)
  expect_equal(ncol(trace_multi[[1]]), 58)
})
