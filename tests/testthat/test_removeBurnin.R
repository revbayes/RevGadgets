context("tests the removeBurnin function")

test_that("removes burnin", {
  # load in the trace file
  file <- system.file("extdata",
                      "sub_models/primates_cytb_GTR_mini.p",
                      package = "RevGadgets")

  one_trace <- readTrace(paths = file, burnin = 0)

  one_trace_burnin <- removeBurnin(trace = one_trace, burnin = 0.1)

  #check the new length of the trace
  expect_equal(dim(one_trace_burnin[[1]])[1], 9)

})
