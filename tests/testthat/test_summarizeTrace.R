context("tests the summarize trace function")

test_that("summarizes continuous correctly", {
  file <- system.file("extdata",
                      "sub_models/primates_cytb_GTR_mini.p",
                      package = "RevGadgets")
  one_trace <- readTrace(paths = file)
  trace_sum <- summarizeTrace(trace = one_trace,
                              vars = c("pi[1]", "pi[2]", "pi[3]", "pi[4]"))

  expect_equal(length(trace_sum), 4)
  expect_equal(class(trace_sum), "list")
  expect_equal(names(trace_sum), c("pi[1]", "pi[2]", "pi[3]", "pi[4]"))

  expect_equal(length(trace_sum[[1]]), 1)
  expect_equal(class(trace_sum[[1]]), "list")
  expect_equal(names(trace_sum[[1]]), "trace_1")

  expect_equal(class(trace_sum[[1]][[1]]), "numeric")
  expect_equal(names(trace_sum[[1]][[1]]),
               c("mean", "median", "MAP", "quantile_2.5", "quantile_97.5"))
  expect_equal(round(trace_sum[[1]][[1]]["mean"], digits = 3),
               c("mean" = 0.255))
  expect_equal(round(trace_sum[[1]][[1]]["quantile_2.5"], digits = 3),
               c("quantile_2.5" = 0.236))
  expect_equal(round(trace_sum[[1]][[1]]["quantile_97.5"], digits = 3),
               c("quantile_97.5" = 0.276))
})
