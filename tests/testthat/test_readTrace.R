context("tests the readTrace function")

test_that("dummy test", {
  t <- 5
  expect_equal(5, t)
})


# test_that("reads one trace", {
#   trace_single <- readTrace(path = "../../data/sub_models/primates_cytb_covariotide.log")
#   expect_equal(length(trace_single), 1)
#   expect_equal(class(trace_single), "list")
#   expect_equal(class(trace_single[[1]]), "data.frame")
#   expect_equal(nrow(trace_single[[1]]), 4002)
#   expect_equal(ncol(trace_single[[1]]), 17)
# })
#
# test_that("reads multiple traces", {
#   trace_multi <- readTrace(path = c("../../data/sub_models/primates_cytb_covariotide_run_1.log",
#                                     "../../data/sub_models/primates_cytb_covariotide_run_2.log"))
#   expect_equal(length(trace_multi), 2)
#   expect_equal(class(trace_multi), "list")
#   expect_equal(class(trace_multi[[1]]), "data.frame")
#   expect_equal(nrow(trace_multi[[1]]), 2001)
#   expect_equal(ncol(trace_multi[[1]]), 16)
# })
