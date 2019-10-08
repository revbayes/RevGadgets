context("tests the plotTrace function")

test_that("reads one trace", {
  file <- system.file("extdata",
                      "sub_models/primates_cytb_covariotide.p",
                      package="RevGadgets")
  one_trace <- readTrace(path = file, burnin = 0)
  plots <- plotTrace(trace = one_trace,
                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
  load("inst/extdata/graphs/plotTrace_pi.RData")
  expect_equal()
