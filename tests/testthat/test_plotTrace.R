context("tests the plotTrace function")

test_that("reads one trace", {
  file_1 <- system.file("extdata",
                      "sub_models/primates_cytb_covariotide.p",
                      package="RevGadgets")
  one_trace <- readTrace(path = file_1, burnin = 0)
  plots_new <- plotTrace(trace = one_trace,
                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
  file_2 <- system.file("extdata",
                        "graphs/plotTrace_pi.RData",
                        package="RevGadgets")
  load(file_2)
  expect_equal(plots[[1]], plots_new[[1]])
  })
