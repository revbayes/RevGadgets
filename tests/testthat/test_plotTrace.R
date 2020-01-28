context("tests the plotTrace function")

test_that("plot pi traces", {
  # load in the trace file
  file_1 <- system.file("extdata",
                      "sub_models/primates_cytb_covariotide.p",
                      package="RevGadgets")
  one_trace <- readTrace(path = file_1)
  # produce the plot pi parameters object
  plots_new <- plotTrace(trace = one_trace,
                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))
  # load the saved plot for comparison
  file_2 <- system.file("extdata",
                        "graphs/plotTrace_pi.Rdata",
                        package="RevGadgets")
  load(file_2) # loads an object called 'plots'
  expect_equal(plots[[1]], plots_new[[1]])
  })
