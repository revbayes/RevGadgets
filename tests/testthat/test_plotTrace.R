context("tests the plotTrace function")

test_that("plot pi traces", {
  # load in the trace file
  file_1 <- system.file("extdata",
                      "sub_models/primates_cytb_covariotide.p",
                      package="RevGadgets")
  one_trace <- readTrace(path = file_1)
  # produce the plot pi parameters object
  plot_new <- plotTrace(trace = one_trace,
                     vars = c("pi[1]","pi[2]","pi[3]","pi[4]"))[[1]]
  # load the saved plot for comparison
  plot_file <- system.file("extdata",
                        "graphs/plotTrace_pi.rds",
                        package="RevGadgets")
  #read in original plot
  plot_orig <- readRDS(plot_file)[[1]]

  # compare plot objects
  for (i in 1:length(plot_new)) {
    if (names(plot_new[i]) != "plot_env") {
      expect_equal(plot_new[[i]], plot_orig[[i]])
    }
  }

})
