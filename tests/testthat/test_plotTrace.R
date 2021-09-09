context("tests the plotTrace function")

test_that("plot pi traces", {
  # load in the trace file
  file_1 <- system.file("extdata",
                        "sub_models/primates_cytb_GTR_mini.p",
                        package = "RevGadgets")
  one_trace <- readTrace(path = file_1, burnin = 0)
  # produce the plot pi parameters object
  plot_new <- plotTrace(trace = one_trace,
                        vars = c("pi[1]", "pi[2]", "pi[3]", "pi[4]"))[[1]]
  # load the saved plot for comparison
  plot_file <- system.file("extdata",
                           "graphs/plotTrace_df.rds",
                           package = "RevGadgets")
  #read in original plot
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))
  # test for errors in plot_new
  expect_error(print(plot_new), NA)
  dev.off()

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

})
