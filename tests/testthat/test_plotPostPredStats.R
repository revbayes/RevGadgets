context("Tests plotPostPredStats function")

test_that("compare processed output from documentation example", {
  file_sim <- system.file("extdata",
                          "PPS/simulated_data_pps_mini.csv", package =
                            "RevGadgets")
  file_emp <- system.file("extdata",
                          "PPS/empirical_data_pps_example.csv", package =
                            "RevGadgets")
  file_old_plot <- system.file("extdata",
                               "graphs/plotPostPredStats_df.rds", package =
                                 "RevGadgets")
  t <- processPostPredStats(path_sim = file_sim,
                            path_emp = file_emp)

  plots <- plotPostPredStats(data = t)
  plot_new <- plots[[1]]
  plot_orig <- readRDS(file = file_old_plot)

  expect_equal(length(plots), dim(t[[1]])[2])

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))
  # test for errors in plot_new
  expect_error(print(plot_new), NA)
  dev.off()

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

})
