context("tests the plotMassExtinctions function")

test_that("tests mass extinctions example", {
  # get files
  mass_extinction_probability_file <-
    system.file("extdata",
                "mass_extinction/crocs_mass_extinction_probabilities_mini.p",
                package = "RevGadgets")

  plot_file <-
    system.file("extdata",
                "graphs/plotMassExtinctions_df.rds",
                package = "RevGadgets")

  mass_extinction_probabilities <-
    readTrace(mass_extinction_probability_file, burnin = 0)

  # prior probability of mass extinction at any time
  prior_n_expected <- 0.1
  n_intervals <- 100
  prior_prob <- prior_n_expected / (n_intervals - 1)

  # times when mass extinctions were allowed
  tree_age <- 243.5
  interval_times <-
    tree_age * seq(1 / n_intervals,
                   (n_intervals - 1) / n_intervals,
                   1 / n_intervals)

  # then plot results:
  plot_new <-
    plotMassExtinctions(
      mass_extinction_trace = mass_extinction_probabilities,
      mass_extinction_times = interval_times,
      mass_extinction_name = "mass_extinction_probabilities",
      prior_prob
    )

  # read original plot object
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))
  # test for errors in plot_new
  expect_error(print(plot_new), NA)
  dev.off()

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

})
