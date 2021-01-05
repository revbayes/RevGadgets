# context("tests the plotMassExtinctionBayesFactor function")
#
# test_that("plots mass extinctions", {
#
#   # get files
#   mass_extinction_probability_file <- system.file("extdata", "mass_extinction/crocs_mass_extinction_probabilities.p", package="RevGadgets")
#   plot_file <- system.file("extdata", "graphs/plotMassExtinctionBayesFactor.rds", package="RevGadgets")
#
#   # make a new plot
#   mass_extinction_probabilities <- readTrace(mass_extinction_probability_file,burnin = 0.25)
#
#   prior_n_expected <- 0.1
#   n_intervals <- 100
#   prior_prob <- prior_n_expected/(n_intervals-1)
#
#   tree_age <- 243.5
#   interval_times <- tree_age * seq(1/n_intervals,(n_intervals-1)/n_intervals,1/n_intervals)
#
#   plot_new <- plotMassExtinctions(mass.extinction.trace=mass_extinction_probabilities,mass.extinction.times=interval_times,mass.extinction.name="mass_extinction_probabilities",prior_prob)
#
#   # read original plot object
#   plot_orig <- readRDS(plot_file)
#
#   # compare plot objects
#   for (i in 1:length(plot_new)) {
#     if (names(plot_new[i]) != "plot_env") {
#       expect_equal(plot_new[[i]], plot_orig[[i]])
#     }
#   }
#
# })
