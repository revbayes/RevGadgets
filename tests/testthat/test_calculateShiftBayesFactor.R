context("Tests Bayes Factor calculations from MRF models")

test_that("compare expected calculations from documentation example", {
  # read in and process primates diversification rate data
  speciation_time_file <- system.file("extdata",
                                      "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
  speciation_rate_file <- system.file("extdata",
                                         "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")

  speciation_rate <- readTrace(speciation_rate_file,burnin = 0.25)
  speciation_times <- readTrace(speciation_time_file,burnin = 0.25)

  expect_equal(5.1857, calculateShiftBayesFactor(speciation_rate,speciation_times,"speciation","interval_times",0.0,40.0,decrease=FALSE))
})
