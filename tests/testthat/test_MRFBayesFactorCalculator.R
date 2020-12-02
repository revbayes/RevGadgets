context("Tests the MRFBayesFactorCalculator")

test_that("compare expected calculations from documentation example", {
  # read in and process primates diversification rate data
  speciation_time_file <- system.file("extdata",
                                      "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
  speciation_rate_file <- system.file("extdata",
                                      "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
  extinction_time_file <- system.file("extdata",
                                      "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
  extinction_rate_file <- system.file("extdata",
                                      "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

  primates <- processDivRates(speciation_time_log = speciation_time_file,
                              speciation_rate_log = speciation_rate_file,
                              extinction_time_log = extinction_time_file,
                              extinction_rate_log = extinction_rate_file,
                              burnin = 0.25)
  expect_equal(5.1857, calculateShiftBayesFactor(primates,"speciation rate","speciation time",0.0,40.0,decrease=FALSE))
})
