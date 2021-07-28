context("Tests Bayes Factor calculations from MRF models")

test_that("compare expected calculations from documentation example", {
  # read in and process primates diversification rate data
  speciation_time_file <-
    system.file("extdata",
                "epi_bd/primates_EBD_speciation_times_mini.p",
                package = "RevGadgets")
  speciation_rate_file <-
    system.file("extdata",
                "epi_bd/primates_EBD_speciation_rates_mini.p",
                package = "RevGadgets")

  speciation_rate <- readTrace(speciation_rate_file, burnin = 0)
  speciation_times <- readTrace(speciation_time_file, burnin = 0)

  expect_equal(
    2.773,
    round(calculateShiftBayesFactor(speciation_rate,
                                    speciation_times,
                                    "speciation",
                                    "interval_times",
                                    0.0,
                                    40.0,
                                    decrease = FALSE),
          digits = 3)
  )

  expect_error(
    calculateShiftBayesFactor(
      speciation_rate,
      speciation_times,
      "fossilization",
      "interval_times",
      0.0,
      40.0,
      decrease = FALSE
    )
  )

})
