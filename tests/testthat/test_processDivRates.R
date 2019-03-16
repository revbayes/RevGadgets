context("tests the processDivRates function")

# relies heavily on readTrace(), so we only test for elements not
# included in the readTrace() testing.

test_that("processes birth-death scripts", {
  primates <- processDivRates(speciation_time_log = "../../data/epi_bd/primates_EBD_speciation_times.log",
                              speciation_rate_log = "../../data/epi_bd/primates_EBD_speciation_rates.log",
                              extinction_time_log = "../../data/epi_bd/primates_EBD_extinction_times.log",
                              extinction_rate_log = "../../data/epi_bd/primates_EBD_extinction_rates.log",
                              burnin = 0.25)
  expect_equal(length(primates), 6)
  expect_equal(class(primates), "list")
  expect_equal(nrow(primates[[1]]), 3750)
})
