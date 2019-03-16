context("tests the plotDivRates function")

test_that("plots speciation rates", {
  primates <- processDivRates(speciation_time_log = "../../data/epi_bd/primates_EBD_speciation_times.log",
                              speciation_rate_log = "../../data/epi_bd/primates_EBD_speciation_rates.log",
                              extinction_time_log = "../../data/epi_bd/primates_EBD_extinction_times.log",
                              extinction_rate_log = "../../data/epi_bd/primates_EBD_extinction_rates.log",
                              burnin = 0.25)

  # check that this throws an error
  t_bad <- try(plotDivRates(primates, fig_types = c("speciation rate"), col = "green"))
  expect_equal(class(t_bad), "try-error")

  t_good <- try(plotDivRates(primates, fig_types = c("speciation rate")))
  expect_null(t_good)

})

test_that("plots extinction rates", {
  primates <- processDivRates(speciation_time_log = "../../data/epi_bd/primates_EBD_speciation_times.log",
                              speciation_rate_log = "../../data/epi_bd/primates_EBD_speciation_rates.log",
                              extinction_time_log = "../../data/epi_bd/primates_EBD_extinction_times.log",
                              extinction_rate_log = "../../data/epi_bd/primates_EBD_extinction_rates.log",
                              burnin = 0.25)

  # check that this throws an error
  t_bad <- try(plotDivRates(primates, fig_types = c("extinction rate"), col = "green"))
  expect_equal(class(t_bad), "try-error")

  t_good <- try(plotDivRates(primates, fig_types = c("extinction rate")))
  expect_null(t_good)

})

test_that("plots net-diversification rates", {
  primates <- processDivRates(speciation_time_log = "../../data/epi_bd/primates_EBD_speciation_times.log",
                              speciation_rate_log = "../../data/epi_bd/primates_EBD_speciation_rates.log",
                              extinction_time_log = "../../data/epi_bd/primates_EBD_extinction_times.log",
                              extinction_rate_log = "../../data/epi_bd/primates_EBD_extinction_rates.log",
                              burnin = 0.25)

  # check that this throws an error
  t_bad <- try(plotDivRates(primates, fig_types = c("net-diversification rate"), col = "green"))
  expect_equal(class(t_bad), "try-error")

  t_good <- try(plotDivRates(primates, fig_types = c("net-diversification rate")))
  expect_null(t_good)

})

test_that("plots relative-extinction rates", {
  primates <- processDivRates(speciation_time_log = "../../data/epi_bd/primates_EBD_speciation_times.log",
                              speciation_rate_log = "../../data/epi_bd/primates_EBD_speciation_rates.log",
                              extinction_time_log = "../../data/epi_bd/primates_EBD_extinction_times.log",
                              extinction_rate_log = "../../data/epi_bd/primates_EBD_extinction_rates.log",
                              burnin = 0.25)

  # check that this throws an error
  t_bad <- try(plotDivRates(primates, fig_types = c("relative-extinction rate"), col = "green"))
  expect_equal(class(t_bad), "try-error")

  t_good <- try(plotDivRates(primates, fig_types = c("relative-extinction rate")))
  expect_null(t_good)

})
