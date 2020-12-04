context("tests the plotDivRates function")

test_that("plots speciation rates", {
 file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
 file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
 file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
 file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

 primates <- processDivRates(speciation_time_log = file_spectimes,
                             speciation_rate_log = file_specrates,
                             extinction_time_log = file_exttimes,
                             extinction_rate_log = file_extrates,
                             burnin = 0.25)

 # check that this throws an error
 t_bad <- try(plotDivRates(primates, fig_types = c("speciation rate"), col = "green"))
 expect_equal(class(t_bad), "try-error")

 t_good <- try(plotDivRates(primates, fig_types = c("speciation rate")))
 expect_null(t_good)

})

test_that("plots extinction rates", {
 file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
 file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
 file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
 file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

 primates <- processDivRates(speciation_time_log = file_spectimes,
                             speciation_rate_log = file_specrates,
                             extinction_time_log = file_exttimes,
                             extinction_rate_log = file_extrates,
                             burnin = 0.25)

 # check that this throws an error
 t_bad <- try(plotDivRates(primates, fig_types = c("extinction rate"), col = "green"))
 expect_equal(class(t_bad), "try-error")

 t_good <- try(plotDivRates(primates, fig_types = c("extinction rate")))
 expect_null(t_good)

})

test_that("plots net-diversification rates", {
 file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
 file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
 file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
 file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

 primates <- processDivRates(speciation_time_log = file_spectimes,
                             speciation_rate_log = file_specrates,
                             extinction_time_log = file_exttimes,
                             extinction_rate_log = file_extrates,
                             burnin = 0.25)

 # check that this throws an error
 t_bad <- try(plotDivRates(primates, fig_types = c("net-diversification rate"), col = "green"))
 expect_equal(class(t_bad), "try-error")

 t_good <- try(plotDivRates(primates, fig_types = c("net-diversification rate")))
 expect_null(t_good)

})

test_that("plots relative-extinction rates", {
 file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
 file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
 file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
 file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

 primates <- processDivRates(speciation_time_log = file_spectimes,
                             speciation_rate_log = file_specrates,
                             extinction_time_log = file_exttimes,
                             extinction_rate_log = file_extrates,
                             burnin = 0.25)

 # check that this throws an error
 t_bad <- try(plotDivRates(primates, fig_types = c("relative-extinction rate"), col = "green"))
 expect_equal(class(t_bad), "try-error")

 t_good <- try(plotDivRates(primates, fig_types = c("relative-extinction rate")))
 expect_null(t_good)

})
