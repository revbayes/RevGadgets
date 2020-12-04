context("tests the plotDivRates function")

test_that("plot works", {
 file_spectimes <- system.file("extdata", "epi_bd/primates_EBD_speciation_times.p", package="RevGadgets")
 file_specrates <- system.file("extdata", "epi_bd/primates_EBD_speciation_rates.p", package="RevGadgets")
 file_exttimes <- system.file("extdata", "epi_bd/primates_EBD_extinction_times.p", package="RevGadgets")
 file_extrates <- system.file("extdata", "epi_bd/primates_EBD_extinction_rates.p", package="RevGadgets")

 primates <- processDivRates(speciation_time_log = file_spectimes,
                             speciation_rate_log = file_specrates,
                             extinction_time_log = file_exttimes,
                             extinction_rate_log = file_extrates,
                             burnin = 0.25)

 expect_silent(t_good <- plotDivRates(primates))

})
