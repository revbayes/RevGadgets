context("tests the plotDivRates function")
#note: does not compare the generated plot to the expectation
test_that("plot works", {
        file_plot_orig <-
                system.file("extdata",
                            "graphs/plotDivRates_df.rds",
                            package = "RevGadgets")
        file_spectimes <-
                system.file("extdata",
                            "epi_bd/primates_EBD_speciation_times_mini.p",
                            package = "RevGadgets")
        file_specrates <-
                system.file("extdata",
                            "epi_bd/primates_EBD_speciation_rates_mini.p",
                            package = "RevGadgets")
        file_exttimes <-
                system.file("extdata",
                            "epi_bd/primates_EBD_extinction_times_mini.p",
                            package = "RevGadgets")
        file_extrates <-
                system.file("extdata",
                            "epi_bd/primates_EBD_extinction_rates_mini.p",
                            package = "RevGadgets")

        primates <- processDivRates(
                speciation_time_log = file_spectimes,
                speciation_rate_log = file_specrates,
                extinction_time_log = file_exttimes,
                extinction_rate_log = file_extrates,
                burnin = 0.25,
                summary = "mean"
        )
        plot_new <- plotDivRates(primates, )
        plot_orig <- readRDS(file_plot_orig)

        tmp <- tempdir()
        pdf(paste0(tmp,"/Rplots.pdf"))

        # check that plot doesn't error out
        expect_error(print(plot_new), NA)

        dev.off()

        # compare plot data objects
        expect_equal(plot_new$data, plot_orig)

})
