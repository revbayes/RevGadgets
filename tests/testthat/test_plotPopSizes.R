context("tests the plotPopSizes function")

test_that("plot population size trajectories", {
        file_plot_orig <-
                system.file("extdata",
                            "graphs/plotPopSizes.rds",
                            package = "RevGadgets")
        
        file_popsizes_CPP <-
                system.file("extdata",
                            "pop_size/horses_CPP_popsizes_mini.p",
                            package = "RevGadgets")
        file_changepoints_CPP <-
                system.file("extdata",
                            "pop_size/horses_CPP_times_mini.p",
                            package = "RevGadgets")
        
        file_popsizes_GMRF <-
                system.file("extdata",
                            "pop_size/horses_GMRF_popsizes_mini.p",
                            package = "RevGadgets")
        file_changepoints_GMRF <-
                system.file("extdata",
                            "pop_size/horses_GMRF_times_mini.p",
                            package = "RevGadgets")
        
        CPP <- processPopSizes(
                population_size_log = file_popsizes_CPP,
                interval_change_points_log = file_changepoints_CPP,
                num_grid_points = 200
        )
        
        GMRF <- processPopSizes(
                population_size_log = file_popsizes_GMRF,
                interval_change_points_log = file_changepoints_GMRF,
                spacing = "equal",
                min_age = 0,
                max_age = 300000
        )
        
        plot_new_1 <- plotPopSizes(CPP)
        plot_new_2 <- plotPopSizes(GMRF, add = TRUE, 
                                   existing_plot = plot_new_1, 
                                   col = "blue")
        
        plot_orig <- readRDS(file_plot_orig)

        tmp <- tempdir()
        pdf(paste0(tmp,"/Rplots.pdf"))

        # check that plot doesn't error out
        expect_error(print(plot_new_2), NA)

        dev.off()

        # compare plot data objects
        expect_equal(plot_new_2$data, plot_orig$data)
        # expect_equal(plot_new_2$scales, plot_orig$scales)
        # expect_equal(plot_new_2$layers, plot_orig$layers)
        
})
