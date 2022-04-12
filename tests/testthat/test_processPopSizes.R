context("tests the processPopSizes function")

test_that("processes population size output", {
        file_popsizes_constant <-
                system.file("extdata",
                            "pop_size/horses_constant_popsizes_mini.p",
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

        constant <- processPopSizes(
                population_size_log = file_popsizes_constant,
        )
        
        CPP <- processPopSizes(
                population_size_log = file_popsizes_CPP,
                interval_change_points_log = file_changepoints_CPP,
                num_grid_points = 200
        )
        
        GMRF <- processPopSizes(
                population_size_log = file_popsizes_GMRF,
                interval_change_points_log = file_changepoints_GMRF,
                max_age = 10000
        )
        
        expect_equal(ncol(constant), 4)
        expect_equal(class(constant), c("tbl_df", "tbl", "data.frame"))
        expect_equal(nrow(constant), 100)
        expect_equal(max(constant$time), 1e5)
        
        expect_equal(ncol(CPP), 4)
        expect_equal(class(CPP), c("tbl_df", "tbl", "data.frame"))
        expect_equal(nrow(CPP), 200)
        expect_equal(max(CPP$time), 466634)
        
        expect_equal(ncol(GMRF), 4)
        expect_equal(class(GMRF), c("tbl_df", "tbl", "data.frame"))
        expect_equal(nrow(GMRF), 100)
        expect_equal(max(GMRF$time), 10000)
})
