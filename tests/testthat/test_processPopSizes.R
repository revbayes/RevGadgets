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
                burnin = 0,
                spacing = "equal",
                min_age = 0,
                max_age = 300000,
                distribution = TRUE
        )
        
        expect_equal(ncol(constant), 4)
        expect_equal(class(constant), c("tbl_df", "tbl", "data.frame"))
        expect_equal(nrow(constant), 100)
        expect_equal(max(constant$time), 1e5)
        expect_equal(constant$value[1], 313240.8)
        
        expect_equal(ncol(CPP), 4)
        expect_equal(class(CPP), c("tbl_df", "tbl", "data.frame"))
        expect_equal(nrow(CPP), 200)
        expect_equal(max(CPP$time), 466634)
        expect_equal(round(CPP$value[c(1,200)]), c(505782, 27198800))
        
        expect_equal(ncol(GMRF), 11)
        expect_equal(class(GMRF), c("matrix", "array"))
        expect_equal(nrow(GMRF), 100)
        expect_equal(GMRF[1,], c(545066, 465357, 603280, 394040, 682751, 535117, 519354, 535540, 502588, 729053, 737713))
        expect_equal(GMRF[100,], c(601042, 375759, 137952, 344189, 312312, 547615, 405185, 347487, 498597, 126805, 635706))
})
