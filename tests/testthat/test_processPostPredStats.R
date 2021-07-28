context("Tests processPostPredStats function")

test_that("compare processed output from documentation example", {
   file_sim <- system.file("extdata",
                           "PPS/simulated_data_pps_mini.csv",
                           package = "RevGadgets")
   file_emp <- system.file("extdata",
                           "PPS/empirical_data_pps_example.csv",
                           package = "RevGadgets")
   t <- processPostPredStats(path_sim = file_sim,
                             path_emp = file_emp)

   expect_equal(class(t), "list")
   expect_equal(length(t), 2)
   expect_equal(names(t), c("simulated", "observed"))

})
