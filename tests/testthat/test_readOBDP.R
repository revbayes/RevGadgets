context("tests the readOBDP function")
test_that("read and format OBDP outputs", {
    start_time_trace_file <- 
          system.file("extdata", "obdp/start_time_trace.p", package="RevGadgets")
    popSize_distribution_matrices_file <- 
          system.file("extdata", "obdp/Kt_trace.p", package="RevGadgets")
    trees_trace_file <- 
          system.file("extdata", "obdp/mcmc_OBDP_trees.p", package="RevGadgets")
      
    capture.output(Kt_mean <- readOBDP( start_time_trace_file=start_time_trace_file, 
                        popSize_distribution_matrices_file=popSize_distribution_matrices_file, 
                        trees_trace_file=trees_trace_file ))
    
    expect_equal(length(Kt_mean), 96)
    expect_equal(class(Kt_mean), "data.frame")
    expect_equal(nrow(Kt_mean), 100)
    
})
