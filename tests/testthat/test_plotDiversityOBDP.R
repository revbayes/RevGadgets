context("tests the plotDiversityOBDP function")
#note: does not compare the generated plot to the expectation
test_that("plot diversity OBDP", {
    file_plot_orig <-
          system.file("extdata", "graphs/plotDiversityOBDP_df.rds", package = "RevGadgets")
    start_time_trace_file <- 
          system.file("extdata", "obdp/start_time_trace.p", package="RevGadgets")
    popSize_distribution_matrices_file <- 
          system.file("extdata", "obdp/Kt_trace.p", package="RevGadgets")
    trees_trace_file <- 
          system.file("extdata", "obdp/mcmc_OBDP_trees.p", package="RevGadgets")
      
    capture.output(Kt_mean <- readOBDP( start_time_trace_file=start_time_trace_file, 
                        popSize_distribution_matrices_file=popSize_distribution_matrices_file, 
                        trees_trace_file=trees_trace_file ))
    
    plot_new <- plotDiversityOBDP( Kt_mean )
    
    plot_orig <- readRDS(file_plot_orig)
    
    tmp <- tempdir()
    pdf(paste0(tmp,"/Rplots.pdf"))
    
    # check that plot doesn't error out
    expect_error(print(plot_new), NA)
    
    dev.off()
    
    # compare plot data objects
    expect_equal(plot_new$data, plot_orig)
})
