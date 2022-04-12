context("tests the geomStepRibbon function")

test_that("geomStepRibbon example works", {

  # make new plot
  huron <- data.frame(year = 1875:1972, level = as.vector(LakeHuron))
  plot_new <- ggplot2::ggplot(huron, ggplot2::aes(year)) +
    geom_stepribbon(ggplot2::aes(ymin = level - 1, ymax = level + 1),
                    fill = "grey70") +
      ggplot2::geom_step(ggplot2::aes(y = level))

  # read in old saved version
  plot_file <- system.file("extdata",
                           "graphs/geomStepRibbon_df.rds",
                           package = "RevGadgets")
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))

  # test for errors in plot_new
  expect_error(print(plot_new), NA)

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

  dev.off()
})

