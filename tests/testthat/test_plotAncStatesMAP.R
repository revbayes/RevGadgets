context("tests the plotAncStatesMAP function")

test_that("plots MAP of ancestral states", {
  # get files
  tree_file <-
    system.file("extdata",
                "comp_method_disc/ase_freeK.tree",
                package = "RevGadgets")
  plot_file <-
    system.file("extdata",
                "graphs/plotAncStatesMAP_df.rds",
                package = "RevGadgets")

  # make a new plot
  example <-
    processAncStates(tree_file,
                     state_labels = c("1" = "Awesome",
                                      "2" = "Beautiful",
                                      "3" = "Cool!"))
  plot_new <- plotAncStatesMAP(t = example)

  # read original plot object
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))

  #  plot new doesn't error out
  expect_error(print(plot_new), NA)

  # compare plot dataobjects
  expect_equal(plot_new$data, plot_orig)

  dev.off()
})

test_that("error messages behave as expected", {
  # get files
  tree_file <-
    system.file("extdata",
                "comp_method_disc/ase_freeK.tree",
                package = "RevGadgets")

  # make a new plot
  example <-
    processAncStates(tree_file,
                     state_labels = c("1" = "Awesome",
                                      "2" = "Beautiful",
                                      "3" = "Cool!"))
  expect_error(plotAncStatesMAP(t = example,
                                geo_units =
                                  list('epochs','periods','years')))

})
