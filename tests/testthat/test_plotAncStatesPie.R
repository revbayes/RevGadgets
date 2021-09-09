context("tests the plotAncStatesPie function")
# note - this test does NOT compare the layers elements of the plot objects.
test_that("plots pies of ancestral states", {
  # get files
  tree_file <-
    system.file("extdata",
                "dec/small_dec.tre",
                package = "RevGadgets")
  plot_file <-
    system.file("extdata",
                "graphs/plotAncStatesPie_df.rds",
                package = "RevGadgets")

  # make a new plot
  # labels that correspond to each region/ possible combination of regions
  labs <- c(
    "1" = "K",
    "2" = "O",
    "3" = "M",
    "4" = "H",
    "5" = "KO",
    "6" = "KM",
    "7" = "OM",
    "8" = "KH",
    "9" = "OH",
    "10" = "MH",
    "11" = "KOM",
    "12" = "KOH",
    "13" = "KMH",
    "14" = "OMH",
    "15" = "KOMH"
  )
  dec_example  <- processAncStates(tree_file , state_labels = labs)
  # Use the state_labels in the returned tidytree object to define color palette
  # These state_labels may be a subset of the labels you provided
  # (not all possible regions may be sampled in the dataset)
  colors <-
    colorRampPalette(colFun(12))(length(dec_example@state_labels))
  names(colors) <- dec_example@state_labels
  # create plot
  plot_new <-
    plotAncStatesPie(
      t = dec_example,
      pie_colors = colors,
      tip_labels_size = 2,
      cladogenetic = TRUE,
      tip_labels_offset = 0.01,
      timeline = F,
      node_pie_size = .5,
      tip_pie_size = .3
    ) +
    ggplot2::scale_x_continuous(limits = c(-0.5, 1)) +
    ggplot2::theme(legend.position = c(0.1, 0.75))

  # read original plot object
  plot_orig <- readRDS(plot_file)

  tmp <- tempdir()
  pdf(paste0(tmp,"/Rplots.pdf"))

  # test for errors in plot_new
  expect_error(print(plot_new), NA)

  dev.off()

  # compare plot data objects
  expect_equal(plot_new$data, plot_orig)

})
