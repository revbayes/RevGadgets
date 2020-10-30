context("tests the plotAncStatesPie function")

test_that("plots pies of ancestral states", {

  # get files
  tree_file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
  plot_file <- system.file("extdata", "graphs/plotAncStatesPie.RData", package="RevGadgets")

  # make a new plot
  # labels that correspond to each region/ possible combination of regions
  labs <- c("1" = "K", "2" = "O", "3" = "M", "4" = "H", "5" = "KO",
            "6" = "KM", "7" = "OM", "8" = "KH", "9" = "OH", "10" = "MH",
            "11" = "KOM", "12" = "KOH", "13" = "KMH", "14" = "OMH", "15" = "KOMH")
  dec_example  <- processAncStates(tree_file , state_labels = labs)
  # Use the state_labels in the returned tidytree object to define color palette
  # These state_labels may be a subset of the labels you provided
  # (not all possible regions may be sampled in the dataset)
  colors <- colorRampPalette(RevGadgets:::.colFun(12))(length(dec_example@state_labels))
  names(colors) <- dec_example@state_labels
  # plot
  plot_new <- plotAncStatesPie(t = dec_example, pie_colors = colors, tip_labels_size = 3,
                        cladogenetic = TRUE, tip_labels_offset = 0.25, timeline = T) +
    ggplot2::theme(legend.position = c(0.1, 0.75))


  # read original plot object
  load(plot_file)

  # compare plot objects
  expect_equal(plot[[1]], plot_new[[1]])

})
