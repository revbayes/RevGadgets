context("tests the plotAncStatesDiscrete summary_statistic=MAPChromosome function")

test_that("plots chromosome ancestral states", {

  # get files
  tree_file = system.file("extdata", "chromo/chromosomes_ancestral_states.tree", package="RevGadgets")
  plot_file = system.file("extdata", "chromo/chromosomes_ancestral_state_plot.rds", package="RevGadgets")

  # make a new plot
  require(treeio)
  t = treeio::read.beast(tree_file)
  p = plotAncStatesDiscrete(t,
                            summary_statistic="MAPChromosome",
                            include_start_states=TRUE,
                            tip_label_offset=0.4,
                            node_label_size=2,
                            shoulder_label_size=2,
                            tip_label_size=3)
  p = p + ggplot2::coord_cartesian(xlim = c(0, 16.5))

  # read original plot object
  original_p = readRDS(plot_file)

  # compare plot objects
  attributes(p$plot_env$tr)$file = ""
  attributes(original_p$plot_env$tr)$file = ""
  expect_equal(isTRUE(all.equal(p, original_p)), TRUE)

})
