context("tests the plotAncStatesDiscrete summary_statistic=PieState function")

test_that("plots PieState ancestral states", {

  # get files
  tree_file = system.file("extdata", "discrete_ancestral_states/ancestral_states_pie.tree", package="RevGadgets")
  plot_file = system.file("extdata", "discrete_ancestral_states/ancestral_state_plot_PieState.rds", package="RevGadgets")

  # make a new plot
  require(treeio)
  t = treeio::read.beast(tree_file)
  p = plotAncStatesDiscrete(t,
                            size=0.8,
                            summary_statistic="PieState",
                            include_start_states=FALSE,
                            tip_label_size=1.5,
                            tip_label_offset=0.2,
                            node_size_range=c(2, 2),
                            tip_label_italics=FALSE,
                            node_label_size=0.0,
                            show_posterior_legend=!FALSE,
                            state_labels=c(as.character(0:2)),
                            node_pie_diameter=0.005,
                            pie_nudge_x=0.02,
                            pie_nudge_y=0.2,
                            alpha=.9)
  p = p + ggplot2::coord_cartesian(xlim = c(0, 20))
  #print(p)
  
  # read original plot object
  original_p = readRDS(plot_file)

  # compare plot objects
  expect_equal(p$data, original_p$data)
  expect_equal(p$mapping, original_p$mapping)
  expect_equal(p$theme, original_p$theme)
  expect_equal(p$coordinates, original_p$coordinates)
  expect_equal(p$labels, original_p$labels)
  expect_equal(p$guides, original_p$guides)
  expect_equal(p$layers[[1]], original_p$layers[[1]])
  expect_equal(p$layers[[2]], original_p$layers[[2]])
  expect_equal(p$layers[[4]], original_p$layers[[4]])

})
