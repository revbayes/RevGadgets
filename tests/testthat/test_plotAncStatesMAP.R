context("tests the plotAncStatesDiscrete summary_statistic=MAP function")

test_that("plots MAP ancestral states", {

  # get files
  tree_file = system.file("extdata", "discrete_ancestral_states/ancestral_states.tree", package="RevGadgets")
  plot_file = system.file("extdata", "discrete_ancestral_states/ancestral_state_plot_MAP.rds", package="RevGadgets")

  # make a new plot
  require(treeio)
  t = treeio::read.beast(tree_file)
  p = plotAncStatesDiscrete(t,
                            size=0.8,
                            summary_statistic="MAP",
                            include_start_states=FALSE,
                            tip_label_size=1.5,
                            tip_label_offset=0.2,
                            node_size_range=c(2, 2),
                            tip_label_italics=FALSE,
                            node_label_size=0.0,
                            show_posterior_legend=FALSE,
                            alpha=.7)
  p = p + ggplot2::coord_cartesian(xlim = c(0, 20))
  
  # read original plot object
  original_p = readRDS(plot_file)

  # compare plot objects
  original_p$layers[[3]]$aes_params = p$layers[[3]]$aes_params
  original_p$scales$scales = p$scales$scales
  original_p$plot_env$p$scales$scales = p$plot_env$p$scales$scales
  attributes(original_p$plot_env$tr)$file = attributes(p$plot_env$tr)$file
  expect_equal(p, original_p)

})
