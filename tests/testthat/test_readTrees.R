context("tests the readTrees function")

test_that("reads single tree", {
  tree_single <- readTrees(path = "../../data/mammals_thinned.tree", format = "nexus")
  expect_equal(length(tree_single), 1)
  expect_equal(class(tree_single), "multiPhylo")
  expect_equal(length(tree_single[[1]]$tip.label), 342)
})

test_that("reads multi tree", {
  tree_multi <- readTrees(path = "../../data/trees_posterior_1.trees", format = "newick")
  expect_equal(length(tree_multi), 86)
  expect_equal(class(tree_multi), "multiPhylo")
  expect_equal(length(tree_multi[[1]]$tip.label), 49)
})
