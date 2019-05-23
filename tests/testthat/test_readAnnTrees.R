context("tests the readAnnTrees function")

test_that("reads single annotated tree", {
  file <- system.file("extdata", "comp_method_disc/ase_freeK.tree", package="RevGadgets")
  tree_single <- readAnnTrees(path = file)
  expect_equal(length(tree_single), 1)
})


test_that("reads multi annotated trees", {
  file <- system.file("extdata", "nexus_multi_ann.nex", package="RevGadgets")
  tree_multi <- readAnnTrees(path = file)
  expect_equal(length(tree_multi), 5)
  expect_equal(class(tree_multi), "beastList")
  expect_equal(length(tree_multi[[1]]@phylo$tip.label), 342)
})
