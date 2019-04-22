context("tests the readAnnTrees function")

test_that("reads single annotated tree", {
  tree_single <- readAnnTrees(path = "../../data/comp_method_disc/ase_freeK.tree")
  expect_equal(length(tree_single), 1)
  expect_equal(class(tree_single), "beastList")
  expect_equal(length(tree_single[[1]]@phylo$tip.label), 342)
})

test_that("reads multi annotated trees", {
  tree_multi <- readAnnTrees(path = "../../data/nexus_multi_ann.nex")
  expect_equal(length(tree_multi), 5)
  expect_equal(class(tree_multi), "beastList")
  expect_equal(length(tree_multi[[1]]@phylo$tip.label), 342)
})
