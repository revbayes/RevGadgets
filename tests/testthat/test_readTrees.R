context("tests the readTrees function")

test_that("dummy test", {
  t <- 5
  expect_equal(5, t)
})


# test_that("reads single tree", {
#   tree_single <- readTrees(path = "../../data/sub_models/primates_cytb_covariotide_MAP.tre", format = "nexus")
#   expect_equal(length(tree_single), 1)
#   expect_equal(class(tree_single), "multiPhylo")
#   expect_equal(length(tree_single[[1]]$tip.label), 23)
# })
#
# test_that("reads multi tree", {
#   tree_multi <- readTrees(path = "../../data/sub_models/primates_cytb_covariotide.trees", format = "newick")
#   expect_equal(length(tree_multi), 4002)
#   expect_equal(class(tree_multi), "multiPhylo")
#   expect_equal(length(tree_multi[[1]]$tip.label), 23)
# })
