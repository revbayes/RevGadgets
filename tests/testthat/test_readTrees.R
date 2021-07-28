context("tests the readTrees function")

test_that("reads single nexus tree", {
  file <-
    system.file("extdata",
                "sub_models/primates_cytb_GTR_MAP.tre",
                package = "RevGadgets")
  tree_single <- readTrees(paths = file)
  expect_equal(length(tree_single), 1)
  expect_equal(class(tree_single[[1]][[1]])[1], "treedata")
  expect_equal(length(tree_single[[1]][[1]]@phylo$tip.label), 23)
})

test_that("reads tree trace", {
  file <-
    system.file("extdata", "sub_models/primates_cytb_GTR_mini.trees", package =
                  "RevGadgets")
  tree_multi <- readTrees(path = file)
  expect_equal(length(tree_multi[[1]]), 10)
  expect_equal(class(tree_multi[[1]]), "list")
  expect_equal(length(tree_multi[[1]][[1]]@phylo$tip.label), 23)
})

test_that("reads single newick", {
  file <-
    system.file("extdata", "bds/primates.tre", package = "RevGadgets")
  tree_new <- readTrees(path = file)
  expect_equal(length(tree_new[[1]]), 1)
  expect_equal(class(tree_new[[1]][[1]])[1], "treedata")
  expect_equal(length(tree_new[[1]][[1]]@phylo$tip.label), 233)
})
