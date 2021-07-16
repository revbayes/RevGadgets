context("tests the dropTip function")

test_that("drops tips", {
  # load in the file
  file <-
    system.file("extdata",
                "sub_models/primates_cytb_GTR_MAP.tre",
                package = "RevGadgets")
  tree <- readTrees(paths = file)
  # root, then drop 1 tip
  tree <- rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")
  tree_dropped <- dropTip(tree, "Otolemur_crassicaudatus")

  # check that tree_rooted is a list of lists of a treedata object
  expect_equal(class(tree_dropped), "list")
  expect_equal(length(tree_dropped), 1)
  expect_equal(class(tree_dropped[[1]]), "list")
  expect_equal(length(tree_dropped[[1]]), 1)
  expect_equal(class(tree_dropped[[1]][[1]])[1], "treedata")

  # check that dropped tip tree has one fewer tip
  expect_equal(length(tree[[1]][[1]]@phylo$tip.label),
               length(tree_dropped[[1]][[1]]@phylo$tip.label)+1)

  # check that data was also dropped
  expect_equal(dim(tree[[1]][[1]]@data)[1],
               dim(tree_dropped[[1]][[1]]@data)[1] + 2)

})
