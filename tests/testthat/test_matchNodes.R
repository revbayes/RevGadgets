context("Tests matchNodes function")

test_that("compare node map from documentation example", {
  treefile <-
    system.file("extdata", "bds/primates.tre", package = "RevGadgets")
  tree <- readTrees(treefile)
  map <- matchNodes(tree[[1]][[1]]@phylo)

  expect_equal(class(map), "data.frame")
  expect_equal(dim(map), c(465, 2))
  expect_equal(map[1, 1], 1)
  expect_equal(map[1, 2], 233)
})
