context("Tests processBranchData function")

test_that("compare annotated tree from documentation example", {
  treefile <-
    system.file("extdata", "bds/primates.tre", package = "RevGadgets")
  logfile <-
    system.file("extdata", "bds/primates_BDS_rates_truncated.p", package =
                  "RevGadgets")

  branch_data <- readTrace(logfile)[[1]]
  tree <- readTrees(paths = treefile)[[1]][[1]]

  annotated_tree <-
    processBranchData(tree, branch_data, summary = "median")

  expect_equal(class(annotated_tree), "list")
  expect_equal(length(annotated_tree), 1)
  expect_equal(class(annotated_tree[[1]]), "list")
  expect_equal(length(annotated_tree[[1]]), 1)
  expect_equal(class(annotated_tree[[1]][[1]])[1], "treedata")
  expect_equal(class(annotated_tree[[1]][[1]]@phylo), "phylo")

})
