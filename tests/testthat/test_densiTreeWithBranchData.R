context("tests densiTree-style plots with branch data")

test_that("plot doesn't error out", {
  trees = lapply(1:5, function(x) ape::rcoal(5))
  data = lapply(1:5, function(x) stats::runif(9, 1, 10))
  #TODO file used here is Beast2 file because RB doesn't output distributions as Nexus
  tree_file = system.file("extdata", "beast2", "msbd.rates.trees", package = "RevGadgets")
  
  expect_silent(densiTreeWithBranchData(trees = trees, data = data))
  expect_silent(densiTreeWithBranchData(trees = trees, data = data, data_intervals = c(0,11)))
  expect_silent(densiTreeWithBranchData(trees = trees, data = data, data_intervals = 0:5))
  expect_silent(densiTreeWithBranchData(tree_file = tree_file, data_name = "lambda"))
})

test_that("invalid inputs are rejected", {
  trees = lapply(1:5, function(x) ape::rcoal(5))
  data = lapply(1:5, function(x) stats::runif(9, 1, 10))
  #TODO file used here is Beast2 file because RB doesn't output distributions as Nexus
  tree_file = system.file("extdata", "beast2", "msbd.rates.trees", package = "RevGadgets")
  
  expect_error(densiTreeWithBranchData(trees = trees))
  expect_error(densiTreeWithBranchData(data = data))
  expect_error(densiTreeWithBranchData(tree_file = tree_file, data_name = "psi"))
})