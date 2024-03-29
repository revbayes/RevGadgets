context("tests the rerootPhylo function")

test_that("reroots tree", {
  # load in the file
  file <-
    system.file("extdata",
                "sub_models/primates_cytb_GTR_MAP.tre",
                package = "RevGadgets")
  tree <- readTrees(paths = file)
  # root with one taxon
  tree_rooted <-
    rerootPhylo(tree = tree, outgroup = "Galeopterus_variegatus")

  # check that tree_rooted is a list of lists of a treedata object
  expect_equal(class(tree_rooted), "list")
  expect_equal(length(tree_rooted), 1)
  expect_equal(class(tree_rooted[[1]]), "list")
  expect_equal(length(tree_rooted[[1]]), 1)
  expect_equal(class(tree_rooted[[1]][[1]])[1], "treedata")

  # if input is unrooted, check that rerooted added a branch to edge matrix
  expect_equal(dim(tree[[1]][[1]]@phylo$edge)[1], 43)
  expect_equal(dim(tree_rooted[[1]][[1]]@phylo$edge)[1], 44)
  
  # check that data are correctly re-associated 
  node <- ape::getMRCA(tree[[1]][[1]]@phylo, c("Callicebus_donacophilus", 
                                               "Pan_paniscus"))
  node_pp <- tree[[1]][[1]]@data[which(tree[[1]][[1]]@data$node == node), 
                                 "posterior"]
  
  node_root <- ape::getMRCA(tree_rooted[[1]][[1]]@phylo, c("Saimiri_sciureus",
                                                           "Cebus_albifrons"))
  
  node_pp_root <- 
    tree_rooted[[1]][[1]]@data[which(tree_rooted[[1]][[1]]@data$node == 
                                       node_root), 
                               "posterior"]
  
  expect_equal(node_pp, node_pp_root)

})
