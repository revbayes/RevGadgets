context("tests the processAncStates function")

test_that("processes ancestral states scripts", {
  # read in and process file
  file <-
    system.file("extdata",
                "comp_method_disc/ase_freeK.tree",
                package = "RevGadgets")
  example <- processAncStates(file,
                              state_labels = c("1" = "Awesome",
                                               "2" = "Beautiful",
                                               "3" = "Cool!"))
  # test file format
  expect_equal(class(example)[1], "treedata")

  # check tree
  expect_equal(example@phylo$Nnode, 341)
  expect_equal(example@phylo$tip.label[1], "Vulpes_macrotis")

  # check data
  expect_equal(class(example@data), c("tbl_df", "tbl", "data.frame"))
  expect_equal(names(table(example@data$anc_state_1)),
               c("Awesome", "Beautiful", "Cool!"))

})
