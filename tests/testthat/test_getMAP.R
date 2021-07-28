context("Tests getMAP() function for traces")

test_that("compare getMAP() for example from documentation", {
  # read in and process example data
  file <- system.file("extdata",
                      "sub_models/primates_cytb_GTR_mini.p",
                      package = "RevGadgets")
  trace <- readTrace(paths = file)

  #some error is to be expected given the function, so round to 3 digits
  expect_equal(0.26, round(getMAP(trace[[1]]$"pi[1]"), digits = 2))
})
